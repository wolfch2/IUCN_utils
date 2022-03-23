import ctypes
import functools
import glob
import json
import math
import multiprocessing
import operator
import os
import random
import shutil
import subprocess  # https://docs.python.org/3/library/subprocess.html#replacing-os-system
import time
from collections import defaultdict
from pathlib import (
    Path,  # https://stackoverflow.com/questions/273192/how-can-i-safely-create-a-nested-directory-in-python
)

import affine
import ee
import fiona
import geojson
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import ray
import rioxarray as rxr  # for the extension to load
import shapely
from dbfread import DBF, FieldParser, InvalidValue
from fiona import transform
from osgeo import gdal, ogr
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
from pyproj import Proj, Transformer, transform
from rasterio import features
from rasterio.crs import CRS
from rasterio.features import bounds as calculate_bounds
from shapely.geometry import shape


class MyFieldParser(FieldParser):
    def parse(self, field, data):
        try:
            return FieldParser.parse(self, field, data)
        except ValueError:
            return InvalidValue(data)


@ray.remote
class MyActor(object):
    def __init__(
        self, elev_mat_id, project_dir, cache_dir, proj
    ):  # TODO -- maybe only pass in one directory lol
        self.elev = rasterio.open(project_dir + "/elev_proj.tif")
        self.elev_mat = elev_mat_id
        self.richness = (0 * self.elev_mat).astype("int16")
        self.cache_dir = cache_dir  # TODO - maybe reorder constructor lines for clarity
        self.proj = proj

    def add_range(self, id_info):
        print(id_info["id"])
        rast_file_name = self.cache_dir + "/" + str(id_info["id"]) + ".npz"
        if os.path.exists(rast_file_name):  # TODO - make caching optional??
            print("using cached raster")
            # TODO - Consider having each actor maintain cached ranges in a separate npz file (create on init)
            # to avoid cluttering the cache folder with many little files. Not sure how this would affect reading
            # though (can look into it)
            # TODO - Consider looking into HDF5 as a replacement option (with parallel writing):
            # https://stackoverflow.com/questions/48603119/incremently-appending-numpy-arrays-to-a-save-file
            range_rast_packed = np.load(rast_file_name)[
                "arr_0"
            ]  # see np.load(rast_file_name).files
            range_rast = np.unpackbits(range_rast_packed)[
                : np.prod(self.elev.shape)
            ].reshape(self.elev.shape)
        else:
            print("rasterizing")
            all_features = fiona.open(id_info["path"])
            shapes = [
                shape(all_features[index]["geometry"]) for index in id_info["rows"]
            ]
            gdf = gpd.GeoDataFrame(
                index=range(len(shapes)), crs=all_features.crs["init"], geometry=shapes
            ).to_crs(self.proj)
            ###
            range_rast = rasterio.features.rasterize(
                gdf["geometry"],
                out_shape=self.elev.shape,
                transform=self.elev.transform,
                default_value=1,
                all_touched=False,
            )
            range_rast[
                np.logical_or(
                    self.elev_mat < id_info["lower"], self.elev_mat > id_info["upper"]
                )
            ] = 0
        ### save range if we are caching; TODO -- look into more efficient format (scipy sparse array?)
        if not os.path.exists(self.cache_dir + "/" + str(id_info["id"]) + ".tif"):
            np.savez_compressed(rast_file_name, np.packbits(range_rast))
        ###
        self.richness += range_rast
        return

    def get_richness(self):
        return self.richness


class SpatialTools:
    def __init__(self, dl_dir, cred_file):
        self.dl_dir = dl_dir  # TODO - remove trailing slash if present
        self.cred_file = cred_file
        return

    def _dl_elevation(self, proj, res):
        if os.path.exists(self.dl_dir + "/elev.tif"):
            return
        ee.Initialize()
        Path(self.dl_dir).mkdir(parents=True, exist_ok=True)
        scale = 0.008333333333333  # https://lta.cr.usgs.gov/sites/default/files/GTOPO30_README_0.doc
        ee.batch.Export.image.toDrive(
            ee.Image("USGS/GTOPO30"),
            description="elev",
            crsTransform=[scale, 0, -180, 0, -scale, 90],
            crs="EPSG:4326",
            dimensions=str(int(360 / scale)) + "x" + str(int(180 / scale)),
            folder="GEE_rasts",
            maxPixels=1e11,
        ).start()
        while True:  # wait for task to complete
            states = pd.DataFrame(map(lambda x: x.status(), ee.batch.Task.list()))[
                "state"
            ]
            if "RUNNING" in states.unique() or "READY" in states.unique():
                print(
                    "Waiting on GEE.  Tasks: "
                    + str(((states == "READY") | (states == "RUNNING")).sum())
                )
                time.sleep(60)
            else:
                break
        ### dl from Google Drive
        gauth = GoogleAuth()
        gauth.LoadCredentialsFile(self.cred_file)
        if gauth.access_token_expired:
            gauth.Refresh()
        drive = GoogleDrive(gauth)
        root_files = drive.ListFile(
            {"q": "'root' in parents and trashed=false"}
        ).GetList()
        root_files_df = pd.DataFrame(root_files)
        folder_id = root_files_df["id"][root_files_df["title"] == "GEE_rasts"].tolist()[
            0
        ]
        files = drive.ListFile(
            {"q": "'{}' in parents and trashed=false".format(folder_id)}
        ).GetList()
        for file in files:
            file.GetContentFile(self.dl_dir + "/elev.tif")
            file.Delete()
        ### reproject
        rast = rxr.open_rasterio(
            self.dl_dir + "/elev.tif"
        ).squeeze()  # rast.rio.resolution() ; rast.rio.bounds()
        rast_proj = rast.rio.reproject(
            dst_crs=CRS.from_string(proj), resolution=res
        )  # to numpy w/ .values
        rast_proj.rio.to_raster(
            self.dl_dir + "/elev_proj.tif", compress="LZW", tiled=True
        )
        return

    def _dl_ranges(self):
        # if os.path.exists(self.dl_dir + 'ranges_vector/MAMMALS_TERRESTRIAL_ONLY.shp'): # TODO - maybe make this work
        #    return
        range_dir = self.dl_dir + "/ranges_vector"
        Path(range_dir).mkdir(parents=True, exist_ok=True)
        subprocess.call(
            f'wget -nc -P "{range_dir}" "https://sdownloads.iucnredlist.org/groups/MAMMALS_TERRESTRIAL_ONLY.zip"',
            shell=True,
        )
        subprocess.call(
            f'unzip -n -d "{range_dir}" "{range_dir}/MAMMALS_TERRESTRIAL_ONLY.zip"',
            shell=True,
        )
        return

    def get_indices(self, project_dir, species_data):  # TODO - underscore prefix
        keep_ids = species_data.taxonid.tolist()
        species_data["lower"] = species_data["elevation_lower"]
        species_data["upper"] = species_data["elevation_upper"]
        ###
        species_data.loc[pd.isna(species_data["lower"]), "lower"] = -1e7
        species_data.loc[pd.isna(species_data["upper"]), "upper"] = 1e7
        ids = []
        for gp in [
            "MAMMALS_TERRESTRIAL_ONLY"
        ]:  # ['MAMMALS','AMPHIBIANS','REPTILES'] + ['BIRDS_' + str(x) for x in range(1,9)]:
            print(gp)
            shp_path = project_dir + "/data/IUCN_spatial/ranges_vector/" + gp + ".shp"
            dbf_path = project_dir + "/data/IUCN_spatial/ranges_vector/" + gp + ".dbf"
            gp_data = pd.DataFrame(
                iter(
                    DBF(
                        dbf_path,
                        encoding="utf-8",
                        char_decode_errors="ignore",
                        parserclass=MyFieldParser,
                    )
                )
            )
            if "BIRDS_" in gp:
                gp_data["presence"] = gp_data["presenc"]
            gp_ids = (
                gp_data.loc[(gp_data["presence"] <= 2) & (gp_data["origin"] <= 2)][
                    "id_no"
                ]
                .unique()
                .tolist()
            )
            gp_ids = list(set(keep_ids) & set(gp_ids))
            for gp_id in gp_ids:
                ids += [
                    {
                        "path": shp_path,
                        "id": int(gp_id),
                        "rows": gp_data.index[gp_data["id_no"] == gp_id].tolist(),
                        "lower": float(
                            species_data.loc[species_data["taxonid"] == gp_id]["lower"]
                        ),
                        "upper": float(
                            species_data.loc[species_data["taxonid"] == gp_id]["upper"]
                        ),
                    }
                ]
        print("Species matched:" + str(len(ids)))
        return ids

    def calculate_richness(self, project_dir, species_df, rast_name, proj, num_cpu):
        if os.path.exists(project_dir + "/" + rast_name + ".tif"):
            return
        Path(project_dir + "/ranges_cache").mkdir(
            parents=True, exist_ok=True
        )  # TODO - only run if caching ranges...
        id_info_list = self.get_indices(project_dir, species_df)
        print(id_info_list)
        elev = rasterio.open(project_dir + "/data//elev_proj.tif")
        elev_mat = elev.read(1)
        #
        ray.init(
            object_store_memory=25 * 1024 * 1024 * 1024
        )  # TODO - change back to 45
        elev_mat_id = ray.put(elev_mat)
        actors = [
            MyActor.remote(
                elev_mat_id, project_dir, project_dir + "/ranges_cache", proj
            )
            for _ in range(num_cpu)
        ]
        pool = ray.util.ActorPool(actors)
        list(
            pool.map_unordered(lambda actor, v: actor.add_range.remote(v), id_info_list)
        )
        #
        results = ray.get([actor.get_richness.remote() for actor in actors])
        ray.shutdown()
        richness = sum(results)
        #
        rast = rasterio.open(
            project_dir + "/" + rast_name + ".tif",
            "w",
            driver="GTiff",
            height=elev.shape[0],
            width=elev.shape[1],
            count=1,
            dtype=rasterio.int16,
            crs=proj,
            compress="lzw",
            transform=elev.transform,
        )
        rast.write(richness, 1)
        rast.close()
        return

    # https://stackoverflow.com/questions/5615648/how-to-call-a-function-within-class
    def run(self, species_df, proj, res, group_name):
        print("Downloading elevation and species range data")
        self._dl_elevation(proj, res)
        self._dl_ranges()
        print("Rasterizing and calculating richness")
        self.calculate_richness(self.dl_dir, species_df, group_name, proj, num_cpu=3)
        return


if __name__ == "__main__":
    species_df = pd.read_excel(
        "/media/data_disk_01/shared/analysis/misc/carbon_diversity/data/species_data/species_data.xlsx"
    )
    ###
    spatial_tools = SpatialTools(
        dl_dir="/mnt/shared/analysis/misc/carbon_diversity/data/IUCN_spatial",
        cred_file="/mnt/shared/analysis/mycreds.txt",
    )
    spatial_tools.run(
        species_df=species_df,
        proj="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs",
        res=1000,
        group_name="threatened forest primates",
    )
