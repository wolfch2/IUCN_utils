import glob
import json
import os
import shutil
import subprocess  # https://docs.python.org/3/library/subprocess.html#replacing-os-system
from pathlib import (
    Path,  # https://stackoverflow.com/questions/273192/how-can-i-safely-create-a-nested-directory-in-python
)

import pandas as pd


class Downloader:
    def __init__(self, api_key_path, dl_dir, class_list):
        with open(api_key_path) as f:
            self.api_key = f.readlines()[0].replace("\n", "")
        self.dl_dir = dl_dir  # TODO - remove trailing slash if present
        self.class_list = class_list
        return

    def _dl_base(self):
        Path(self.dl_dir + "/pages").mkdir(parents=True, exist_ok=True)
        page_df_list = []
        page = 0
        while True:
            print(f"Downloading page {page}")
            subprocess.call(
                f'wget -nc -P "{self.dl_dir}/pages" "http://apiv3.iucnredlist.org/api/v3/species/page/{page}?token={self.api_key}"',
                shell=True,
            )
            with open(f"{self.dl_dir}/pages/{page}?token={self.api_key}") as f:
                page_json = json.load(f)
            if len(page_json["result"]) == 0:
                break
            page_df_list += [pd.DataFrame(page_json["result"])]
            page += 1
        self.base_df = pd.concat(page_df_list, axis=0).query(
            "infra_rank.isna()", engine="python"
        )  # remove subspecies and population
        if self.class_list != None:
            self.base_df = self.base_df[self.base_df.class_name.isin(self.class_list)]
        return

    # download a particular type of data (e.g., habitats, countries, narratives)
    # "main" is just the base factsheet page (have to adjust url accordingly)
    def _dl_type(self, type_name):
        if os.path.exists(f"{self.dl_dir}/{type_name}.json"):
            print("Already complete")
            return
        Path(self.dl_dir + "/" + type_name).mkdir(parents=True, exist_ok=True)
        for (
            id
        ) in (
            self.base_df.taxonid
        ):  # TODO - skip already downloaded files listed in the .json
            url = f"http://apiv3.iucnredlist.org/api/v3/{type_name}/species/id/{id}?token={self.api_key}"
            if type_name == "main":
                url = url.replace("main/", "")  # this is just the base url
            subprocess.call(
                f'wget -nc --tries=0 -P "{self.dl_dir}/{type_name}" "{url}"', shell=True
            )
        # TODO - verify that we downloaded all files successfully
        json_list = []
        data_files = glob.glob(f"{self.dl_dir}/{type_name}/*")
        for data_file in data_files:
            with open(data_file) as f:
                json_list += [json.load(f)]
        with open(f"{self.dl_dir}/{type_name}.json", "w") as f:
            json.dump(json_list, f)
        [os.remove(x) for x in data_files]
        return

    # TODO -- probably switch to just deleting files that have the api key in their names (seems safer)
    def clear_cache(self, prompt=True):
        if prompt:
            answer = input(f"Clear all files in '{self.dl_dir}' (y/n)? ")
            if answer != "y":
                return
        if os.path.exists(self.dl_dir):
            shutil.rmtree(self.dl_dir)
        return

    # build full dataset
    def _combine_info(self):
        ### habitat codes
        with open(f"{self.dl_dir}/habitats.json") as f:
            habitats_json = json.load(f)  # ignore suitability, importance, etc. for now
        codes = pd.DataFrame(
            [
                {
                    "taxonid": x["id"],
                    "codes": " ".join([y["code"] for y in x["result"]]),
                }
                for x in habitats_json
            ]
        )
        codes.taxonid = codes.taxonid.astype("int64")
        ### elevations
        with open(f"{self.dl_dir}/main.json") as f:
            main_json = json.load(f)
        other = pd.DataFrame([x["result"][0] for x in main_json])
        ### join all
        return codes.merge(other, on="taxonid")

    # https://stackoverflow.com/questions/5615648/how-to-call-a-function-within-class
    def run(self):
        print("Downloading species data")
        self._dl_base()
        ###
        self._dl_type("main")
        self._dl_type("habitats")
        ###
        return self._combine_info()


if __name__ == "__main__":
    downloader = Downloader(
        api_key_path="/media/data_disk_01/shared/analysis/misc/carbon_diversity/api_key.txt",
        dl_dir="/media/data_disk_01/shared/analysis/misc/carbon_diversity/data/species_data",
        class_list=["MAMMALIA"],
    )
    species_df = downloader.run()
