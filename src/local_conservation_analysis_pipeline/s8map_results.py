import json
from pathlib import Path

import pandas as pd
from openpyxl import Workbook, load_workbook

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools


def get_failure_map(json_files):
    failure_map = {}
    for json_file in json_files:
        with open(json_file, "r") as f:
            json_dict = json.load(f)
        if "critical_error" in json_dict:
            failure_map[json_dict['reference_index']] = json_dict["critical_error"]
    return failure_map


def get_json_map(json_files):
    json_map = {}
    for json_file in json_files:
        with open(json_file, "r") as f:
            json_dict = json.load(f)
        json_map[json_dict["reference_index"]] = json_file
    return json_map


def check_json(json_file):
    with open(json_file, "r") as f:
        json_dict = json.load(f)
    if "reference_index" in json_dict:
        return True
    else:
        return False


def get_image_map(json_files, image_score_key):
    image_map = {}
    for json_file in json_files:
        with open(json_file, "r") as f:
            json_dict = json.load(f)
        if f"multilevel_plot_file-{image_score_key}" in json_dict:
            # file = Path(json_dict["multilevel_plot_file-property_entropy"]).resolve().relative_to(Path.cwd())
            file = Path(json_dict[f"multilevel_plot_file-{image_score_key}"]).resolve().relative_to(Path.cwd())
            image_map[json_dict["reference_index"]] = rf'=HYPERLINK("./{file}")'
    return image_map


def find_motif_regex(og: group_tools.ConserGene, regex):
    hit_sequence = og.hit_sequence
    matches = list(tools.get_regex_matches(regex, hit_sequence))
    if len(matches) == 0:
        return None, None, None
    if len(matches) > 1:
        return None, None, None
    m = matches[0]
    matchseq = m[0]
    matchst = m[1]
    matchend = m[2]
    return matchseq, matchst, matchend


def main(search_dir, table_file, image_score_key):
    json_files = Path(search_dir).rglob("*.json")
    checked_jsons = [i for i in json_files if check_json(i)]
    table_file = table_file.replace(".csv", "_original_reindexed.csv")
    table_df = pd.read_csv(table_file)
    failure_map = get_failure_map(checked_jsons)
    json_map = get_json_map(checked_jsons)
    table_df["fail_reason"] = table_df["reference_index"].map(failure_map)
    table_df["json_file"] = table_df["reference_index"].map(json_map)
    image_map = get_image_map(checked_jsons, image_score_key)
    table_df["image_file"] = table_df["reference_index"].map(image_map)







    table_df.to_csv(table_file, index=False)











    # from openpyxl import Workbook
    # table_df.to_excel(table_file.replace(".csv", ".xlsx"), index=False)
    # workbook = Workbook()
    # sheet = workbook.active
    # for index, row in table_df.iterrows():
    #     sheet.append(row.to_dict())
    # workbook.save('your_excel_file.xlsx')
