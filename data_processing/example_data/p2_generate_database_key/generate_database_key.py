import json
from pathlib import Path

database_dir = Path("../../../data/example_orthogroup_database/human_odb_groups")
JSON_DIR = database_dir / "info_jsons"
OUTPUT_FILE = database_dir / "database_key.json"


def main(json_dir, output_file):
    json_files = json_dir.glob("*.json")
    database_key = {}
    for json_file in json_files:
        with open(json_file, "r") as f:
            json_dict = json.load(f)
        odb_gene_id = json_dict["query_odb_gene_id"]
        level = json_dict["oglevel"]
        if odb_gene_id not in database_key:
            database_key[odb_gene_id] = {}
        database_key[odb_gene_id]["query_uniprot_id"] = json_dict["query_uniprot_id"]
        if level in database_key[odb_gene_id]:
            raise ValueError(f"level {level} already in database_key for {odb_gene_id}")
        database_key[odb_gene_id][level] = {}
        lvl_dict = database_key[odb_gene_id][level]
        lvl_dict["alignment_file"] = json_dict["alignment_clustered_ldos_file"]
        # lvl_dict["conservation_scores"] = json_dict["conservation_scores"]
    with open(output_file, "w") as f:
        json.dump(database_key, f, indent=4)


if __name__ == "__main__":
    main(JSON_DIR, OUTPUT_FILE)
