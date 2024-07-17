```bash
mamba activate odb_groups_x86
cd ./p1_initial_database_preparation/
bash create_database.sh
cd ..

cd ./p2_alignment_conservation_scores/
mamba activate slim_conservation_scoring
python property_entropy_scores.py
cd ..
cd ./p3_generate_database_key
python ./generate_database_key.py 
```