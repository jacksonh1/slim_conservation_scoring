```bash
mamba activate slim_conservation_orthogroup_generation
cd ./p1_initial_database_preparation/
bash create_database.sh
```
```bash
cd ..
cd ./p2_generate_database_key
python ./generate_database_key.py 
```

<!-- cd ./p2_alignment_conservation_scores/
mamba activate slim_conservation_scoring
python property_entropy_scores.py
cd .. -->