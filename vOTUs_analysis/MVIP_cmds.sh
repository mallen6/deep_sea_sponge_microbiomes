conda install mvip

source activate mvip

cd /path/to/folder/containing/metagenomes/

mvip MVP_00_set_up_MVP -i /MVP_results -m metadata.txt

mvip MVP_01_run_genomad_checkv -i /MVP_results -m metadata.txt --threads 36

mvip MVP_02_filter_genomad_checkv -i /MVP_results -m metadata.txt

mvip MVP_03_do_clustering -i /MVP_results -m metadata.txt --read-type long --threads 24

mvip MVP_04_do_read_mapping -i /MVP_results -m metadata.txt --read_type long --threads 24 --delete_files
