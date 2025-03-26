import time
import argparse
import yaml
import uproot
import pandas as pd

import sys
sys.path.append('utils')
import utils as utils

start_time = time.time()

############# define config file #############
parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', help='path to the YAML file with configuration.', default='')
args = parser.parse_args()

if args.config_file == '':
    print('** No config file provided. Exiting. **')
    exit()

############# Open config files #############
config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)
input_file_names = config['input_files']
input_tree_name = config['input_tree']
output_path = config['output_path']
output_file_key = config['output_file_key']

############# Read trees #############
print(f'')
print(f'############# MC #############')
print('Retreiving candidate trees from ROOT files...', end='\r')
dataframes = []
for file in input_file_names:
  root_file = uproot.open(file)
  tree = root_file.get(input_tree_name)
  df = tree.arrays(library="pd")
  dataframes.append(df)
  del tree, df
  del root_file
print('Retreiving candidate trees from ROOT files: Done')

############# Merge all found dataframes #############
if dataframes:
    print(f'Number of trees found in MC files:', len(dataframes))
    df_merged = pd.concat(dataframes, ignore_index=True)
    print(f'Length of merged MC dataframe:', len(df_merged))
    print(f'Saving parquet-file to {output_path}...', end='\r')
    df_merged.to_parquet(f'{output_path}/{output_file_key}.parquet.gzip', compression='gzip')
    print(f'Saving parquet-file to {output_path}: Done')
else:
    print(f'No trees named {input_tree_name} found in MC files')

print(f'')
print(f'Program completed.')
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total runtime: {elapsed_time / 60:.2f} minutes")
print(f'')
