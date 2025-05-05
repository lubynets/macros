import argparse
import yaml
import pandas as pd
import numpy as np
import uproot
from hipe4ml.model_handler import ModelHandler
import sys
sys.path.append('utils')
import utils as utils

#------------ Define config files ------------
parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--config-file-sel', dest='config_file_selections', help='path to the YAML file with selections.', default='')
parser.add_argument('--input-file', dest='input_file', help='input file in parquet.gzip format.', default='')
parser.add_argument('--tree-name', dest='tree_name', help='tree name inside input file if it is provided in .root format.', default='')
parser.add_argument('--model-file', dest='model_file', help='BDT model file in .pkl format.', default='')
parser.add_argument('--output-directory', dest='output_directory', help='destination directory for pandas file with candidates with assigned bdt score values.', default='')
parser.add_argument('--pT-interval', nargs=2, type=float, metavar=('pT_min', 'pT_max'), help='Specify the pT interval as two floats: min max')
args = parser.parse_args()
## Check if config files are provided
if args.config_file_selections == '':
    print('** No selections config file provided. Exiting. **')
    exit()

output_file_name='appliedBdt'

#------------ Open config files ------------
input_file = args.input_file
tree_name = args.tree_name
model_file = args.model_file
output_directory = args.output_directory
config_file_selections = open(args.config_file_selections, 'r')
config_selections = yaml.full_load(config_file_selections)
pT_min, pT_max = args.pT_interval
## Read config file

#------------ Define selections ------------
selections = config_selections['selection']
selections_string = utils.convert_sel_to_string(selections)
if selections_string != '':
    selections_string = f'((fPt >= {pT_min}) & (fPt < {pT_max})) & ' + selections_string
else:
    selections_string = f'((fPt >= {pT_min}) & (fPt < {pT_max}))'

#------------ Apply BDT to data ------------
## Load model
model_hdl = ModelHandler()
model_hdl.load_model_handler(filename=model_file)

## Define input files
BDT_variables = ['fKFChi2PrimProton',
                 'fKFChi2PrimKaon',
                 'fKFChi2PrimPion',
                 'fKFChi2Geo',
                 'fKFChi2Topo',
                 'fKFDecayLengthNormalised',
                 'fLiteNSigTpcPr',
                 'fLiteNSigTpcKa',
                 'fLiteNSigTpcPi',
                 'fKFT',
                 'fKFPt',
                 'fKFMassInv',
                 'fKFSigBgStatus']

## Loop over input files
applied_dfs = []
## Read pandas DataFrame

with uproot.open(f"{input_file}:{tree_name}") as tree:
    df = tree.arrays(BDT_variables, library="pd")

## Remove candidates with infinite and NaN values
if df.isna().any().any():
    print(f"Columns with NaN values in {input_file}:", df.columns[df.isna().any()])
    df.dropna(inplace=True)
elif np.isinf(df).any().any():
    print(f"Columns with infinities in {input_file}:", df.columns[np.isinf(df).any()])
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(inplace=True)

if selections:  # if selections dictionary is not empty
    for index, (variable, selection) in enumerate(selections.items()): # variable - key, selection - value, index - natural numbers from 0 (Cf. std::iota)
        df.query(selection, inplace=True)   # applies selection on it

df.query(f'((fKFPt >= {pT_min}) & (fKFPt < {pT_max}))', inplace=True)

## Calculate model prediction
prediction = model_hdl.predict(df, output_margin=False)
column_names =['bkg_score', 'prompt_score', 'non_prompt_score']
for i_class in range(3):
    df[column_names[i_class]] = prediction[:, i_class]

applied_df = df[['fKFSigBgStatus', 'fLiteNSigTpcPr', 'fLiteNSigTpcKa', 'fLiteNSigTpcPi', 'fKFT', 'fKFPt', 'fKFMassInv', 'bkg_score', 'prompt_score', 'non_prompt_score']]
applied_dfs.append(applied_df)

## Free memory
del df, applied_df

# Concatenate all dataframes
applied_dataframe = pd.concat(applied_dfs, ignore_index=True)
print(f'Number of Lc candidates in data: {len(applied_dataframe)}')
del applied_dfs

data_dict = {col: np.array(applied_dataframe[col]) for col in applied_dataframe.columns}

#------------ Save output ------------
print(f'Saving output file...')

output_file = f"{output_directory}/{output_file_name}.root"
tree_name = "plainTree"  # Name of the TTree in the ROOT file

with uproot.recreate(output_file) as f:
    f[tree_name] = data_dict

print(f'Saving output file: Done.')
