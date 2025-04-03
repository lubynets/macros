import argparse
import yaml
import pandas as pd
import numpy as np
from hipe4ml.model_handler import ModelHandler
import sys
sys.path.append('utils')
import utils as utils

#------------ Define config files ------------
parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', help='path to the YAML file with MC configuration.', default='')
parser.add_argument('--config-file-sel', dest='config_file_selections', help='path to the YAML file with selections.', default='')
args = parser.parse_args()
## Check if config files are provided
if args.config_file == '':
    print('** No config file provided. Exiting. **')
    exit()
if args.config_file_selections == '':
    print('** No selections config file provided. Exiting. **')
    exit()

#------------ Open config files ------------
config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)
config_file_selections = open(args.config_file_selections, 'r')
config_selections = yaml.full_load(config_file_selections)
## Read config file
input_files = config['input_files']
if isinstance(input_files, str):
    input_files = [input_files]
model_path = config['model_path']
output_directory = config['output_directory']
output_file_name = config['output_file_name']
pT_interval = config['pT_interval']

#------------ Define selections ------------
selections = config_selections['selection']
selections_string = utils.convert_sel_to_string(selections)
if selections_string != '':
    selections_string = f'((fPt >= {pT_interval[0]}) & (fPt < {pT_interval[1]})) & ' + selections_string
else:
    selections_string = f'((fPt >= {pT_interval[0]}) & (fPt < {pT_interval[1]}))'

#------------ Apply BDT to data ------------
## Load model
model_hdl = ModelHandler()
model_hdl.load_model_handler(filename=model_path)

## Define input files
BDT_variables = ['fCandidateSelFlag', 'fSign',
                 'fP', 'fPt', 'fY', 'fEta', 'fPhi', 'fCt',
                 'fPtXi', 'fPtPi0', 'fPtPi1', 'fPtSumPi0Pi1',
                 'fM', 'fInvMassXi', 'fInvMassXiPi0', 'fInvMassXiPi1',
                 'fCpa', 'fCpaXY', 'fCpaXi', 'fCpaXYXi', 'fCpaLam', 'fCpaXYLam', 'fCpaLamToXi', 'fCpaXYLamToXi',
                 'fDecayLength', 'fDecayLengthNormalised', 'fDecayLengthXY', 'fDecayLengthXYNormalised',
                 'fImpactParameterXi', 'fImpactParameterNormalisedXi', 'fImpactParameterPi0', 'fImpactParameterNormalisedPi0', 'fImpactParameterPi1', 'fImpactParameterNormalisedPi1', 'fMaxNormalisedDeltaIP',
                 'fChi2Sv', 'fChi2XiVtx', 'fChi2LamVtx', 'fChi2TopoXicPlusToPVBeforeConstraint', 'fChi2TopoXicPlusToPV', 'fChi2TopoXiToXicPlusBeforeConstraint', 'fChi2TopoXiToXicPlus',
                 'fDcaPi0Pi1', 'fDcaXYPi0Pi1', 'fDcaPi0Xi', 'fDcaXYPi0Xi', 'fDcaPi1Xi', 'fDcaXYPi1Xi', 'fDcaXiDaughters',
                 'fPPi0', 'fPPi1', 'fPBachelorPi', 'fPPiFromLambda', 'fPPrFromLambda',
                 'fNSigTpcPiFromXicPlus0', 'fNSigTpcPiFromXicPlus1', 'fNSigTpcBachelorPi', 'fNSigTpcPiFromLambda', 'fNSigTpcPrFromLambda',
                 'fNSigTofPiFromXicPlus0', 'fNSigTofPiFromXicPlus1', 'fNSigTofBachelorPi', 'fNSigTofPiFromLambda', 'fNSigTofPrFromLambda']

## Loop over input files
applied_dfs = []
for file in input_files:
    ## Read pandas DataFrame
    df = pd.read_parquet(file, columns=BDT_variables)

    ## Remove candidates with infinite and NaN values
    if df.isna().any().any():
        print(f"Columns with NaN values in {file}:", df.columns[df.isna().any()])
        df.dropna(inplace=True)
    elif np.isinf(df).any().any():
        print(f"Columns with infinities in {file}:", df.columns[np.isinf(df).any()])
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df.dropna(inplace=True)

    ## Calculate model prediction
    prediction = model_hdl.predict(df, output_margin=False)
    column_names =['bkg_score', 'prompt_score', 'non_prompt_score']
    for i_class in range(3):
        df[column_names[i_class]] = prediction[:, i_class]
    
    applied_df = df[['fSign', 'fPt', 'fM', 'bkg_score', 'prompt_score', 'non_prompt_score']]
    applied_dfs.append(applied_df)

    ## Free memory
    del df, applied_df

# Concatenate all dataframes
applied_dataframe = pd.concat(applied_dfs, ignore_index=True)
print(f'Number of XicPlus candidates in data: {len(applied_dataframe)}')
del applied_dfs

#------------ Save output ------------
print(f'Saving output file...')
applied_dataframe.to_parquet(f'{output_directory}/{output_file_name}.parquet.gzip', compression='gzip')
print(f'Saving output file: Done.')