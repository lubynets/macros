import time
import argparse
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
from hipe4ml import plot_utils
from hipe4ml import analysis_utils
import sys
sys.path.append('utils')
import utils as utils

start_time = time.time()

## Define labels
LegLabels = ['Background', 'Prompt', 'Non-prompt']
OutputLabels = ['Bkg', 'Prompt', 'NonPrompt']

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
input_files_mc = config['input_files_mc']
input_files_data = config['input_files_data']
output_directory = config['output_directory']
model_directory = config['model_directory']
hyperpar_study_file = config["hyperpar_study_file"]
hyperpar_study_ntrials = config["hyperpar_study_ntrials"]
hyperpar_study_timeout = config["hyperpar_study_timeout"]
save_hyperpar_study = config["save_hyperpar_study"]
model_version = config['model_version']
pT_interval = config['pT_interval']
resonance_weights = config['resonance_weights']
training_variables = config['training_variables']
print(f'ML analysis for XicPlus candidates with {int(pT_interval[0])} < pT < {int(pT_interval[1])} GeV/c (model name: {model_version})')
## Retrieve selections
selections = config_selections['selection']
selections_string = utils.convert_sel_to_string(selections)
print('The following selections will be applied:')
print(selections_string)

#------------ Open input files ------------
sel_variables = ['total', 'finite']
signal_counts = np.zeros(len(selections) + 2)
bkg_counts = np.zeros(len(selections) + 2)

VarsToSkip = ['fDebugMcRec',
              'fInvMassXi',
              'fP', 'fPhi',
              'fPtXi', 'fPPi0', 'fPtPi0', 'fPPi1', 'fPtPi1', 'fPBachelorPi', 'fPPiFromLambda', 'fPPrFromLambda',
              'fDecayLengthNormalised', 'fDecayLengthXYNormalised',
              'fDcaXiDaughters']
## MC
print('Reading MC dataframes and applying selections...', end='\r')
McDfs = []
for input_file_mc in input_files_mc:
    ## read input file
    df = pd.read_parquet(input_file_mc)
    df.drop(columns = VarsToSkip, inplace = True)
    df.query(f'((fPt >= {pT_interval[0]}) & (fPt < {pT_interval[1]}))', inplace=True)
    signal_counts[0] += len(df)
    ## Check for infinities and NaNs
    if df.isna().any().any():
        print(f"Columns with NaN values in {input_file_mc}:", df.columns[df.isna().any()])
        df.dropna(inplace=True)
    elif np.isinf(df).any().any():
        print(f"Columns with infinities in {input_file_mc}:", df.columns[np.isinf(df).any()])
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df.dropna(inplace=True)
    signal_counts[1] += len(df)
    ## Apply selections
    if selections:
        for index, (variable, selection) in enumerate(selections.items()):
            sel_variables.append(variable)
            df.query(selection, inplace=True)
            signal_counts[index+2] += len(df)
    McDfs.append(df)
    del df
McDf = pd.concat(McDfs, ignore_index=True)
del McDfs
print('Reading MC dataframes and applying selections: Done.')
print(f'MC: {signal_counts[-1]} candidates selected out of {signal_counts[0]} candidates '
      + f'with {int(pT_interval[0])} < pT < {int(pT_interval[1])} GeV/c in total.')

## Data
print('Reading data dataframes and applying selections...', end='\r')
BkgDfs = []
for input_file_data in input_files_data:
    ## read input file
    df = pd.read_parquet(input_file_data)
    df.drop(columns = VarsToSkip, inplace = True)
    df.query(f'((fPt >= {pT_interval[0]}) & (fPt < {pT_interval[1]})) & ((fM < 2.411) | (fM > 2.525))', inplace=True)
    bkg_counts[0] += len(df)
    ## Check for infinities and NaNs
    if df.isna().any().any():
        print(f"Columns with NaN values in {input_file_data}:", df.columns[df.isna().any()])
        df.dropna(inplace=True)
    elif np.isinf(df).any().any():
        print(f"Columns with infinities in {input_file_data}:", df.columns[np.isinf(df).any()])
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df.dropna(inplace=True)
    bkg_counts[1] += len(df)
    ## Apply selections
    if selections:
        for index, (variable, selection) in enumerate(selections.items()):
            df.query(selection, inplace=True)
            bkg_counts[index+2] += len(df)
    BkgDfs.append(df)
    del df
BkgDf = pd.concat(BkgDfs, ignore_index=True)
del BkgDfs
print('Reading data dataframes and applying selections: Done.')
print(f'Data: {bkg_counts[-1]} candidates selected out of {bkg_counts[0]} candidates '
      + f'with {int(pT_interval[0])} < pT < {int(pT_interval[1])} GeV/c and M <2.411 GeV/c² or M> 2.525 GeV/c² in total.')

## Plot counter histogram
signal_efficiency = signal_counts / signal_counts[0]
bkg_efficiency = bkg_counts / bkg_counts[0]
PreselEff, ax  = plt.subplots(figsize=[9,6])
ax.step(sel_variables, signal_efficiency, where='mid', label='Signal', color='red')
ax.step(sel_variables, bkg_efficiency, where='mid', label='Background', color='blue')
ax.set_ylabel("Preselection efficiency")
ax.set_ylim(0, 1.1)
plt.xticks(rotation=90)
ax.legend()
PreselEff.tight_layout()
PreselEff.savefig(f'{output_directory}/PreselEff_pT_{int(pT_interval[0])}_{int(pT_interval[1])}.png', dpi=300, bbox_inches='tight')

#------------ Define prompt/non-prompt/background sets ------------
## Split MC set into prompt/non-promt data sets
PromptDirDf = McDf.query('(fOriginRec == 1) & (abs(fFlagMcMatchRec) == 1)')
PromptResoDf = McDf.query('(fOriginRec == 1) & (abs(fFlagMcMatchRec) == 2)')
NonPromptDirDf = McDf.query('(fOriginRec == 2) & (abs(fFlagMcMatchRec) == 1)')
NonPromptResoDf = McDf.query('(fOriginRec == 2) & (abs(fFlagMcMatchRec) == 2)')
del McDf

## print number of candidates
nPromptDir = len(PromptDirDf)
nPromptReso = len(PromptResoDf)
nNonPromptDir = len(NonPromptDirDf)
nNonPromptReso = len(NonPromptResoDf)
nBkg = len(BkgDf)
print(f'Number of prompt non-resonant candidates after preselection: {nPromptDir}')
print(f'Number of prompt resonant candidates after preselection: {nPromptReso}')
print(f'Total number of prompt candidates after preselection: {nPromptDir + nPromptReso}')
print(f'Number of non-prompt non-resonant candidates after preselection: {nNonPromptDir}')
print(f'Number of non-prompt resonant candidates after preselection: {nNonPromptReso}')
print(f'Total number of non-prompt candidates after preselection: {nNonPromptDir + nNonPromptReso}')
print(f'Number of background candidates after preselection: {nBkg}')

## Drop all columns which are not needed
VarsToDrop = ['fFlagMcMatchRec',  'fOriginRec', 'fCandidateSelFlag', 'fSign',
              'fPt', 'fY', 'fEta']
PromptDirDf.drop(columns = VarsToDrop, inplace = True)
PromptResoDf.drop(columns = VarsToDrop, inplace = True)
NonPromptDirDf.drop(columns = VarsToDrop, inplace = True)
NonPromptResoDf.drop(columns = VarsToDrop, inplace = True)
BkgDf.drop(columns = VarsToDrop, inplace = True)

#------------ Create training and test set ------------
## Combine signal and background data frames
print('Combining prompt and background data frame...', end='\r')
nBkg = (nPromptDir + nPromptReso + nNonPromptDir + nNonPromptReso)*2
CombDf = pd.concat([BkgDf[:nBkg], PromptDirDf, PromptResoDf, NonPromptDirDf, NonPromptResoDf], sort=True)
LabelsArray = [0 for iCand in range(nBkg)] + [1 for iCand in range(nPromptDir + nPromptReso)] + [2 for iCand in range(nNonPromptDir + nNonPromptReso)]
WeightsArray = ([1 for iCand in range(nBkg)] + [1 for iCand in range(nPromptDir)] + [resonance_weights for iCand in range(nPromptReso)]
                + [1 for iCand in range(nNonPromptDir)] + [resonance_weights for iCand in range(nNonPromptReso)])
print('Combining prompt and background data frame: Done.')

## Split Data frame into random training set and testing set
TrainSet, TestSet, yTrain, yTest, wTrain, wTest = train_test_split(CombDf, LabelsArray, WeightsArray, test_size=0.4, random_state=42)
## random_state: Controls the shuffling applied to the data before applying the split. Pass an int for reproducible output across multiple function calls.
TrainTestData = [TrainSet, yTrain, TestSet, yTest]

## Define variables for training
if isinstance(training_variables, list) and training_variables:
    TrainVars = training_variables
else:
    TrainVars = list(CombDf.columns)
    VarsToDropForTraining = ['fM', 'fInvMassXiPi0', 'fInvMassXiPi1',
                             'fNSigTofPiFromXicPlus0', 'fNSigTofPiFromXicPlus1', 'fNSigTofBachelorPi', 'fNSigTofPiFromLambda', 'fNSigTofPrFromLambda']
    for x in VarsToDropForTraining:
        TrainVars.remove(x)

# --------------------------------------------
#            Training and testing 
# --------------------------------------------
InputModel = xgb.XGBClassifier()
model_hdl = ModelHandler(InputModel, TrainVars)

## Optimize hyperparameters
HypeRanges = {'max_depth': (1, 3),
              'learning_rate': (0.01, 0.1),
              'n_estimators': (100, 1000),
              'min_child_weight': (1, 10), 
              'subsample': (0.8, 1.), 
              'colsample_bytree': (0.8, 1.)}

if (hyperpar_study_file != '') and (save_hyperpar_study == True):
    save_hyperpar_study = hyperpar_study_file
elif (hyperpar_study_file != '') and (save_hyperpar_study == False):
    save_hyperpar_study = None
elif (hyperpar_study_file == '') and (save_hyperpar_study == True):
    hyperpar_study_file = None
    save_hyperpar_study = f'{model_directory}/hyperpar_study_file_pT_{int(pT_interval[0])}_{int(pT_interval[1])}_{model_version}.txt'
elif (hyperpar_study_file == '') and (save_hyperpar_study == False):
    hyperpar_study_file = None
    save_hyperpar_study = None

model_hdl.optimize_params_optuna(TrainTestData, HypeRanges, cross_val_scoring="roc_auc_ovr", nfold = 5, direction="maximize",
                                 resume_study=hyperpar_study_file, save_study=save_hyperpar_study,
                                 n_trials=hyperpar_study_ntrials, timeout=hyperpar_study_timeout)

## Train model
rawoutput = False
rocaucAverage = 'macro' # 'macro' or 'weighted'
multiClassOpt = 'ovr' # 'ovo' or 'ovr'
model_hdl.train_test_model(TrainTestData, False, rawoutput, average=rocaucAverage, multi_class_opt=multiClassOpt, sample_weight=wTrain)

## Calculate predictions
yTestPred = []
yTrainPred = []
if rawoutput:
    yTrainPred = model_hdl.predict(TrainTestData[0], True)
    yTestPred = model_hdl.predict(TrainTestData[2], True)
else:
    yTrainPred = model_hdl.predict(TrainTestData[0], False)
    yTestPred = model_hdl.predict(TrainTestData[2], False)

## Save model handler
model_hdl.dump_model_handler(f'{model_directory}/BDTmodel_pT_{int(pT_interval[0])}_{int(pT_interval[1])}_{model_version}.pkl')

# --------------------------------------------
#                  Plotting 
# --------------------------------------------
# --------------- Plot BDT output probability ----------------
print('Plotting BDT ouput probability...', end='\r')
BDTprob = plot_utils.plot_output_train_test(model_hdl, TrainTestData, 100, rawoutput, LegLabels, logscale=True, density=True)
for (fig, label) in zip(BDTprob, OutputLabels):
    fig.savefig(f'{output_directory}/BDTprob{label}_pT_{int(pT_interval[0])}_{int(pT_interval[1])}.png', dpi=300, bbox_inches='tight')
print('Plotting BDT ouput probability: Done.')

# --------------- Plot ROC-AUC curve ----------------
print('Plotting ROC-AUC curve...', end='\r')
ROCcurve = plot_utils.plot_roc_train_test(TrainTestData[3], yTestPred, TrainTestData[1], yTrainPred, labels=LegLabels, average=rocaucAverage,
                                          multi_class_opt=multiClassOpt)
ROCcurve.savefig(f'{output_directory}/ROCcurve_pT_{int(pT_interval[0])}_{int(pT_interval[1])}.png', dpi=300, bbox_inches='tight')
print('Plotting ROC-AUC curve: Done.')

# --------------- Plot BDT efficiency --------------------
BDTEffArrays, threshold = analysis_utils.bdt_efficiency_array(TrainTestData[3], yTestPred, n_points=101)
BDTEff = plt.figure()
for i in range(len(LegLabels)):
    plt.plot(threshold, BDTEffArrays[i], label=f'{LegLabels[i]} efficiency')
plt.legend()
plt.xlabel(f'{LegLabels[0]}/{LegLabels[1]}/{LegLabels[2]} score')
plt.ylabel('Efficiency')
plt.title('Efficiency vs score')
plt.grid()
BDTEff.savefig(f'{output_directory}/BDTeff_pT_{int(pT_interval[0])}_{int(pT_interval[1])}.png', dpi=300, bbox_inches='tight')

# --------------- Plot feature importance ----------------
print('Plotting feature importance...', end='\r')
FeatImp = plot_utils.plot_feature_imp(TrainTestData[0], TrainTestData[1], model_hdl, labels=LegLabels, n_sample=10000, approximate=True)
## shap violin plots
for (shap, label) in zip(FeatImp[:-1], OutputLabels):
    shap.savefig(f'{output_directory}/shap_{label}_pT{int(pT_interval[0])}_{int(pT_interval[1])}.png', dpi=300, bbox_inches='tight')
## shap summary plot
shapSummary = FeatImp[-1]
shapSummary.savefig(f'{output_directory}/shapSummary_pT{int(pT_interval[0])}_{int(pT_interval[1])}.png', dpi=300, bbox_inches='tight')
print('Plotting feature importance: Done.')

# --------------- Combine direct and resonant dataframes ----------------
PromptDf = pd.concat([PromptDirDf, PromptResoDf.sample(n=int(np.round(resonance_weights*nPromptReso)), random_state=42)])
NonPromptDf = pd.concat([NonPromptDirDf, NonPromptResoDf.sample(n=int(np.round(resonance_weights*nNonPromptReso)), random_state=42)])

# --------------- Plot variable distributions ----------------
DrawVarsDict = {'Basic': ['fM', 'fInvMassXiPi0', 'fInvMassXiPi1', 'fPtSumPi0Pi1', 'fCt', 'fDecayLength', 'fDecayLengthXY'],
                'Chi2': ['fChi2Sv', 'fChi2XiVtx', 'fChi2LamVtx', 'fChi2TopoXicPlusToPVBeforeConstraint', 'fChi2TopoXicPlusToPV',
                         'fChi2TopoXiToXicPlusBeforeConstraint', 'fChi2TopoXiToXicPlus'],
                'CPA': ['fCpa', 'fCpaXY', 'fCpaXi', 'fCpaXYXi', 'fCpaLam', 'fCpaXYLam', 'fCpaLamToXi', 'fCpaXYLamToXi'],
                'DCA': ['fDcaPi0Pi1', 'fDcaXYPi0Pi1', 'fDcaPi0Xi', 'fDcaXYPi0Xi', 'fDcaPi1Xi', 'fDcaXYPi1Xi'],
                'ImpPar': ['fImpactParameterXi', 'fImpactParameterNormalisedXi', 'fImpactParameterPi0', 'fImpactParameterNormalisedPi0',
                           'fImpactParameterPi1', 'fImpactParameterNormalisedPi1', 'fMaxNormalisedDeltaIP'],
                'NSigTPC': ['fNSigTpcPiFromXicPlus0', 'fNSigTpcPiFromXicPlus1', 'fNSigTpcBachelorPi', 'fNSigTpcPiFromLambda', 'fNSigTpcPrFromLambda'],
                'NSigTOF': ['fNSigTofPiFromXicPlus0', 'fNSigTofPiFromXicPlus1', 'fNSigTofBachelorPi', 'fNSigTofPiFromLambda', 'fNSigTofPrFromLambda']}
print('Plotting variable distributions...', end='\r')
for label, VarList in DrawVarsDict.items():
    VarDist = plot_utils.plot_distr([BkgDf, PromptDf, NonPromptDf], VarList, 100, LegLabels, log=True, figsize=(11,7), alpha=0.3, grid=False, density=True)
    plt.tight_layout()
    plt.savefig(f'{output_directory}/{label}_distributions_pT_{int(pT_interval[0])}_{int(pT_interval[1])}.png', dpi=300, bbox_inches='tight')
print('Plotting variable distributions: Done.')

# --------------- Plot correlation matrices ----------------
print('Plotting correlation matrices...', end='\r')
CorrVars = ['fM'] + sorted([x for x in TrainVars], key=str.lower)
CorrMatr = plot_utils.plot_corr([BkgDf, PromptDf, NonPromptDf], CorrVars, LegLabels)
for (fig, label) in zip(CorrMatr, OutputLabels):
    fig.tight_layout()
    fig.savefig(f'{output_directory}/CorrMatr{label}_pT_{int(pT_interval[0])}_{int(pT_interval[1])}.png', dpi=300, bbox_inches='tight')
print('Plotting correlation matrices: Done.')

## delete dataframes to release memory
plt.close('all')
del BkgDf, PromptDf, PromptDirDf, PromptResoDf, NonPromptDf, NonPromptDirDf, NonPromptResoDf, TrainSet, TestSet, yTrain, yTest

print(f'Program completed.')
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total runtime: {elapsed_time / 60:.2f} minutes")