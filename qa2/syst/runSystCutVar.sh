LOG_FILE=systCutVar.log
rm $LOG_FILE
exec > >(tee -a "$LOG_FILE") 2>&1

source /usr/local/install/bin/thisroot.sh

O2PHYSICS_DIR=/scratch/alice/lubynets/alice2/O2Physics
WORK_DIR=/lustre/alice/users/lubynets/syst/cutVar
EXE_NAME=configCutVarWriter
CONFIG_TEMPLATE=config_cutvar.json

LEFT_RANGES=(2.12 2.14) # 2.16 2.18)
RIGHT_RANGES=(2.42 2.40) # 2.38 2.36)
REBIN_FACTORS=(1 2) # 4 6 8)
BG_FUNCTIONS=(2 5)

CT_LO=1
CT_HI=5

smoothFitsLocation=/tmp/lubynets

for lera in ${LEFT_RANGES[@]}; do
for rira in ${RIGHT_RANGES[@]}; do
for refa in ${REBIN_FACTORS[@]}; do
for bgfu in ${BG_FUNCTIONS[@]}; do

  dirPath=lera_$lera/rira_$rira/refa_$refa/bgfu_$bgfu

  echo $dirPath
  date
  echo

  mkdir -p $WORK_DIR/$dirPath
  cd $WORK_DIR/$dirPath

  ln -s $WORK_DIR/$CONFIG_TEMPLATE $CONFIG_TEMPLATE

  $WORK_DIR/$EXE_NAME $CONFIG_TEMPLATE $dirPath/RawYields_Lc

  for ICT in `seq $CT_LO $CT_HI`; do
    python3 $O2PHYSICS_DIR/PWGHF/D2H/Macros/compute_fraction_cutvar.py config_cutvar_ct${ICT}.json
  done

  /u/lubynets/macros_on_git/qa2/postscripts/mergeIndividualCutVarOutputs.sh $CT_LO $CT_HI

done; done; done; done
