LOG_FILE=systMultiTrial.log
rm $LOG_FILE
exec > >(tee -a "$LOG_FILE") 2>&1

source /lustre/alice/users/lubynets/soft/qa2/bin/qa2Config.sh

LEFT_RANGES=(2.12 2.14) # 2.16 2.18)
RIGHT_RANGES=(2.42 2.40) # 2.38 2.36)
REBIN_FACTORS=(1 2) # 4 6 8)
BG_FUNCTIONS=(2 5)
N_TRIALS=20

fitsLocation=/lustre/alice/users/lubynets/runMassFit/outputs/HL/data/HF_LHC23_pass4_Thin_2P3PDstar/574294/ctbin2/syst
outputDir=/lustre/alice/users/lubynets/syst/multiFit/rawy

for lera in ${LEFT_RANGES[@]}; do
for rira in ${RIGHT_RANGES[@]}; do
for refa in ${REBIN_FACTORS[@]}; do
for bgfu in ${BG_FUNCTIONS[@]}; do

  dirPath=lera_$lera/rira_$rira/refa_$refa/bgfu_$bgfu

  cd /tmp/lubynets

  mkdir -p $dirPath/trials
  cd $dirPath/trials
  cp $fitsLocation/$dirPath/*tar .
  for I in `seq 1 $N_TRIALS`; do
    tar -xvf RawYields_Lc.$I.tar
  done
  cd ..
  multifit_qa

  cp -r smooth/RawYields_Lc $outputDir

done; done; done; done

