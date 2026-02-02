PLOTS=(
'Basic_distributions'
'BDTeff'
'BDTprobBkg'
'BDTprobNonPrompt'
'BDTprobPrompt'
'Chi2Geo_distributions'
'Chi2Prim_distributions'
'CorrMatrBkg'
'CorrMatrNonPrompt'
'CorrMatrPrompt'
'DCA_distributions'
'ImpactParameter_distributions'
'NTpcTofSigma_distributions'
'Others_distributions'
'PreselEff'
'ROCcurve'
'shap_Bkg'
'shap_NonPrompt'
'shap_Prompt'
'shapSummary'
)

N=8

for PL in ${PLOTS[@]}; do
  echo "processing " $PL
  input=()
  for I in `seq 1 $N`; do
    input+=( "$I/${PL}*" )
  done
  convert "${input[@]}" "$PL.pdf"
done

for I in `seq 1 $N`; do
  echo "processing pT interval " $I
  cd $I
  convert *png complexQA.Pt_$I.pdf
  mv *pdf ../
  cd ../
done
