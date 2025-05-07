PLOTS=(
        'Basic_distributions'
        'BDTeff'
        'BDTprobBkg'
        'BDTprobNonPrompt'
        'BDTprobPrompt'
        'Chi2_distributions'
        'CorrMatrBkg'
        'CorrMatrNonPrompt'
        'CorrMatrPrompt'
        'NTpcSigma_distributions'
        'PreselEff'
        'ROCcurve'
        'shap_Bkg'
        'shap_NonPrompt'
        'shap_Prompt'
        'shapSummary'
)

for PL in ${PLOTS[@]}; do
  echo "processing " $PL
  convert 1/${PL}* 2/${PL}* 3/${PL}* 4/${PL}* 5/${PL}* $PL.pdf
done

for I in `seq 1 5`; do 
  echo "processing pT interval " $I
  cd $I
  convert *png complexQA.Pt_$I.pdf
  mv *pdf ../
  cd ../
done
