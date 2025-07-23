FILE_NAME_TEMPLATE="CutVarLc"

FILE_NAME_PREFIXES=("Distr_" "Eff_" "Frac_" "Unc_" "CovMatrix_" "")

NBINS=5

for PRE in "${FILE_NAME_PREFIXES[@]}"; do
  INPUT_FILES=""
  for index in `seq 1 $NBINS`; do
    INPUT_FILES="${INPUT_FILES}${PRE}${FILE_NAME_TEMPLATE}_bin_${index}.pdf "
  done
  pdftk ${INPUT_FILES}cat output ${PRE}${FILE_NAME_TEMPLATE}.merged.pdf
  rm $INPUT_FILES
done

if [ -f $FILE_NAME_TEMPLATE.merged.root ]; then
rm $FILE_NAME_TEMPLATE.merged.root
fi
hadd $FILE_NAME_TEMPLATE.merged.root ${FILE_NAME_TEMPLATE}*bin*root
