PT_MAX=7
CT_MAX=0
DO_CREATE_BLANKS=true

for CC in 'ccCorrBg' 'ccBoth'; do
  for PT in `seq 0 $PT_MAX`; do
    for CT in `seq 0 $CT_MAX`; do
      if [ -f ${CC}_pt${PT}_ct${CT}.pdf ]; then
        pdftoppm -png -cropbox -singlefile ${CC}_pt${PT}_ct${CT}.pdf ${CC}_pt${PT}_ct${CT}
      elif [ $DO_CREATE_BLANKS == "true" ]; then
        convert -size 1182x1182 xc:white ${CC}_pt${PT}_ct${CT}.png
      else
        echo "Missing ${CC}_pt${PT}_ct${CT}.pdf"
      fi
    done
  done
done
