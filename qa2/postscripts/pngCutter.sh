COMPONENT=('X' 'Y' 'Z')
VERTEX=(
#         'pv'
        'sv'
)
PULL_RES=('pull' 'res')

for VX in ${VERTEX[@]}; do
for PR in ${PULL_RES[@]}; do
for CO in ${COMPONENT[@]}; do
FILE=mc_qa2diff_${CO}${VX}_${PR}
num_pages=$(pdftk $FILE.pdf dump_data | grep NumberOfPages | awk '{print $2}')
pdftoppm -png -cropbox -singlefile -l $(($num_pages-2)) -f $(($num_pages-2)) $FILE.pdf $FILE.meanband
pdftoppm -png -cropbox -singlefile -l $(($num_pages-1)) -f $(($num_pages-1)) $FILE.pdf $FILE.meanstat
pdftoppm -png -cropbox -singlefile -l $num_pages -f $num_pages $FILE.pdf $FILE.width
done
./merger $VX $PR
mv out.png mc_qa2diff.$VX.$PR.png
done; done
