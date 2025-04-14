DATATYPES=(
           'background'
           'prompt'
           'nonPrompt'
)

PTS=(
     '0_2'
     '2_5'
     '5_8'
     '8_12'
     '12_20'
)

pageNum=46

for DT in ${DATATYPES[@]}; do
  FILEOUT=corrPlot.$DT.pdf
  for PT in ${PTS[@]}; do
    FILEIN=$DT.pT_$PT.pdf
    ls -l $FILEIN
    pdftk $FILEIN cat $pageNum output temp.pdf
    if [ ! -f $FILEOUT ]; then
      mv temp.pdf $FILEOUT
    else
      mv $FILEOUT temp_merged.pdf
      pdftk temp_merged.pdf temp.pdf cat output $FILEOUT
      rm temp*pdf
    fi
  done
done
