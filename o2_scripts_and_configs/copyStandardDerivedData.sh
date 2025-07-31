##########################################################################################################
# To run this script:                                                                                    #
# - click on the "Copy all output directories" button on the "Submitted jobs" tab of the train page,     #
#   and paste the contents into a file named "HyperloopOutputDirectories.txt" in the folder where the    #
#   AO2Ds are supposed to be stored                                                                      #
##########################################################################################################

export HYPERLOOP_OUTPUT_DIRECTORIES="$1"

if [[ -z "$HYPERLOOP_OUTPUT_DIRECTORIES" ]]; then
   echo "Not enough arguments! Please use"
   echo "./copyStandardDerivedData.sh HyperloopOutputDirectories.txt"
   exit 1
fi

source /lustre/alice/users/lubynets/.export_tokens.sh

apptainer shell -B /lustre -B /scratch /lustre/alice/containers/singularity_base_o2compatibility.sif << \EOF
alienv -w /scratch/alice/lubynets/alice/sw enter O2Physics/latest

for path in `sed 's/,/\n/g' ${HYPERLOOP_OUTPUT_DIRECTORIES}`; do
   # Extract last directory name from $path
   last_dir=$(basename $path)

   if [[ $(alien_find -r $path AOD/[0-9]*/AO2D.root) ]]; then
      echo "found merged files in AOD subfolder for $path, starting copy"
      alien_cp -f -parent 0 -name ends_AOD/[0-9]*/AO2D.root  $path file:.
      alien_cp -f -parent 0 -name ends_AOD/[0-9]*/AnalysisResults.root  $path file:.
   else
      echo "No AOD subfolder found in $path, falling back to parent directory"
      alien_cp -f -parent 0 -name ends_[0-9]*/AO2D.root $path file:.
      alien_cp -f -parent 0 -name ends_[0-9]*/AnalysisResults.root $path file:.
   fi
done

# cleanup existing AO2D list if present
if [ -f "localAO2DList.txt" ]; then 
   echo "removing previous AO2D list"
   rm localAO2DList.txt
fi
if [ -f "localAnalysisResultsList.txt" ]; then
   echo "removing previous AnalysisResults list"
   rm localAnalysisResultsList.txt
fi

echo "creating new AO2D list"
find ./ -name "AO2D.root" -exec readlink -f {} \; > localAO2DList.txt


echo "creating new AnalysisResults list"
find ./ -name "AnalysisResults.root" -exec readlink -f {} \; > localAnalysisResultsList.txt
EOF
