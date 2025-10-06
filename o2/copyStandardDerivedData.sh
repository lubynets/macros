##########################################################################################################
# To run this script:                                                                                    #
# - click on the "Copy all output directories" button on the "Submitted jobs" tab of the train page,     #
#   and paste the contents into a file named "HyperloopOutputDirectories.txt" in the folder where the    #
#   AO2Ds are supposed to be stored                                                                      #
##########################################################################################################

HYPERLOOP_OUTPUT_DIRECTORIES="$1"
COPY_SELECTION="${2:-2}"

if [[ -z "$HYPERLOOP_OUTPUT_DIRECTORIES" ]]; then
   echo "Not enough arguments! Please use"
   echo "./copyStandardDerivedData.sh HyperloopOutputDirectories.txt (N=2)"
   echo "N=0 for copy AO2D.root only, 1 for copy AnalysisResults.root only and 2 (default) for copy both"
   exit 1
fi

source /lustre/alice/users/lubynets/.export_tokens.sh

# Define patterns to search and copy
file_name=("AO2D.root" "AnalysisResults.root")

if [[ "$COPY_SELECTION" -eq 0 ]]; then
  start=0; end=0
elif [[ "$COPY_SELECTION" -eq 1 ]]; then
  start=1; end=1
elif [[ "$COPY_SELECTION" -eq 2 ]]; then
  start=0; end=1
else
  echo "COPY_SELECTION=${COPY_SELECTION}, while it must be 0, 1 or 2"
  exit 1
fi

for i in `seq $start $end`; do
   echo "Processing ${file_name[i]}s"

   patterns_aod="AOD/[0-9]*/"${file_name[i]}
   patterns_numerated="[0-9]*/"${file_name[i]}
   patterns_single_file=${file_name[i]}

   for path in `sed 's/,/\n/g' ${HYPERLOOP_OUTPUT_DIRECTORIES}`; do
      if [[ $path == *"/AOD" ]]; then
         path="${path%"/AOD"}"
      fi

      parent_dir=$(dirname "$path")

      if [[ $(alien_find -r "$path" ${patterns_single_file}) ]]; then
         echo "Found single (slim) file at $path - starting copy"
         pattern=$patterns_single_file
      elif [[ $(alien_find -r "$path" ${patterns_aod}) ]]; then
         echo "Found merged files under AOD at $path - starting copy"
         pattern=$patterns_aod
      elif [[ $(alien_find -r "$path" ${patterns_numerated}) ]]; then
         echo "Found unmerged files in numerated directories at $path - starting copy"
         pattern=$patterns_numerated
      else
         echo "Nothing was found at $path. Switch to the next parent_dir"
      fi

      # Get list of remote files
      remote_files=$(alien_find -r "$path" "$pattern")
      while read -r remote; do
      if [[ -z "$remote" ]]; then
         continue
      fi
      # Extract basename (filename)
      fname="${remote#"$parent_dir"/}"
      if [[ -f "$fname" ]]; then
         echo "Skipping $fname - already exists"
      else
         echo "Copying $fname from $remote"
         alien_cp $remote file:./$fname
      fi
      done <<< "$remote_files"
   done
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
