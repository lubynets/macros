#!/bin/bash

START_TIME=$SECONDS

source /lustre/alice/users/lubynets/.export_tokens.sh

# log file where the terminal output will be saved
LOGFILE="stdout.log"

INPUT_FILE="/lustre/alice/users/lubynets/ao2ds/tpc/AO2D.lite.root"

# command line options of O2 workflows
# execute the mini task workflow and its dependencies
# shellcheck disable=SC2086 # Ignore unquoted options.
o2-analysis-pid-tpc-skimscreation -b --configuration json://../configuration.json | \
o2-analysis-dq-v0-selector -b --configuration json://../configuration.json | \
o2-analysis-multiplicity-table -b --configuration json://../configuration.json | \
o2-analysis-trackselection -b --configuration json://../configuration.json | \
o2-analysis-pid-tof-merge -b --configuration json://../configuration.json | \
o2-analysis-pid-tpc -b --configuration json://../configuration.json | \
o2-analysis-lf-strangenessbuilder -b --configuration json://../configuration.json | \
o2-analysis-event-selection -b --configuration json://../configuration.json | \
o2-analysis-track-propagation -b --configuration json://../configuration.json | \
o2-analysis-ft0-corrected-table -b --configuration json://../configuration.json | \
o2-analysis-pid-tpc-base -b --configuration json://../configuration.json | \
o2-analysis-occ-table-producer -b --configuration json://../configuration.json | \
o2-analysis-timestamp -b --configuration json://../configuration.json \
--aod-file $INPUT_FILE \
--aod-writer-json ../OutputDirector.json \
> "$LOGFILE" 2>&1


# report status
rc=$?
FINISH_TIME=$SECONDS
if [ $rc -eq 0 ]; then
  echo "No problems!"
  echo "elapsed time " $(($(($FINISH_TIME-$START_TIME))/60)) "m " $(($(($FINISH_TIME-$START_TIME))%60)) "s"
else
  echo "Error: Exit code $rc"
  echo "Check the log file $LOGFILE"
  echo "elapsed time " $(($(($FINISH_TIME-$START_TIME))/60)) "m " $(($(($FINISH_TIME-$START_TIME))%60)) "s"
  exit $rc
fi

