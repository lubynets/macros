#!/bin/bash

# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

# @brief Bash script to produce derived data AnalysisResults_trees.root with 2-prong mini skims from Run 3 real-data input for the D0 mini task
#
# The input AO2D.root file is expected in the working directory.
#
# @author Vít Kučera <vit.kucera@cern.ch>, Inha University
# @date 2023-10-25

source /lustre/alice/users/lubynets/.export_tokens.sh

# log file where the terminal output will be saved
LOGFILE="stdout.log"

# directory of this script
DIR_THIS="$(dirname "$(realpath "$0")")"

# O2 configuration file (in the same directory)
JSON="$DIR_THIS/dpl-config_tasklc.json"

# command line options of O2 workflows
OPTIONS="-b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --resources-monitoring 2 --aod-parent-access-level 1 --aod-writer-keep AOD/HFCAND3PBASE/0,AOD/HFCAND3PMCGEN/0,AOD/HFCAND3PMCREC/0,AOD/HFSELLC/0,DYN/HFCAND3PEXT/0,AOD/HFCANDLCFULL/0,AOD/HFCANDLCLITE/0,AOD/HFCOLLIDLCLITE/0,AOD/HFCANDLCFULLEV/0,AOD/HFCANDLCFULLP/0,AOD/HFCAND3PKF/0,AOD/HFCANDLCKF/0,AOD/HFCANDLCMC/0"
#
# execute the mini task workflow and its dependencies
# shellcheck disable=SC2086 # Ignore unquoted options.
# o2-analysis-hf-task-lc $OPTIONS | \
o2-analysis-hf-tree-creator-lc-to-p-k-pi $OPTIONS | \
o2-analysis-multiplicity-table $OPTIONS | \
o2-analysis-centrality-table $OPTIONS | \
o2-analysis-hf-candidate-selector-lc $OPTIONS | \
o2-analysis-pid-tpc $OPTIONS | \
o2-analysis-pid-tpc-base $OPTIONS | \
o2-analysis-pid-tof-full $OPTIONS | \
o2-analysis-pid-tof-base $OPTIONS | \
o2-analysis-hf-pid-creator $OPTIONS | \
o2-analysis-hf-candidate-creator-3prong $OPTIONS | \
o2-analysis-timestamp $OPTIONS | \
o2-analysis-event-selection $OPTIONS | \
o2-analysis-mccollision-converter $OPTIONS | \
o2-analysis-track-propagation $OPTIONS \
> "$LOGFILE" 2>&1

# OPTIONS="-b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --resources-monitoring 2"
# o2-analysis-hf-task-lc $OPTIONS > "$LOGFILE" 2>&1

# report status
rc=$?
if [ $rc -eq 0 ]; then
  echo "No problems!"
else
  echo "Error: Exit code $rc"
  echo "Check the log file $LOGFILE"
  exit $rc
fi
