#!/usr/bin/env bash
shopt -s extglob

for dir in hy_*; do
  [[ -d "$dir" ]] || continue

  for FILE in AO2D.root AnalysisResults.root; do
    N=0

    FILES=$(ls "$dir"/"$FILE" 2>/dev/null)
    if [[ -n $FILES ]]; then
      N=$(($N+1))
    fi

    FILES=$(ls "$dir"/AOD/*/"$FILE" 2>/dev/null)
    if [[ -n $FILES ]]; then
      N=$(($N+1))
    fi

    FILES=$(ls "$dir"/+([0-9])/"$FILE" 2>/dev/null)
    if [[ -n $FILES ]]; then
      N=$(($N+1))
    fi

    if (( N != 1 )); then
      echo "Check dir = ${dir}; FILE = ${FILE}; N = ${N}"
    fi
  done
done
