#!/bin/bash

FILE_HASH=qa2Version.txt
FILE_DIFF=qa2Diff.patch
for F in $FILE_HASH $FILE_DIFF; do
  if [ -f $F ]; then
    rm $F
  fi
done

SRC_DIR=${1}
cd $SRC_DIR

if [ -d "../../.git" ]; then
  GITCOMMIT=$(git rev-parse HEAD)
  GITSTATUS=$(git status --porcelain)
  echo "qa2 commit hash = ${GITCOMMIT}" >> $FILE_HASH
  if [ -z "${GITSTATUS}" ]; then
    echo "qa2 commit original = true" >> $FILE_HASH
  else
    echo "qa2 commit original = false" >> $FILE_HASH
    git diff >> $FILE_DIFF
  fi
else
  echo "Not a git repo" >> $FILE_HASH
fi
cd -

mv $SRC_DIR/$FILE_HASH .
if [ -f $SRC_DIR/$FILE_DIFF ]; then
  mv $SRC_DIR/$FILE_DIFF .
fi
