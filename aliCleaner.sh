#!/bin/bash

LABLES_TO_SAVE=('latest' 'latest-KF-o2')

cd sw
TEMPDIR=tempdir
cd slc9_x86-64

for dir in */ ; do
cd $dir
mkdir ../../$TEMPDIR/$dir
for LABEL in ${LABLES_TO_SAVE[@]}; do
cp -r $LABEL ../../$TEMPDIR/$dir
cp -r `readlink $LABEL` ../../$TEMPDIR/$dir
done
done
