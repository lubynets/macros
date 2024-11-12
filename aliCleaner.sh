#!/bin/bash

LABLES_TO_SAVE=('latest' 'latest-KF-o2')

cd sw
TEMPDIR=tempdir
cd slc9_x86-64

for dir in */ ; do
echo $dir
cd $dir
mkdir -p ../../$TEMPDIR/$dir
for LABEL in ${LABLES_TO_SAVE[@]}; do
cp -r $LABEL ../../$TEMPDIR/$dir
REALNAME=`readlink $LABEL`
if [ ! -d ../../$TEMPDIR/$dir/$REALNAME ]; then
cp -r $REALNAME ../../$TEMPDIR/$dir
fi
done
cd ..
done
