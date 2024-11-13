#!/bin/bash

LABLES_TO_SAVE=('latest' 'latest-KF-o2')

cd sw
TEMPDIR=tempdir
TARGET_DIR=slc9_x86-64
cd $TARGET_DIR

for dir in */ ; do
echo $dir
cd $dir
mkdir -p ../../$TEMPDIR/$dir
for LABEL in ${LABLES_TO_SAVE[@]}; do
REALNAME=`readlink $LABEL`
mv $LABEL ../../$TEMPDIR/$dir
echo $REALNAME
if [ ! -d ../../$TEMPDIR/$dir/$REALNAME ]; then
echo "mv $REALNAME ../../$TEMPDIR/$dir"
mv $REALNAME ../../$TEMPDIR/$dir
fi
done
cd ..
done

cd ..

rm -r $TARGET_DIR
mv $TEMPDIR $TARGET_DIR
