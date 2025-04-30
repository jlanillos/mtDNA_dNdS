#!/bin/bash
set -o errexit
set -o nounset
#parallel -j 30 --colsep ',' -a files.txt "bash getWGSCoveragegcloud_bucket.sh {1} {2}"


tumorid=$1
file=$2
idxfile=$file.crai


TMPFOLDER=/home/crodriguez/data500gb

BUCKET=gs://ddr225
PROJECT=mtdna-mtx-pancan # ${projectId} #
outbucket=gs://ddr225/COVERAGE
OUTFILE=$tumorid.coverage.txt

# Checking that outfile exist: gsutil -u mtdna-mtx-pancan cp gs://ddr225/COVERAGE/*txt data500gb/.
checkfile=$TMPFOLDER/$OUTFILE

if [[ -f "$checkfile" ]]; then
  echo "$checkfile exist"
else
  gsutil -u mtdna-mtx-pancan cp $file $BUCKET/CRAMS/$tumorid.cram   # Copy files from bucket of origin
  gsutil -u mtdna-mtx-pancan cp $idxfile $BUCKET/CRAMS/$tumorid.cram.crai
  samtools depth $BUCKET/CRAMS/$tumorid.cram  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > $TMPFOLDER/$OUTFILE
  if [ -s "$TMPFOLDER/$OUTFILE" ]; then
    echo "Output file is empty"
    rm $TMPFOLDER/$OUTFILE
    gsutil -u mtdna-mtx-pancan rm -f $BUCKET/CRAMS/$tumorid.cram
    gsutil -u mtdna-mtx-pancan rm -f $BUCKET/CRAMS/$tumorid.cram.crai
  else
    gsutil -u mtdna-mtx-pancan mv $TMPFOLDER/$OUTFILE $outbucket/$OUTFILE
    gsutil -u mtdna-mtx-pancan rm -f $BUCKET/CRAMS/$tumorid.cram
    gsutil -u mtdna-mtx-pancan rm -f $BUCKET/CRAMS/$tumorid.cram.crai
    echo "Finished: $tumorid"
  fi
fi
