#!/bin/bash

set -o errexit
set -o nounset

# Get positional argument $1. It's the input CRAM file
file=$1

# This is necessary to do before hand and make Samtools work:
#export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
#gcloud auth application-default login
#export GCS_REQUESTER_PAYS_PROJECT=mtdna-mtx-pancan

# Project
PROJECT=mtdna-mtx-pancan # ${projectId} #
# Get output bucket
BUCKET=gs://mtdna-hmf # ${my_bucket} # gs://mtdna-hmf
# Get sample ID
SAMPLEID=$(echo ${file} | rev | cut -f1 -d '/' | rev | cut -f1 -d '_')


#Extract mitochondrial tumor CRAM and CRAI
# Get paths to sample's CRAM and CRAI
samplecopy=$BUCKET/$SAMPLEID.cram
samplecraicopy=$BUCKET/$SAMPLEID.cram.crai
# Copy files to personal bucket
gsutil -u $PROJECT cp ${file} $samplecopy
gsutil -u $PROJECT cp ${file}.crai $samplecraicopy

# Set output file name
outcram=$BUCKET/$SAMPLEID.mt.cram

# Extract mitochondrial DNA aligned reads from CRAM files
samtools view -C -o ${outcram} ${samplecopy} MT
rm $SAMPLEID.cram.crai # A CRAI file is created by default, I manually remove it

# Remove input files from Bucket
gsutil -u $PROJECT rm $samplecopy
gsutil -u $PROJECT rm $samplecraicopy
