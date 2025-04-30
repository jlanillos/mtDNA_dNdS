#!/bin/bash

set -o errexit
set -o nounset


# This is necessary to make Samtools work:
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
gcloud auth application-default login
export GCS_REQUESTER_PAYS_PROJECT=mtdna-mtx-pancan

# Project
PROJECT=${projectId} # mtdna-mtx-pancan
# Get output bucket
BUCKET=${my_bucket} # gs://mtdna-hmf
# Get sample ID
SAMPLE_ID=${sampleId} # Extract from manifest.json


# TUMOR SAMPLE
#Extract mitochondrial tumor CRAM and CRAI
# Get paths to tumor CRAM and CRAI
tumorcopy=$BUCKET/$SAMPLEID.tumor.cram
tumorcraicopy=$BUCKET/$SAMPLEID.tumor.cram.crai
# Copy files to personal bucket
gsutil -u $PROJECT cp ${TUMOR_CRAM} $tumorcopy
gsutil -u $PROJECT cp ${TUMOR_CRAI} $tumorcraicopy

# Set output file name
outtumorcram=$BUCKET/$SAMPLE_ID.tumormt.cram

# Extract mitochondrial DNA aligned reads from CRAM files
# Delete all lines up to and including "FORMAT"
"samtools view " "${tumorcopy}" " MT -C > " ${outtumorcram}"

# Remove input files from Bucket
gsutil -u $PROJECT rm $tumorcopy
gsutil -u $PROJECT rm $tumorcraicopy



# NORMAL SAMPLE

#Extract mitochondrial normal CRAM and CRAI
# Get paths to normal CRAM and CRAI
normalcopy=$BUCKET/$SAMPLEID.normal.cram
normalcraicopy=$BUCKET/$SAMPLEID.normal.cram.crai
# Copy files to personal bucket
gsutil -u $PROJECT cp ${NORMAL_CRAM} $normalcopy
gsutil -u $PROJECT cp ${NORMAL_CRAI} $normalcraicopy

# Set output file name
outnormalcram=$BUCKET/$SAMPLE_ID.normalmt.cram

# Extract mitochondrial DNA aligned reads from CRAM files
# Delete all lines up to and including "FORMAT"
"samtools view " "${normalcopy}" " MT -C > " ${outnormalcram}"

# Remove input files from Bucket
gsutil -u $PROJECT rm $normalcopy
gsutil -u $PROJECT rm $normalcraicopy
