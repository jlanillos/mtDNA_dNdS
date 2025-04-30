# mtDNA_dNdS

This repository contains scripts to compute the **dN/dS ratio** (non-synonymous to synonymous mutation rate) for **somatic mitochondrial DNA (mtDNA) mutations** using **whole-genome sequencing (WGS)** data from **metastatic tumors**, specifically the Hartwig Medical Foundation dataset.

The goal of this pipeline is to assess selection pressures acting on mtDNA mutations across cancer samples by combining somatic mutation calls, mitochondrial copy number estimation, and coverage normalization.

The data source is the Hartwig Medical Fountation dataset, which is a set of more than 5000 pan-cancer metastatic tumor samples. This data is stored in Google buckets, that researchers can access upon access request. Since only mtDNA reads were ncessary, the strategy was to extract smaller BAM files containing such data into a personal bucket and download it (ca 1TB). Then, the rest of the analysis (variant calling, coverage, QC) was done in a workstation.

---

## 📁 Repository Structure

### Main Scripts

- **`calculate_dNdS.py`**  
  Computes the dN/dS ratio for mtDNA somatic variants by classifying mutations as synonymous or nonsynonymous and normalizing by the respective substitution space.

- **`getMTmutDB.py`**  
  Parses VCF or mutation call files to generate a clean database of mitochondrial somatic mutations suitable for downstream analysis.

- **`getmtdnaCNVs.py`**  
  Estimates mtDNA copy number variations from WGS data. This step is essential for normalizing mutation burden by mitochondrial genome abundance.

### Coverage Extraction Utilities

These scripts extract per-sample mitochondrial chromosome (chrM) coverage using `samtools` and Google Cloud infrastructure:

- **`samtools_get_chrM_sample_ids.sh`**  
  Extracts chrM coverage from BAM/CRAM files for a list of sample IDs. Meant for serial execution.

- **`parallel_samtools_get_chrM_sample_ids.sh`**  
  A parallelized version of the above script for more efficient batch processing using `GNU parallel`.

- **`getWGSCoveragegcloud.sh`**  
  Downloads and processes WGS files stored in Google Cloud using `gcloud` CLI and extracts mtDNA coverage.

- **`getWGSCoveragegcloud_bucket.sh`**  
  A bucket-specific version of the previous script, designed to handle public or shared GCP buckets containing WGS data.

---

## 📦 Dependencies

- Python ≥ 3.7  
- `pysam`, `pandas`, `numpy` for Python scripts  
- `samtools`, `gcloud` CLI, and `GNU parallel` for shell scripts

Install Python dependencies with:

```bash
pip install -r requirements.txt
