# HmnQc

[![Github Version](https://img.shields.io/github/v/release/guillaume-gricourt/HmnQc?display_name=tag&sort=semver)](version) [![Conda Release](https://img.shields.io/conda/vn/bioconda/hmnqc.svg)](https://anaconda.org/bioconda/hmnqc)  
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)  
[![GitHub Super-Linter](https://github.com/guillaume-gricourt/HmnQc/workflows/Tests/badge.svg)](https://github.com/marketplace/actions/super-linter) [![DOI](https://zenodo.org/badge/582107760.svg)](https://zenodo.org/badge/latestdoi/582107760)  

Compute differents metrics about quality, check identity and coverage from high-throughput sequencing provided by targeted NGS.

## Install

```sh
conda install -c bioconda hmnqc
```

## Running

Software is available by `hmnqc <command> <parameters>`

## Quality Metrics

Extract statistics to produce a JSON file

Compute raw statistics from FASTQ, BAM or VCF files

```sh
hmnqc quality \
    --input-sample-name <name of the sample> \
    --output-hmnqc-json <Json file> \
    --input-fastq-forward-raw <FASTQ file> \
    --input-fastq-reverse-raw <FASTQ file> \
    --input-fastq-forward-trim <FASTQ file> \
    --input-fastq-forward-trim <FASTQ file> \
    --input-sample-bam <BAM file> \
    --input-sample-bed <BED file> \
    --input-sample-vcf <VCF file>
```

## Coverage Metrics

### Position not covered

Extract position not covered under customizable cut off

```sh
hmnqc depthmin \
    --input-sample-bam <BAM file> \
    --input-sample-bed <BED file> \
    --parameter-cut-off <int> \
    --output-hmnqc-xlsx <XLXS file>
```

### Coverage of bed file

Compute statistics of coverage from a bed file

```sh
hmnqc depthtarget \
    --input-sample-bam <BAM file> \
    --input-sample-bed <BED file> \
    --parameter-mode target \
    --ouput-hmnqc-xlsx <XLSX file>
```-

## Check identiy

### Infer sexe of samples

Infer sexe from BAM files and BED file to produce XLSX file.

```sh
hmnqc infersexe \
    --input-sample-bam <BAM file> \
    --input-sample-bed <BED file> \
    --output-hmnqc-xlsx <XLSX file>
```

### Extract SNPs

Extract SNPs in VCF file from BAM files.

```sh
hmnqc extractvcf \
    --input-sample-bam <BAM file> \
    --input-reference-vcf <VCF file> \
    --output-hmnqc-xlsx <XLSX file>
```

## Built with these main libraries

* [CNVkit](https://github.com/etal/cnvkit) - Powerful library
* [pysam](https://github.com/pysam-developers/pysam) - Essential library to work with BAM and VCF files
* [biopython](https://github.com/biopython/biopython) - Essential library to work with FASTQ files
* [pandas](https://github.com/pandas-dev/pandas) - Essential dataframe object
