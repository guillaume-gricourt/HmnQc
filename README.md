# HmnQc

## Introduction

A command-line toolkit to compute differents metrics about quality, check identity and coverage from high-throughput sequencing provided by targeted NGS.

## Install

Installation is done with `pip`

```sh
git clone git@github.com:guillaume-gricourt/HmnQc.git
pip install HmnQc
```

## Running

Software is available by `HmnQc <command> <parameters>`

## Quality Metrics

Extract statistics to produce a JSON file

Compute raw statistics from FASTQ, BAM or VCF files

```sh
HmnQc quality \
    --name <name of the sample> \
    --output <Json file> \
    --fastq-forward-before <FASTQ file> \
    --fastq-reverse-before <FASTQ file> \
    --bam <BAM file> \
    --bed <BED file> \
    --vcf <VCF file>
```

## Coverage Metrics

### Position not covered

Extract position not covered under customizable cut off

```sh
HmnQc depthmin \
    -i <BAM file> \
    -b <BED file> \
    --cut-off <int> \
    -o <XLXS file>
```

### Coverage of bed file

Compute statistics of coverage from a bed file

```sh
HmnQc depthtarget \
    -i <BAM file> \
    -b <BED file> \
    -m target \
    -o <XLSX file>
```

## Check identiy

### Infer sexe of samples

Infer sexe from BAM files and BED file to produce XLSX file.

```sh
HmnQc infersexe \
    -i <BAM file> \
    -b <BED file> \
    -o <XLSX file>
```

### Extract SNPs

Extract SNPs in VCF file from BAM files.

```sh
HmnQc extractvcf \
    -i <BAM file> \
    --vcf-reference <VCF file> \
    -o <XLSX file>
```

## Built with these main libraries

* [CNVkit](https://github.com/etal/cnvkit) - Powerful library
* [pysam](https://github.com/pysam-developers/pysam) - Essential library to work with BAM and VCF files
* [biopython](https://github.com/biopython/biopython) - Essential library to work with FASTQ files
* [pandas](https://github.com/pandas-dev/pandas) - Essential dataframe object

## Authors

* **Guillaume Gricourt**
