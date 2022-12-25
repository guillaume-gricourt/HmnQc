import argparse
import logging
import os
import re
import shutil
import sys
import tempfile
from collections import namedtuple
from concurrent import futures
from multiprocessing import Pool

import cnvlib
import numpy as np
import pandas as pd

from hmnqc import context, depth, extractvcf, formatvcf, infersexe, quality
from hmnqc._version import __app_name__, __version__
from hmnqc.utils import abort, read_json

AP = argparse.ArgumentParser(
    description="Quality information for targeted DNA analysis",
    epilog="See online documentation: ",
)
AP_subparsers = AP.add_subparsers(help="Sub-commnands (use with -h for more info)")

##  QUALITY
def _parallel_quality_fq(args):
    """Wrapper for parallel."""
    return list(quality.fastq(*args))


def _cmd_quality(args):

    isFastqCompute, isFastqTrimmed, isBamCompute, isBedCompute, isVcfCompute = (
        False,
        False,
        False,
        False,
        False,
    )
    data = {}
    data["files"] = {}

    # Fastq
    data["files"]["fastq"] = {}

    if args.fastq_forward_before:
        data["files"]["fastq"]["forward_before"] = args.fastq_forward_before
    if args.fastq_forward_after:
        data["files"]["fastq"]["forward_after"] = args.fastq_forward_after
    if args.fastq_reverse_before:
        data["files"]["fastq"]["reverse_before"] = args.fastq_reverse_before
    if args.fastq_reverse_after:
        data["files"]["fastq"]["reverse_after"] = args.fastq_reverse_after
    if len(data["files"]["fastq"].keys()) > 0:
        isFastqCompute = True
    ks = list(data["files"]["fastq"].keys())
    if ("forward_before" in ks and "forward_after" in ks) or (
        "reverse_before" in ks and "reverse_after" in ks
    ):
        isFastqTrimmed = True
    del ks

    # Bam
    BamFile = namedtuple("BamFile", ["filename", "mode"])
    file_bam_mode = None
    if args.bam:
        data["files"]["bam"] = args.bam
        if os.path.basename(data["files"]["bam"]).endswith(".sam"):
            data["files"]["bam_format"] = "r"
        elif os.path.basename(data["files"]["bam"]).endswith(".bam"):
            data["files"]["bam_format"] = "rb"
        else:
            abort("File input doesn't recognize extension '.sam' or '.bam'")
        file_bam_mode = BamFile(data["files"]["bam"], data["files"]["bam_format"])
        isBamCompute = True
    coverage_cut_off = args.coverage_cut_off.split(",")
    if args.bed:
        data["files"]["bed"] = args.bed
        isBedCompute = True
    # Vcf
    if args.vcf:
        data["files"]["vcf"] = args.vcf
        isVcfCompute = True

    data["name"] = args.name
    data["files"]["output"] = args.output

    logging.info("Name sample : %s" % (data["name"]))
    logging.info(
        "Compute FASTQ files : %s - BAM file : %s with BED file : %s - VCF file : %s"
        % (isFastqCompute, isBamCompute, isBedCompute, isVcfCompute)
    )
    logging.info("Threads used : %s" % (args.threads))

    data["statistics"] = {}

    # FASTQ
    logging.info("Compute FASTQ files")
    results = []

    if isFastqCompute:
        data["statistics"]["fastq"] = {}

        # with Pool(processes=args.threads) as pool:
        #    results = pool.starmap(quality.fastq, [ (origine, filename) for origine, filename in data["files"]["fastq"].items()] )

        with futures.ProcessPoolExecutor(args.threads) as pool:
            args_iter = (
                (origine, filename)
                for origine, filename in data["files"]["fastq"].items()
            )
            for result in pool.map(_parallel_quality_fq, args_iter):
                data["statistics"]["fastq"][result[0]] = result[1]

        # for result in results:
        #    data["statistics"]["fastq"][result[0]] = result[1]
        # Trimming
        if isFastqTrimmed:
            total_before = data["statistics"]["fastq"].get("forward_before", {}).get(
                "total", 0
            ) + data["statistics"]["fastq"].get("reverse_before", {}).get("total", 0)
            total_after = data["statistics"]["fastq"].get("forward_after", {}).get(
                "total", 0
            ) + data["statistics"]["fastq"].get("reverse_after", {}).get("total", 0)
            drop = total_before - total_after
            surviving_percent = int(
                100 - float(total_after) * 100 / float(total_before)
            )
            data["statistics"]["trimming"] = {
                "total": total_before,
                "surviving": total_after,
                "drop": drop,
                "surviving_percent": surviving_percent,
            }

    # Bed
    logging.info("Compute BED file")
    regions = []
    if isBedCompute:
        regions = quality.open_bed(data["files"]["bed"])

    # BAM
    logging.info("Compute BAM file")
    if isBamCompute:
        data["statistics"]["bam"] = {}
        logging.info("Compute BAM file - Flagstat")
        data["statistics"]["bam"]["flagstat"] = quality.flagstat(
            file_bam_mode, args.threads
        )
        logging.info("Compute BAM file - Insert")
        data["statistics"]["bam"]["insert"] = quality.insert(file_bam_mode)
        if len(regions) > 0:
            logging.info("Compute BAM file - Coverage")
            data["statistics"]["bam"]["coverage"] = quality.coverage(
                file_bam_mode, regions, coverage_cut_off
            )

    # VCF
    logging.info("Compute VCF file")
    if isVcfCompute:
        data["statistics"]["vcf"] = quality.vcf(data["files"]["vcf"])

    # Write output
    logging.info("Write output")
    quality.write(args.output, data)


P_quality = AP_subparsers.add_parser("quality", help=_cmd_quality.__doc__)

P_quality.add_argument("--name", required=True, help="Name of sample")
P_quality.add_argument("--threads", type=int, default=1, help="Threads")
P_quality.add_argument("--output", required=True, help="Output json file")

# Fastq.
P_quality_fastq = P_quality.add_argument_group("Quality of fastq files")
P_quality_fastq.add_argument(
    "-ffb", "--fastq-forward-before", help="Fastq forward before trimming"
)
P_quality_fastq.add_argument(
    "-ffa", "--fastq-forward-after", help="Fastq forward after trimming"
)
P_quality_fastq.add_argument(
    "-frb", "--fastq-reverse-before", help="Fastq reverse before trimming"
)
P_quality_fastq.add_argument(
    "-fra", "--fastq-reverse-after", help="Fastq reverse after trimming"
)

# Bam.
P_quality_bam = P_quality.add_argument_group("Quality of bam file")
P_quality_bam.add_argument("--bam", help="Bam file to analyse")
P_quality_bam.add_argument("--bed", help="Bed file to compute coverage")
P_quality_bam.add_argument(
    "--coverage-cut-off", help="Cut off coverage to compute", default="10,20,30,40"
)

# Vcf.
P_quality_vcf = P_quality.add_argument_group("Quality of vcf file")
P_quality_vcf.add_argument("--vcf", help="Vcf file to analyse")

P_quality.set_defaults(func=_cmd_quality)

## IDENTITY
# Check Sexe.
def _cmd_infer_sexe(args):
    """Create Html Report with Fastq"""
    df = infersexe.analysis(args.input, args.bed, args.threads)
    if args.verification:
        df = infersexe.verification(df)
    if args.simple:
        df.drop(columns=["chr autosome", "chr x", "chr y"], inplace=True)
    infersexe.write(args.output, df, args.simple)
    logging.info("Analysis is finished")


P_infer_sexe = AP_subparsers.add_parser("infersexe", help=_cmd_infer_sexe.__doc__)
P_infer_sexe.add_argument("-i", "--input", nargs="+", help="Bam files")
P_infer_sexe.add_argument("-b", "--bed", required=True, help="Bed File")
P_infer_sexe.add_argument("-o", "--output", required=True, help="Xls output")
P_infer_sexe.add_argument("-t", "--threads", default=1, type=int, help="Threads")
P_infer_sexe.add_argument(
    "--verification",
    action="store_true",
    help="If sexe can be guess with name sample, it will be checked again coverage",
)
P_infer_sexe.add_argument(
    "--simple",
    action="store_true",
    help="Mean coverage for ChrAutosome, ChrX and ChrY are shown by default",
)

P_infer_sexe.set_defaults(func=_cmd_infer_sexe)

# Check Snps.
def _cmd_extract_vcf(args):
    logging.info("Start analysis")
    vcf_reference = args.vcf_reference
    if not os.path.isfile(vcf_reference):
        abort("No file of reference provided")

    df = None
    if args.mode == "snp":
        df = extractvcf.snp(args.input, vcf_reference)
    elif args.mode == "vaf":
        df = extractvcf.vaf(args.input, vcf_reference, args.variant_caller)
    else:
        raise ValueError("Mode is unknown")
    logging.info("Write data to : %s" % (args.output,))
    extractvcf.write(args.output, df, args.mode)


P_extract_vcf = AP_subparsers.add_parser("extractvcf", help=_cmd_extract_vcf.__doc__)
P_extract_vcf.add_argument("-i", "--input", nargs="+", help="Vcf files")
P_extract_vcf.add_argument(
    "-m",
    "--mode",
    choices=extractvcf.MODES,
    default=extractvcf.MODES[0],
    help="Which information to extract",
)
P_extract_vcf.add_argument(
    "--vcf-reference", required=True, help="Vcf file with SNPs to extract"
)
P_extract_vcf.add_argument("-o", "--output", required=True, help="Xls output")

P_extract_vaf = P_extract_vcf.add_argument_group("Extract VAF")
P_extract_vaf.add_argument(
    "--variant-caller",
    choices=extractvcf.CALLERS,
    default=extractvcf.CALLERS[0],
    help="From which caller VCF samples are created",
)

P_extract_vcf.set_defaults(func=_cmd_extract_vcf)

##  DEPTH
# Depth Min.
def _cmd_depth_min(args):
    logging.info("Start analysis")
    depth.min(args.input, args.bed, args.output, args.cut_off)


P_depth_min = AP_subparsers.add_parser("depthmin", help=_cmd_depth_min.__doc__)
P_depth_min.add_argument("-i", "--input", required=True, help="Bam file")
P_depth_min.add_argument("-b", "--bed", required=True, help="Bed File")
P_depth_min.add_argument("-o", "--output", required=True, help="Xls output")
P_depth_min.add_argument(
    "-c", "--cut-off", default=30, type=int, help="Minimal CutOff Depth to report"
)

P_depth_min.set_defaults(func=_cmd_depth_min)

# Depth by target.
def _parallel_depth_target(args):
    """Wrapper for parallel."""
    return list(depth.target(*args))


def _cmd_depth_target(args):
    def _get_basename(filename):
        name = os.path.basename(filename)
        name = name.replace(".bam", "")
        return name

    logging.info("Start analysis")
    logging.info("Threads used : %s" % (args.threads,))

    # Read bed.
    regions = depth.read_bed(args.bed)

    # Run process.
    df = pd.DataFrame()
    with futures.ProcessPoolExecutor(args.threads) as pool:
        args_iter = (
            (filename, _get_basename(filename), regions, args.mode)
            for filename in args.input
        )
        for x in pool.map(_parallel_depth_target, args_iter):
            df = pd.concat([df, x[0]], axis=1)
    # Write output
    depth.target_write(args.output, df)


P_depth_target = AP_subparsers.add_parser("depthtarget", help=_cmd_depth_target.__doc__)
P_depth_target.add_argument("-i", "--input", nargs="+", help="Bam file")
P_depth_target.add_argument("-b", "--bed", required=True, help="Bed File")
P_depth_target.add_argument("-o", "--output", required=True, help="Xls output")
P_depth_target.add_argument("-t", "--threads", type=int, default=1, help="Threads used")
P_depth_target.add_argument(
    "-m",
    "--mode",
    choices=["target", "gene"],
    default="target",
    help="Concatenate between genes ?",
)
P_depth_target.set_defaults(func=_cmd_depth_target)

## Format vcf
def _parallel_format_vcf(args):
    """Wrapper for parallel."""
    return list(formatvcf.build(*args))


def _cmd_format_vcf(args):
    if not (len(args.input_annovar_txt) == len(args.input_gatk_vcf) == len(args.name)):
        abort("Arguments must be same size")
    for f in args.input_annovar_txt + args.input_gatk_vcf:
        if not os.path.isfile(f):
            abort("File %s is not a regular file" % (f,))

    logging.info("Build recurrences")
    recurrences = formatvcf.build_recurrence(args.input_annovar_txt)

    logging.info("Start analysis")
    dfs = []
    with futures.ProcessPoolExecutor(args.threads) as pool:
        args_iter = (
            (n, g, a)
            for n, g, a in zip(args.name, args.input_gatk_vcf, args.input_annovar_txt)
        )
        for x in pool.map(_parallel_format_vcf, args_iter):
            dfs.append(x)

    # Write output
    formatvcf.write(args.outdir, dfs)


P_format_vcf = AP_subparsers.add_parser("format-vcf", help=_cmd_format_vcf.__doc__)
P_format_vcf.add_argument(
    "-i", "--input-annovar-txt", nargs="+", help="Annovar tab file"
)
P_format_vcf.add_argument("-g", "--input-gatk-vcf", nargs="+", help="Vcf file - gatk")
P_format_vcf.add_argument("-n", "--input-name", nargs="+", help="Name to target")
P_format_vcf.add_argument("-o", "--outdir", required=True, help="Outdir")
P_format_vcf.add_argument("-t", "--threads", type=int, default=1, help="Threads used")

P_format_vcf.set_defaults(func=_cmd_format_vcf)


# Version.
def print_version(_args):
    """Display this program's version"""
    print(__version__)


P_version = AP_subparsers.add_parser("version", help=print_version.__doc__)
P_version.set_defaults(func=print_version)

# Help.
def print_help():
    """Display this program's help"""
    print(AP_subparsers.help)
    AP.exit()


# Main.
def parse_args(args=None):
    """Parse the command line"""
    return AP.parse_args(args=args)
