import argparse
import logging
import os
from collections import namedtuple
from concurrent import futures

import pandas as pd
from hmnqc import depth, extractvcf, infersexe, quality
from hmnqc._version import __version__
from hmnqc.utils import abort

AP = argparse.ArgumentParser(
    description="Quality information for targeted DNA analysis",
    epilog="See online documentation: ",
)
AP_subparsers = AP.add_subparsers(help="Sub-commnands (use with -h for more info)")


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

    if args.input_fastq_forward_raw:
        data["files"]["fastq"]["forward_raw"] = args.input_fastq_forward_raw
    if args.input_fastq_forward_trim:
        data["files"]["fastq"]["forward_trim"] = args.input_fastq_forward_trim
    if args.input_fastq_reverse_raw:
        data["files"]["fastq"]["reverse_raw"] = args.input_fastq_reverse_raw
    if args.input_fastq_reverse_trim:
        data["files"]["fastq"]["reverse_trim"] = args.input_fastq_reverse_trim
    if len(data["files"]["fastq"].keys()) > 0:
        isFastqCompute = True
    ks = list(data["files"]["fastq"].keys())
    if ("forward_raw" in ks and "forward_trim" in ks) or (
        "reverse_raw" in ks and "reverse_trim" in ks
    ):
        isFastqTrimmed = True
    del ks

    # Bam
    BamFile = namedtuple("BamFile", ["filename", "mode"])
    file_bam_mode = None
    if args.input_sample_bam:
        data["files"]["bam"] = args.input_sample_bam
        if os.path.basename(data["files"]["bam"]).endswith(".sam"):
            data["files"]["bam_format"] = "r"
        elif os.path.basename(data["files"]["bam"]).endswith(".bam"):
            data["files"]["bam_format"] = "rb"
        else:
            abort("File input doesn't recognize extension '.sam' or '.bam'")
        file_bam_mode = BamFile(data["files"]["bam"], data["files"]["bam_format"])
        isBamCompute = True
    coverage_cut_off = args.parameter_coverage_cut_off.split(",")
    if args.input_sample_bed:
        data["files"]["bed"] = args.input_sample_bed
        isBedCompute = True
    # Vcf
    if args.input_sample_vcf:
        data["files"]["vcf"] = args.input_sample_vcf
        isVcfCompute = True

    data["name"] = args.input_sample_name
    data["files"]["output"] = args.output_hmnqc_json

    logging.info("Name sample : %s" % (data["name"]))
    logging.info(
        "Compute FASTQ files : %s - BAM file : %s with BED file : %s - VCF file : %s"
        % (isFastqCompute, isBamCompute, isBedCompute, isVcfCompute)
    )
    logging.info("Threads used : %s" % (args.parameter_threads))

    data["statistics"] = {}

    # FASTQ
    logging.info("Compute FASTQ files")

    if isFastqCompute:
        data["statistics"]["fastq"] = {}

        with futures.ProcessPoolExecutor(args.parameter_threads) as pool:
            args_iter = (
                (origine, filename)
                for origine, filename in data["files"]["fastq"].items()
            )
            for result in pool.map(_parallel_quality_fq, args_iter):
                data["statistics"]["fastq"][result[0]] = result[1]

        # Trimming
        if isFastqTrimmed:
            total_before = data["statistics"]["fastq"].get("forward_raw", {}).get(
                "total", 0
            ) + data["statistics"]["fastq"].get("reverse_raw", {}).get("total", 0)
            total_after = data["statistics"]["fastq"].get("forward_trim", {}).get(
                "total", 0
            ) + data["statistics"]["fastq"].get("reverse_trim", {}).get("total", 0)
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
            file_bam_mode, args.parameter_threads
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
    quality.write(args.output_hmnqc_json, data)


P_quality = AP_subparsers.add_parser("quality", help=_cmd_quality.__doc__)
P_quality.add_argument("--input-sample-name", required=True, help="Name of sample")
P_quality.add_argument("--parameter-threads", type=int, default=1, help="Threads")
P_quality.add_argument("--output-hmnqc-json", required=True, help="Output json file")
# Fastq.
P_quality_fastq = P_quality.add_argument_group("Quality of fastq files")
P_quality_fastq.add_argument(
    "-iffr", "--input-fastq-forward-raw", help="Fastq forward before trimming"
)
P_quality_fastq.add_argument(
    "-ifrr", "--input-fastq-reverse-raw", help="Fastq reverse before trimming"
)
P_quality_fastq.add_argument(
    "-ifft", "--input-fastq-forward-trim", help="Fastq forward after trimming"
)
P_quality_fastq.add_argument(
    "-ifrt", "--input-fastq-reverse-trim", help="Fastq reverse after trimming"
)
# Bam.
P_quality_bam = P_quality.add_argument_group("Quality of bam file")
P_quality_bam.add_argument("--input-sample-bam", help="Bam file to analyse")
P_quality_bam.add_argument("--input-sample-bed", help="Bed file to compute coverage")
P_quality_bam.add_argument(
    "--parameter-coverage-cut-off",
    help="Cut off coverage to compute",
    default="10,20,30,40",
)
# Vcf.
P_quality_vcf = P_quality.add_argument_group("Quality of vcf file")
P_quality_vcf.add_argument("--input-sample-vcf", help="Vcf file to analyse")
P_quality.set_defaults(func=_cmd_quality)


# Check Sexe.
def _cmd_infer_sexe(args):
    """Create Html Report with Fastq"""
    df = infersexe.analysis(
        args.input_sample_bam, args.input_sample_bed, args.parameter_threads
    )
    if args.parameter_verification:
        df = infersexe.verification(df)
    if args.parameter_simple:
        df.drop(columns=["chr autosome", "chr x", "chr y"], inplace=True)
    infersexe.write(args.output_hmnqc_xlsx, df, args.paramter_simple)
    logging.info("Analysis is finished")


P_infer_sexe = AP_subparsers.add_parser("infersexe", help=_cmd_infer_sexe.__doc__)
P_infer_sexe.add_argument("-i", "--input-sample-bam", nargs="+", help="Bam files")
P_infer_sexe.add_argument("-b", "--input-sample-bed", required=True, help="Bed File")
P_infer_sexe.add_argument(
    "-o", "--output-hmnqc-xlsx", required=True, help="Excel output"
)
P_infer_sexe.add_argument(
    "-t", "--parameter-threads", default=1, type=int, help="Threads"
)
P_infer_sexe.add_argument(
    "--parameter-verification",
    action="store_true",
    help="If sexe can be guess with name sample, it will be checked again coverage",
)
P_infer_sexe.add_argument(
    "--parameter-simple",
    action="store_true",
    help="Mean coverage for ChrAutosome, ChrX and ChrY are shown by default",
)
P_infer_sexe.set_defaults(func=_cmd_infer_sexe)


# Check Snps.
def _cmd_extract_vcf(args):
    logging.info("Start analysis")
    vcf_reference = args.input_reference_vcf
    if not os.path.isfile(vcf_reference):
        abort("No file of reference provided")

    df = None
    if args.parameter_mode == "snp":
        df = extractvcf.snp(args.input_sample_vcf, vcf_reference)
    elif args.parameter_mode == "vaf":
        df = extractvcf.vaf(
            args.input_sample_vcf, vcf_reference, args.parameter_variant_caller
        )
    else:
        raise ValueError("Mode is unknown")
    logging.info("Write data to : %s" % (args.output_hmnqc_xlsx,))
    extractvcf.write(args.output_hmnqc_xlsx, df, args.parameter_mode)


P_extract_vcf = AP_subparsers.add_parser("extractvcf", help=_cmd_extract_vcf.__doc__)
P_extract_vcf.add_argument("-i", "--input-sample-vcf", nargs="+", help="Vcf files")
P_extract_vcf.add_argument(
    "-m",
    "--parameter-mode",
    choices=extractvcf.MODES,
    default=extractvcf.MODES[0],
    help="Which information to extract",
)
P_extract_vcf.add_argument(
    "--input-reference-vcf", required=True, help="Vcf file with SNPs to extract"
)
P_extract_vcf.add_argument(
    "-o", "--output-hmnqc-xlsx", required=True, help="Excel output"
)
P_extract_vaf = P_extract_vcf.add_argument_group("Extract VAF")
P_extract_vaf.add_argument(
    "--parameter-variant-caller",
    choices=extractvcf.CALLERS,
    default=extractvcf.CALLERS[0],
    help="From which caller VCF samples are created",
)
P_extract_vcf.set_defaults(func=_cmd_extract_vcf)


# Depth Min.
def _cmd_depth_min(args):
    logging.info("Start analysis")
    depth.min(
        args.input_sample_bam,
        args.input_sample_bed,
        args.output_hmnqc_xlsx,
        args.parameter_cut_off,
    )


P_depth_min = AP_subparsers.add_parser("depthmin", help=_cmd_depth_min.__doc__)
P_depth_min.add_argument("-i", "--input-sample-bam", required=True, help="Bam file")
P_depth_min.add_argument("-b", "--input-sample-bed", required=True, help="Bed File")
P_depth_min.add_argument(
    "-o", "--output-hmnqc-xlsx", required=True, help="Excel output"
)
P_depth_min.add_argument(
    "-c",
    "--parameter-cut-off",
    default=30,
    type=int,
    help="Minimal CutOff Depth to report",
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
    logging.info("Threads used : %s" % (args.parameter_threads,))

    # Read bed.
    regions = depth.read_bed(args.input_sample_bed)

    # Run process.
    df = pd.DataFrame()
    with futures.ProcessPoolExecutor(args.parameter_threads) as pool:
        args_iter = (
            (filename, _get_basename(filename), regions, args.parameter_mode)
            for filename in args.input_sample_bam
        )
        for x in pool.map(_parallel_depth_target, args_iter):
            df = pd.concat([df, x[0]], axis=1)
    # Write output
    depth.target_write(args.output_hmnqc_xlsx, df)


P_depth_target = AP_subparsers.add_parser("depthtarget", help=_cmd_depth_target.__doc__)
P_depth_target.add_argument("-i", "--input-sample-bam", nargs="+", help="Bam file")
P_depth_target.add_argument("-b", "--input-sample-bed", required=True, help="Bed File")
P_depth_target.add_argument(
    "-o", "--output-hmnqc-xlsx", required=True, help="Xlsx output"
)
P_depth_target.add_argument(
    "-t", "--parameter-threads", type=int, default=1, help="Threads used"
)
P_depth_target.add_argument(
    "-m",
    "--parameter-mode",
    choices=["target", "gene"],
    default="target",
    help="Concatenate between genes ?",
)
P_depth_target.set_defaults(func=_cmd_depth_target)


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
