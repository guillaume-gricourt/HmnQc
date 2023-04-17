import logging

import pandas as pd
import pysam
from pysam import VariantFile

MODES = ["snp", "vaf"]
# Callers :
#    ls : leaves
#    hc : haplotypecaller (gatk)
CALLERS = ["ls", "hc"]


def _vcf_record_equal(a, b):
    same_obj = isinstance(a, pysam.libcbcf.VariantRecord)
    contig = a.contig == b.contig
    pos = a.pos == b.pos
    ref = a.ref == b.ref
    alt = a.alts == b.alts
    return same_obj and contig and pos and ref and alt


def _vcf_record_is_snp(x):
    if x.ref is None or x.alts is None:
        return None
    if len(x.ref) == 1 and len(x.alts[0]) == 1:
        return True
    return False


def _vcf_record_type(x):
    if len(x.ref) == 1 and len(x.alts[0]) == 1:
        return "snp"
    elif len(x.ref) > len(x.alts[0]):
        return "del"
    return "ins"


def _genotype_to_label(record, gt):
    res = []  # "%s/%s" % (x.ref, x.alts[0])

    if gt[0] == 0:
        res.append(record.ref)
    elif gt[0] == 1:
        res.append(record.alts[0])
    else:
        res.append(".")

    if gt[1] == 0:
        res.append(record.ref)
    elif gt[1] == 1:
        res.append(record.alts[0])
    else:
        res.append(".")

    res = "/".join(res)
    """
    x[1]:
        if x[0] == 0:
            res = "Homozygote-Reference"
        else:
            res = "Homozygote-Alternative"
    else:
        res = "Heterozygote"
    """
    return res


def _record_to_vaf(record, name, variant_caller):
    vaf = 0
    if variant_caller == "ls":
        for label in ["ARR", "FARR"]:
            if label in record.samples[name].keys():
                vaf = record.samples[name][label]
                break
    elif variant_caller == "hc":
        vaf = (
            float(record.samples[name]["AD"][1])
            / float(record.samples[name]["DP"])
            * 100.0
        )
    return vaf


def snp(finputs, vcf_reference):
    logging.info("Parsing reference file")
    references = [x for x in VariantFile(vcf_reference).fetch()]
    logging.info("Get all samples names")
    # Get all names
    names = []
    for finput in finputs:
        vf = VariantFile(finput)
        names += vf.header.samples
    df = pd.DataFrame(index=[x + 1 for x in range(len(references))], columns=names)
    logging.info("Parsing vcf input file")
    for finput in finputs:
        vf = VariantFile(finput)
        names = vf.header.samples

        references_found = []
        # Find record present in references
        for record in vf.fetch():
            if _vcf_record_is_snp(record) is None:
                continue
            for ix, reference in enumerate(references):
                if _vcf_record_equal(reference, record):
                    for name in names:
                        df.at[ix + 1, name] = _genotype_to_label(
                            record, record.samples[name]["GT"]
                        )
                    references_found.append(ix)
                    break
        # Fill if record isn't found in reference
        references_found = set(references_found)
        references_pos = set([x for x in range(len(references))])
        references_notfound = list(references_pos.difference(references_found))
        for ix in references_notfound:
            for name in names:
                df.at[ix + 1, name] = "./."
    return df


def vaf(finputs, vcf_reference, variant_caller):
    logging.info("Parsing reference file")
    references = [x for x in VariantFile(vcf_reference).fetch()]
    logging.info("Get all samples names")
    # Get all names
    names = []
    for finput in finputs:
        vf = VariantFile(finput)
        names += vf.header.samples
    df = pd.DataFrame(index=[x.id for x in references], columns=names)
    logging.info("Parsing vcf input file")
    for finput in finputs:
        vf = VariantFile(finput)
        names = vf.header.samples

        references_found = []
        # Find record present in references
        for record in vf.fetch():
            for ix, reference in enumerate(references):
                if _vcf_record_equal(reference, record):
                    for name in names:
                        df.at[reference.id, name] = _record_to_vaf(
                            record, name, variant_caller
                        )
                    references_found.append(ix)
                    break
        # Fill if record isn't found in reference
        references_found = set(references_found)
        references_pos = set([x for x in range(len(references))])
        references_notfound = list(references_pos.difference(references_found))
        for ix in references_notfound:
            for name in names:
                df.at[references[ix].id, name] = "."
    return df


def write(foutput, df, mode):
    writer = pd.ExcelWriter(foutput, engine="xlsxwriter")
    df.to_excel(writer, sheet_name=mode, index=True, header=True)
    writer.close()
