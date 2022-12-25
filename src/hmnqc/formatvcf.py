import argparse
import gzip
import json
import logging
import os
import sys
from collections import Iterable, namedtuple

import pandas as pd
import pysam
import xlsxwriter

from hmnqc.extractvcf import _record_to_vaf, _vcf_record_is_snp
from hmnqc.utils import *


def load_annovar_txt(path):
    return pd.read_csv(path, sep="\t")


def _format_cell(data):
    # Check if data is void
    if data == [] or data == ["."] or data == ".":
        return np.nan
    if isinstance(data, str):
        return data.replace(",", ";")
    if isinstance(data, Iterable):
        return ";".join([str(x) for x in data])
    return data


def load_gatk(path, name_to_target):
    vcf = pysam.VariantFile(path)
    name = [x for x in vcf.header.samples][0]
    assert name == name_to_target, "Name in file %s not matched with name given %s" % (
        name,
        name_to_target,
    )

    df = pd.DataFrame()

    formats = vcf.header.formats.keys()
    infos = vcf.header.info.keys()

    for record in vcf.fetch():
        ix = len(df.index)
        df.at[ix, "chr"] = record.chrom
        df.at[ix, "start"] = record.start
        df.at[ix, "ref"] = record.ref
        df.at[ix, "alt"] = record.alts[0]
        # Filters
        filters = record.filter.keys()
        filters = [record.filter.get(x).name for x in filters]
        df.at[ix, "filter"] = _format_cell(filters)
        # Info
        for info in infos:
            label = info
            if info in formats:
                label += "_info"
            df.at[ix, label] = _format_cell(record.info[info])
        # Format
        for form in formats:
            label = form
            if form in infos:
                label += "_format"
            df.at[ix, label] = _format_cell(record.samples[name].get(form))
        df.at[ix, "vaf"] = _record_to_vaf(record, name, "hc")
        df.at[ix, "variation_type"] = _vcf_record_is_snp(record)
    return df


def define_ident(x, mode="annovar"):
    fmt = "%s:%s[%s>%s]"
    if mode == "annovar":
        return fmt % (x["Chr"], x["Start"], x["Ref"], x["Alt"])
    else:
        return fmt % (x["chr"], x["start"], x["ref"], x["alt"])


def build_recurrence(finputs):
    recurrences = {}
    for finput in finputs:
        df = load_annovar(finput)
        df["id"] = df.apply(ident_annovar, axis=0)
        for ident, value in df["id"].value_counts().to_dict().items():
            recurrences[ident] = recurrences.get(ident, 0) + value
    return recurrences


"""
Chr    Start    End    Ref    Alt    
Func.refGene    Gene.refGene    GeneDetail.refGene    ExonicFunc.refGene    AAChange.refGene    

1000g2015aug_all    1000g2015aug_afr    1000g2015aug_amr    1000g2015aug_eur    1000g2015aug_eas    1000g2015aug_sas    
avsnp147    
cosmic70    

SIFT_score    SIFT_pred    
Polyphen2_HDIV_score    Polyphen2_HDIV_pred    Polyphen2_HVAR_score    Polyphen2_HVAR_pred    
LRT_score    LRT_pred    
MutationTaster_score    MutationTaster_pred    
MutationAssessor_score    MutationAssessor_pred    
FATHMM_score    FATHMM_pred    
PROVEAN_score    PROVEAN_pred    
VEST3_score    
CADD_raw    CADD_phred    
DANN_score    
fathmm-MKL_coding_score    fathmm-MKL_coding_pred    
MetaSVM_score    MetaSVM_pred    
MetaLR_score    MetaLR_pred    
integrated_fitCons_score    integrated_confidence_value    
GERP++_RS    
phyloP7way_vertebrate    
phyloP20way_mammalian    
phastCons7way_vertebrate    
phastCons20way_mammalian    
SiPhy_29way_logOdds    
ExAC_ALL    ExAC_AFR    ExAC_AMR    ExAC_EAS    ExAC_FIN    ExAC_NFE    ExAC_OTH    ExAC_SAS    
ExAC_nontcga_ALL    ExAC_nontcga_AFR    ExAC_nontcga_AMR    ExAC_nontcga_EAS    ExAC_nontcga_FIN    ExAC_nontcga_NFE    ExAC_nontcga_OTH    ExAC_nontcga_SAS    
CLINSIG    CLNDBN    CLNACC    CLNDSDB    CLNDSDBID    
SIFT_score    SIFT_pred    
Polyphen2_HDIV_score    Polyphen2_HDIV_pred    
Polyphen2_HVAR_score    Polyphen2_HVAR_pred    
LRT_score    LRT_pred    
MutationTaster_score    MutationTaster_pred    
MutationAssessor_score    MutationAssessor_pred    
FATHMM_score    FATHMM_pred    
RadialSVM_score    RadialSVM_pred    
LR_score    LR_pred    
VEST3_score    
CADD_raw    
CADD_phred    
GERP++_RS    
phyloP46way_placental    
phyloP100way_vertebrate    
SiPhy_29way_logOdds    
PopFreqMax    
1000G_ALL    1000G_AFR    1000G_AMR    1000G_EAS    1000G_EUR    1000G_SAS    
ExAC_ALL    ExAC_AFR    ExAC_AMR    ExAC_EAS    ExAC_FIN    ExAC_NFE    ExAC_OTH    ExAC_SAS    
ESP6500siv2_ALL    ESP6500siv2_AA    ESP6500siv2_EA    
CG46    
gnomAD_exome_ALL    gnomAD_exome_AFR    gnomAD_exome_AMR    gnomAD_exome_ASJ    gnomAD_exome_EAS    gnomAD_exome_FIN    gnomAD_exome_NFE    gnomAD_exome_OTH    gnomAD_exome_SAS    
Otherinfo1    Otherinfo2    Otherinfo3    Otherinfo4    Otherinfo5    Otherinfo6    Otherinfo7    Otherinfo8    Otherinfo9    Otherinfo10    Otherinfo11    Otherinfo12    Otherinfo13
"""
"""
identifiant
    chromosome
    position
    reference
    allele    
    variation
mutation
    gene
    nm
    c_dot
    p_dot
    type (synonyme, faux sens, non sens, mut d’épissage)
qualite
    profondeur
    VAF
    Statut (homozygote ou hétérozygote)
    Nb de patients dans le run avec ce variant
pathogenicite
?    Classification ACGM 
    sift
    polyphen_2_hvar
    polyphen_2_hdiv
    mutation_taster
    cadd
?    HGMD
?    OMIM
    clinvar
frequence
    gnomAD : nb homozygotes
    gnomAD All
    gnomAD Afr
    gnomAD SEA
    les autres populations gnomAD…       
"""
# Build df

PATHOGENICITY_SCORE = {
    "sift": {"D": "Deletere", "T": "Tolere"},
    "fathmm": {"D": "Deletere", "T": "Tolere"},
    "meta": {"D": "Deletere", "T": "Tolere"},
    "polyphen": {"D": "Probable dommage", "P": "Possible dommage", "B": "Begnin"},
    "lrt": {"D": "Deletere", "N": "Neutre", "U": "Inconnu"},
    "mutation_taster": {
        "A": "Cause de la maladie - certain",
        "D": "Cause de la maladie",
        "N": "Polymorphisme",
        "P": "Polymorphisme - certain",
    },
    "mutation_assessor": {
        "H": "Eleve",
        "M": "Moyen",
        "L": "Faible",
        "N": "Neutre",
        "H/M": "Fonctionnel",
        "L/N": "Non-Fonctionnel",
    },
}
COLUMNS_NAME = {
    # core
    "Chr": "chr",
    "Start": "start",
    "Ref": "ref",
    "Alt": "alt",
    # technical
    "vaf": "vaf",
    "mutation_type": "mutation_type",
    "DP": "profondeur",
    "filter": "fiter",
    "recurrence": "recurrence"
    # databases
}


def build(name, input_gak_vcf, input_annovar_txt, recurrences):
    Frame = namedtuple("Frame", ["name", "fiables", "candidats"])

    # Load.
    df = load_annovar_txt(input_annovar_txt)
    df_gatk = load_gatk(input_gatk_vcf)

    # Identifier.
    df["id"] = df.apply(define_ident, axis=0)
    df_gatk["id"] = df_gatk.apply(define_ident, axis=0, args=("gatk",))

    # Merge.
    df = pd.merge(df, df_gatk, on="id", how="outer")

    # Recurrence
    df["recurrence"] = df["id"]
    df.replace({"recurrence": recurrences}, inplace=True)

    # Rename.
    df.rename(columns=COLUMNS_NAME, inplace=True)

    # Drop.
    df.drop(
        columns=[x for x in df.columns if x in list(COLUMNS_NAME.values())],
        inplace=True,
    )

    # Annotate
    def _refgene(x):
        gene, nm, c_dot, p_dot = [np.nan] * 4
        if x["AAChange.refGene"] != ".":
            datas = x["AAChange.refGene"]
            datas = datas.split(";")
            gene = ["%s (%s)" % (e.split(":")[0], e.split(":")[2]) for e in datas]
            nm = [e.split(":")[1] for e in datas]
            c_dot = [e.split(":")[3] for e in datas]
            p_dot = [e.split(":")[4] for e in datas]
        else:
            gene = [e for e in x["Gene.refGene"].split(",") if e != "NONE"]
            datas = x["GeneDetail.refGene"]
            datas = datas.split(";")
            nm, c_dot = [], []
            for data in datas:
                if data == ".":
                    nm, c_dot = np.nan, np.nan
                    break
                elif "dist" in data:
                    continue
                else:
                    e = data.split(":")
                    nm.append(e[0])
                    c_dot.append(e[1])

        gene = _format_cell(gene)
        nm = _format_cell(nm)
        c_dot = _format_cell(c_dot)
        p_dot = _format_cell(p_dot)

        return gene, nm, c_dot, p_dot

    genes, nms, c_dots, p_dots = [], [], [], []
    if (
        "AAChange.refGene" in df.columns
        and "Gene.refGene" in df.columns
        and "GeneDetail.refGene" in df.columns
    ):
        df[["gene", "nm", "c_dot", "p_dot"]] = df.apply(
            _refgene, axis=1, result_type="expand"
        )
    if "ExonicFunc.refGene" in df.columns:
        df["mutation_type"] = df["ExonicFunc.refGene"].apply(_format_cell)
    if "Func.refGene" in df.columns:
        df["genomic_pos"] = df["Func.refGene"].apply(_format_cell)

    ## Pathogenicity
    def _score(x, db, db_score):
        kpred = db + "_pred"
        kscore = db + "_score"
        # Format Pred
        pred = x[kpred]
        if pred != ".":
            pred = PATHOGENICITY_SCORE[db_score][pred]
        # Format Score
        score = x[kscore]
        isScore = False
        if score != ".":
            isScore = True
        if isScore:
            score = "%.2f" % float(score)

        # Format output
        res = ""
        if isScore:
            res = "%s (%s)" % (pred, score)
        else:
            res = pred
        return res

    if "SIFT_pred" in df.columns and "SIFT_score" in df.columns:
        df["sift"] = df.apply(_score, args=("SIFT", "sift"), axis=1)
    if "Polyphen2_HVAR_pred" in df.columns and "Polyphen2_HVAR_score" in df.columns:
        df["polyphen_2_hvar"] = df.apply(
            _score, args=("Polyphen2_HVAR", "polyphen"), axis=1
        )
    if "Polyphen2_HDIV_pred" in df.columns and "Polyphen2_HDIV_score" in df.columns:
        df["polyphen_2_hdiv"] = df.apply(
            _score, args=("Polyphen2_HDIV", "polyphen"), axis=1
        )
    if "MutationTaster_pred" in df.columns and "MutationTaster_score" in df.columns:
        df["mutation_taster"] = df.apply(
            _score, args=("MutationTaster", "mutation_taster"), axis=1
        )
    if "CADD_raw" in df.columns and "CADD_phred" in df.columns:
        df["cadd"] = df.apply(_cadd, axis=1)
    if "CLINSIG" in df.columns and "CLNDBN" in dfa.columns and "CLNACC" in df.columns:
        df["clinvar"] = df.apply(_clinvar, axis=1)

    # Return
    return df
