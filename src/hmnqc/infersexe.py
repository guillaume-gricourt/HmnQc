import logging
import os
import re
import sys

import cnvlib
import numpy as np
import pandas as pd

from hmnqc.utils import isNaN

##################
##    VAR    ##
##################

REGEXPS = {"name": [r"^([-_\w]+_S\d+).*"], "check": r".*[-_]([MF])_S\d+.*"}

LABELS = {
    "f": {"name": "Femme", "color": "#d7bde2"},
    "m": {"name": "Homme", "color": "#a9cce3"},
    "u": {"name": "Undefined", "color": "#d5dbdb"},
}
RESULTS = {
    "good": {"name": "ok", "color": "#a9dfbf"},
    "bad": {"name": "erreur", "color": "#f5b7b1"},
}
##########################
##    FUNCTION    ##
##########################


def _translate_sexe(isXX):
    sexe = None
    if isXX is None:
        sexe = LABELS["u"]["name"]
    else:
        if isXX:
            sexe = LABELS["f"]["name"]
        else:
            sexe = LABELS["m"]["name"]
    return sexe


def _guess_name(filename):
    basename = os.path.basename(filename)
    name = basename
    for regexp in REGEXPS["name"]:
        match = re.search(regexp, name)
        if match:
            name = match.group(1)
            break
    if name == basename:
        logging.warn("Sample name format not found")
    return name


##########################
##    Analysis    ##
##########################


def analysis(finputs, fbed, threads=1):
    logging.info(
        "Start Analysis with samples : %s and bed file : %s"
        % (", ".join(finputs), fbed)
    )
    # First round
    COLS = ["sample", "chr autosome", "chr x", "chr y", "sexe"]
    df = pd.DataFrame(columns=COLS, index=[x for x in range(len(finputs))])
    count = True
    min_mapq = 0

    for ix, finput in enumerate(finputs):
        cnar = cnvlib.do_coverage(fbed, finput, True, min_mapq, threads)
        isXx = cnar.guess_xx()
        autosome = int(cnar.autosomes()["depth"].mean())
        chrx = cnar[cnar.chromosome == cnar._chr_x_label]["depth"].mean()
        if isNaN(chrx):
            chrx = 0
        chrx = int(chrx)
        chry = cnar[cnar.chromosome == cnar._chr_y_label]["depth"].mean()
        if isNaN(chry):
            chry = 0
        chry = int(chry)

        # Store
        df.at[ix, "sample"] = _guess_name(finput)
        df.at[ix, "chr autosome"] = autosome
        df.at[ix, "chr x"] = chrx
        df.at[ix, "chr y"] = chry
        df.at[ix, "sexe"] = _translate_sexe(isXx)

    return df


def verification(df):
    def _extract_sexe_name(x):
        verif = ""
        m = re.search(REGEXPS["check"], x["sample"])
        if m:
            sample_sexe = m.group(1)
            sample_sexe = LABELS[sample_sexe.lower()]["name"]
            if sample_sexe == x["sexe"]:
                verif = RESULTS["good"]["name"]
            else:
                verif = RESULTS["bad"]["name"]

        return verif

    df["verification"] = df.apply(_extract_sexe_name, axis=1)
    return df


############
## Output ##
############


def write(foutput, df, simple):
    logging.info("Write output : %s" % (foutput,))
    # From : https://stackoverflow.com/questions/32957441/putting-many-python-pandas-dataframes-to-one-excel-worksheet
    def color(vals):
        col = vals.name
        colors = []
        if col in ["sexe", "verification"]:
            dico = None
            if col == "sexe":
                dico = LABELS
            elif col == "verification":
                dico = RESULTS

            for val in vals:
                color = "white"
                for k, v in dico.items():
                    if v["name"] == val:
                        color = v["color"]
                colors.append("background-color: %s" % color)
        else:
            colors = [""] * len(vals)
        return colors

    writer = pd.ExcelWriter(foutput, engine="xlsxwriter")
    workbook = writer.book
    worksheet = workbook.add_worksheet("EchantillonSexe")
    writer.sheets["EchantillonSexe"] = worksheet

    df.style.apply(color).to_excel(
        writer, sheet_name="EchantillonSexe", na_rep="", index=False
    )

    writer.close()
