import logging
import os
import sys
from io import StringIO

import pandas as pd
import pysam
from skgenome import tabio

from hmnqc import wrapper_depth


def read_bed(filename):
    logging.info("Read bed file : %s" % (filename,))
    regions = tabio.read_auto(filename)
    return regions


def min(fbam, fbed, foutput, cut_off):
    logging.info("Start Analysis")

    # Read bed
    regions = read_bed(fbed)

    writer = pd.ExcelWriter(foutput, engine="xlsxwriter")
    row = 0
    depth = wrapper_depth.WrapperDepth(fbam)

    for ix in regions.data.index:
        # Get data.
        chrom = regions.data.loc[ix, "chromosome"]
        start = regions.data.loc[ix, "start"]
        end = regions.data.loc[ix, "end"]
        gene = regions.data.loc[ix, "gene"]

        # Compute.
        depth.setRegion(chrom, start, end)
        depth.depth()
        depth.subset_by_depth(mode="gt", value=cut_off)
        df = depth.getMatrix()
        df["gene"] = [gene] * len(df.index)

        # Write data.
        if ix == 0:
            df.to_excel(
                writer,
                sheet_name="PositionProfondeur",
                startrow=row,
                startcol=0,
                index=False,
                header=True,
            )
            row += 1  # Add space for header.
        else:
            df.to_excel(
                writer,
                sheet_name="PositionProfondeur",
                startrow=row,
                startcol=0,
                index=False,
                header=False,
            )
        # Add space if df not empty.
        if len(df.index) > 0:
            row += len(df.index)

    logging.info("Write output : %s" % (foutput,))
    writer.save()


# Depth gene.
def target(fbam, name, regions, mode="target"):
    logging.info("Analysis : %s" % (name,))
    dfs = pd.DataFrame(columns=["gene", "roi", "mean", "min", "max"])
    depth = wrapper_depth.WrapperDepth(fbam)

    for ix in regions.data.index:

        # Get data
        chrom = regions.data.loc[ix, "chromosome"]
        start = regions.data.loc[ix, "start"]
        end = regions.data.loc[ix, "end"]
        gene = regions.data.loc[ix, "gene"]
        roi = "%s:%s-%s (%s)" % (chrom, start, end, gene)

        # Compute.
        depth.setRegion(chrom, start, end)
        depth.depth()
        # df = depth.getMatrix()
        # df["gene"] = [gene] * len(df.index)

        # Filter by default : BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP
        # try:
        #    depth = pysam.depth(*["-aa", '-r', region, fbam], split_lines=False)
        # except pysam.utils.SamtoolsError:
        #    continue
        # depth = pd.read_csv(StringIO(depth), sep='\t', header=None, names=["chromosome", "position", "profondeur"])

        dfs.at[ix, "gene"] = gene
        dfs.at[ix, "roi"] = roi
        dfs.at[ix, "mean"] = depth.statistics(mode="mean")
        dfs.at[ix, "min"] = depth.statistics(mode="min")
        dfs.at[ix, "max"] = depth.statistics(mode="max")

    dfs["mean"] = dfs["mean"].astype(float)

    # Build output.
    res = pd.DataFrame(
        columns=pd.MultiIndex.from_product([[name], ["moyenne", "minimum", "maximum"]])
    )
    mean, mini, maxi = None, None, None

    if mode == "gene":
        gr = dfs.groupby("gene")
        mean = gr["mean"].mean().apply(lambda x: int(round(x)))
        mini = gr["min"].min()
        maxi = gr["max"].max()
    else:
        dfs.set_index("roi", inplace=True)
        mean = dfs["mean"]  # .apply(lambda x: int(round(x)))
        mini = dfs["min"]
        maxi = dfs["max"]

    res[(name, "moyenne")] = mean
    res[(name, "minimum")] = mini
    res[(name, "maximum")] = maxi

    return (res,)


def target_write(foutput, df):
    logging.info("Write output : %s" % (foutput,))
    df.sort_index(inplace=True)
    writer = pd.ExcelWriter(foutput)
    sheet_name = "CouvertureGene"
    # writer headers as data frame
    df.columns.to_frame().transpose().to_excel(writer, sheet_name)
    # writer table body without headers
    df.to_excel(writer, sheet_name, header=False, startrow=1)
    writer.save()
    # writer.close()
