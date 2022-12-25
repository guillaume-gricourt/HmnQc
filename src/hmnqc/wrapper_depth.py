import logging
from io import StringIO

import pandas as pd
import pysam


class WrapperDepth:
    COLS = ["contig", "pos", "depth"]

    def __init__(self, filename, contig="", start=0, end=0):
        self._filename = filename
        self._contig = contig
        self._start = start
        self._end = end
        self._matrix = pd.DataFrame(columns=self.COLS)
        self._alignment_file = pysam.AlignmentFile(self._filename)
        if not self._alignment_file.has_index():
            logging.warn("Index %s" % (self._alignment_file,))
            pysam.index(self._alignment_file)

    def getMatrix(self):
        return self._matrix

    def setMatrix(self, matrix):
        self._matrix = matrix

    def getRegion(self):
        return "%s:%s-%s" % (self._contig, self._start, self._end)

    def setRegion(self, contig="", start=0, end=0):
        self._contig = contig
        self._start = start
        self._end = end

    def isInit(self):
        return (
            self._contig == ""
            or self._start == 0
            or self._end == 0
            or self._matrix.empty
        )

    def clean(self):
        self._contig = ""
        self._start = 0
        self._end = 0
        self._matrix = pd.DataFrame(columns=self.COLS)

    def depth(self, select={}):
        base_quality_threshold = select.get("base_quality_threshold", 0)
        mapping_quality_threshold = select.get("base_quality_threshold", 0)
        depth = pysam.depth(
            "-aa",
            "-r",
            self.getRegion(),
            "-q",
            str(base_quality_threshold),
            "-Q",
            str(mapping_quality_threshold),
            "-d",
            "1000000",
            self._filename,
            split_lines=False,
        )
        self._matrix = pd.read_csv(
            StringIO(depth), sep="\t", header=None, names=self.COLS
        )

    def pileup(self, select={}, group_by={}):
        # Init.
        base_quality_threshold = select.get("base_quality_threshold", 0)
        mapping_quality_threshold = select.get("base_quality_threshold", 0)
        isWarn = False

        for pileupcolumn in self._alignment_file.pileup(region=self.getRegion()):
            count = set()
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_refskip:
                    # Select
                    if (
                        not pileupread.is_del
                        and pileupread.alignment.get_forward_qualities()[
                            pileupread.query_position
                        ]
                        < base_quality_threshold
                    ):
                        continue
                    if pileupread.alignment.mapping_quality < mapping_quality_threshold:
                        continue

                    # Group_by
                    item = None
                    if len(group_by.keys()) == 0:
                        item = hash(pileupread.alignment)
                    elif "tag" in group_by.keys():
                        if not pileupread.alignement.has_tag(group_by["tag"]):
                            if not isWarn:
                                logging.warn(
                                    "Tag %s doesn't exist in aligned read"
                                    % (group_by["tag"],)
                                )
                                isWarn = True
                        else:
                            item = pileupread.alignement.get_tag(group_by["tag"])

                    if item is None:
                        raise ValueError("Unable to found group_by method")
                    count.add(item)

            its = dict(
                contig=pileupcolumn.reference_name,
                pos=pileupcolumn.reference_pos,
                depth=len(count),
            )
            self._matrix = self._matrix.append(its, ignore_index=True)

    def subset_by_depth(self, mode="gt", value=30):
        if mode == "gt":
            self._matrix = self._matrix[self._matrix["depth"] <= value]
        elif mode == "lt":
            self._matrix = self._matrix[self._matrix["depth"] >= value]

    def subset_by_coord(self, start=0):
        self._matrix = self._matrix[self._matrix["pos"] >= start]

    def statistics(self, mode="mean"):
        values = self._matrix["depth"]
        if mode == "min":
            return values.min()
        elif mode == "max":
            return values.max()
        return values.mean()

    def __del__(self):
        self._alignment_file.close()
