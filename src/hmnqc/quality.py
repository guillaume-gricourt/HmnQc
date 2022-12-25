import gzip
import json
import re
from io import StringIO

import pandas as pd
import pysam
from Bio import SeqIO

from hmnqc.utils import calculate_percentile


class Region:
    def __init__(self, contig, start, stop, gene="", exon=""):
        self.contig = contig
        self.start = start
        self.stop = stop
        self.gene = gene
        self.exon = exon
        self.sequences = 0
        self.bases = pd.DataFrame()

    def checkInt(self, mystr):
        if type(mystr) == str:
            return int(mystr)
        else:
            return mystr

    def getRoiDepth(self, ct):
        ct = self.checkInt(ct)
        return 1 if self.sequences >= ct else 0

    def getRoiBaseTotal(self):
        return len(self.bases["profondeur"])

    def getRoiBaseCt(self, ct):
        ct = self.checkInt(ct)
        sub = self.bases[self.bases["profondeur"] >= ct]
        return len(sub.index)

    def getRoiBaseMean(self):
        return self.bases["profondeur"].mean()

    def to_str(self):
        return "%s:%s-%s" % (self.contig, self.start, self.stop)


def _calculate_stat_per_position(tab):
    baseQual = {}
    # List to dict
    for i in range(len(tab)):
        if tab[i] != 0:
            baseQual[i] = tab[i]
    # Min
    mini = sorted(baseQual.keys())[0]
    # Maxi
    maxi = sorted(baseQual.keys(), reverse=True)[0]
    # Mean
    total = 0
    for k, v in baseQual.items():
        total += k * v
    count = sum(baseQual.values())

    return {
        "min": mini,
        "max": maxi,
        "lowest": calculate_percentile(baseQual, 10),
        "highest": calculate_percentile(baseQual, 90),
        "mean": total / count,
        "median": calculate_percentile(baseQual, 50),
        "lowerQuartile": calculate_percentile(baseQual, 25),
        "upperQuartile": calculate_percentile(baseQual, 75),
    }


def _calculate_mean(dico):
    total = sum([valeur * frequence for valeur, frequence in dico.items()])
    mean = total / sum(dico.values())
    return mean


def fastq(origine, filename, threads=1):
    bases, tiles, gcs = {}, {}, {}
    total, len_max, len_min = 0, 0, -1
    len_means = {}

    fid = ""
    if filename.endswith("gz"):
        fid = gzip.open(filename, "rt")
    else:
        fid = open(filename)

    for record in SeqIO.parse(fid, "fastq"):

        seq = str(record.seq)
        quality = record.letter_annotations["phred_quality"]
        len_seq = len(quality)
        name = record.id

        total += 1
        # Length
        if len_max < len_seq:
            len_max = len_seq
        if len_min == -1:
            len_min = len_seq
        elif len_min > len_seq:
            len_min = len_seq
        len_means[len_seq] = len_means.get(len_seq, 0) + 1

        # GC
        seq_gc = seq.upper()
        count_gc = 0
        for c in seq_gc:
            if c == "C" or c == "G":
                count_gc += 1
        gc = int(round(float(count_gc) * 100 / len_seq))
        gcs[gc] = gcs.get(gc, 0) + 1

        # Tile
        tile_name = None
        splitID = name.split(":")  # If there are 7 or more fields then it's a 1.8+ file
        if len(splitID) >= 7:
            tile_name = splitID[4]
        elif len(splitID) >= 5:
            tile_name = splitID[2]
        if tile_name and tile_name not in tiles.keys():
            tiles[tile_name] = {}

        # Keep data
        for j, qual in enumerate(quality):
            # Tile
            if tile_name:
                tab = tiles.get(j, [0] * 45)
                tab[qual] += 1
                tiles[tile_name][j] = tab
            # Base
            tab = bases.get(j, [0] * 45)
            tab[qual] += 1
            bases[j] = tab
    fid.close()

    # Length stat sum up.
    len_mean = _calculate_mean(len_means)

    # Sum up GC
    gc = _calculate_mean(gcs)

    # Sum up bases stats
    for position, tab in bases.items():
        stats = _calculate_stat_per_position(tab)
        bases[position] = {}
        bases[position]["quality"] = tab
        bases[position]["statistics"] = stats
    for tile in tiles.keys():
        for position, tab in tiles[tile].items():
            stats = _calculate_stat_per_position(tab)
            tiles[tile][position] = {}
            tiles[tile][position]["quality"] = tab
            tiles[tile][position]["statistics"] = stats

    statistics = {
        "total": total,
        "length": {"minimum": len_min, "maximum": len_max, "mean": len_mean},
        "quality_per_base": bases,
        "quality_per_tile_per_base": tiles,
        "gc": gc,
    }
    return (origine, statistics)


def open_bed(filename):
    regions = []
    with open(filename) as fid:
        for line in fid:
            tab = line.split("\t")
            if len(tab) == 3:
                regions.append(Region(tab[0], float(tab[1]), float(tab[2])))
            elif len(tab) == 4:
                regions.append(Region(tab[0], float(tab[1]), float(tab[2]), tab[3]))
            elif len(tab) >= 5:
                regions.append(
                    Region(tab[0], float(tab[1]), float(tab[2]), tab[3], tab[4])
                )
            else:
                raise ValueError("Format bed doesn't recognize")
    return regions


def _open_bam_file(filename, read):
    # Load bam file
    bam = pysam.AlignmentFile(filename, read)

    # Check index
    if not bam.has_index():
        pysam.index(filename)
        bam = pysam.AlignmentFile(filename, read)

    return bam


def flagstat(file_bam_mode, threads=1):
    data = {}
    flagstat = pysam.flagstat(*["--threads", str(threads), file_bam_mode.filename])
    flagstat = flagstat.splitlines()

    for ix, line in enumerate(flagstat):
        value = int(line.split()[0])
        m = re.search(r"^\d+\s\+\s\d+\s([\w\s]+)(.*)", line)
        item = m.group(1)
        item = " ".join(item.split())
        key = None

        if "in total" == item:
            key = "reads"
        elif "secondary" == item:
            key = "secondary"
        elif "supplementary" == item:
            key = "supp"
        elif "duplicates" == item:
            key = "dup"
        elif "mapped" == item:
            key = "mapped"
        elif "paired in sequencing" == item:
            key = "pair_all"
        elif "read1" == item:
            key = "read1"
        elif "read2" == item:
            key = "read2"
        elif "properly paired" == item:
            key = "pair_good"
        elif "with itself and mate mapped" == item:
            key = "pair_map"
        elif "singletons" == item:
            key = "sgltn"
        elif "with mate mapped to a different chr" == item:
            qual = m.group(2)
            if qual == "":
                key = "diffchr"
            else:
                key = "diffhigh"
        else:
            raise ValueError("Error when parsing flagstat")

        data[key] = value
    return data


def insert(file_bam_mode):
    bam = _open_bam_file(file_bam_mode.filename, file_bam_mode.mode)
    insert = {}
    for ix, read in enumerate(bam.fetch()):
        if (
            read.is_paired
            and read.is_proper_pair
            and not read.is_unmapped
            and read.is_read1
        ):
            # Insert size
            template_length = abs(read.template_length)
            insert[template_length] = insert.get(template_length, 0) + 1
    bam.close()

    return insert


def coverage(file_bam_mode, regions, cut_offs):
    bam = _open_bam_file(file_bam_mode.filename, file_bam_mode.mode)
    # Init data
    coverage = {}
    coverage["cutoff"] = cut_offs

    # Statistics
    coverage["on_target"] = 0

    coverage["roidepth"] = {}
    coverage["roidepth"]["total"] = len(regions)
    coverage["roidepth"]["cutoff"] = {}

    coverage["roibase"] = {}
    coverage["roibase"]["nb_base"] = 0
    coverage["roibase"]["mean"] = 0
    coverage["roibase"]["cutoff"] = {}

    for ct in cut_offs:
        coverage["roidepth"]["cutoff"][ct] = 0
        coverage["roibase"]["cutoff"][ct] = 0

    coverage["roigene"] = {}
    coverage["min_per_base"] = {}
    # Get genes
    coverage["roigene"]["genes"] = {}
    genes = [x.gene for x in regions]
    for gene in genes:
        coverage["roigene"]["genes"][gene] = {}
        coverage["roigene"]["genes"][gene]["number"] = 0
        coverage["roigene"]["genes"][gene]["mean_cumulate"] = 0
        coverage["roigene"]["genes"][gene]["mean"] = 0

    # Parsing file
    for region in regions:
        # Compute by sequences
        region.sequences = bam.count(
            contig=region.contig, start=region.start, stop=region.stop
        )
        # Compute by position
        depth = pysam.depth(
            *["-aa", "-d", "0", "-r", region.to_str(), file_bam_mode.filename],
            split_lines=False
        )
        depth = pd.read_csv(
            StringIO(depth),
            sep="\t",
            header=None,
            names=["chromosome", "position", "profondeur"],
        )
        region.bases = depth.copy()

    # Calculate Stats:
    for region in regions:
        # Statistics
        coverage["on_target"] += region.sequences
        # Roi base
        coverage["roibase"]["mean"] += region.getRoiBaseMean()
        coverage["roibase"]["nb_base"] += region.getRoiBaseTotal()
        # Roi gene
        coverage["roigene"]["genes"][region.gene][
            "mean_cumulate"
        ] += region.getRoiBaseMean()
        coverage["roigene"]["genes"][region.gene]["number"] += 1

        for ct in cut_offs:
            # Roi depth
            coverage["roidepth"]["cutoff"][ct] += region.getRoiDepth(ct)
            # Roi base
            coverage["roibase"]["cutoff"][ct] += region.getRoiBaseCt(ct)
    # Calcul after loop
    coverage["roibase"]["mean"] = coverage["roibase"]["mean"] / len(regions)
    genes = [x for x in coverage["roigene"]["genes"].keys()]
    for gene in genes:
        coverage["roigene"]["genes"][gene]["mean"] = (
            coverage["roigene"]["genes"][gene]["mean_cumulate"]
            / coverage["roigene"]["genes"][gene]["number"]
        )
        del coverage["roigene"]["genes"][gene]["mean_cumulate"]
    bam.close()
    # Store data
    return coverage


def vcf(filename):
    vcf = pysam.VariantFile(filename)

    total = 0
    snp, indel, ins, delet = 0, 0, 0, 0
    ti, tv = 0, 0
    for variant in vcf.fetch():
        total += 1

        ref = variant.ref
        alt = variant.alts[0]

        # Simple stat
        if len(ref) != 1 or len(alt) != 1:
            if len(ref) > 1:
                delet += 1
            else:
                ins += 1
            indel += 1
        else:
            snp += 1
            # Titv
            ref = ref[0].upper()
            alt = alt[0].upper()
            if (
                (ref == "A" and alt == "G")
                or (ref == "G" and alt == "A")
                or (ref == "C" and alt == "T")
                or (ref == "T" and alt == "C")
            ):
                ti += 1
            else:
                tv += 1
    vcf.close()
    statistics = {
        "total": total,
        "snp": snp,
        "indel": indel,
        "ins": ins,
        "del": delet,
        "ti": ti,
        "tv": tv,
    }
    return statistics


def write(foutput, data):
    """
    Format output

    - files
      - fastq
        - forward_before | forward_after | reverse_before | reverse_after : filename
      - bam : filename
      - bam_format : sam or bam
      - bed : filename
      - vcf : filename
      - output : filename
    - name : name sample
    - rules
      - fastq : list of rules
      - bam : list of rules
      - vcf : list of rules
    - statistics
      - fastq
        - forward_before | forward_after | reverse_before | reverse_after
          - total : int
          - length
        - minimum : int
        - maximum : int
        - mean : int
          - quality_per_base
        - position
          - quality : int
          - quality_per_tile_per_base
        - tile
          - position
            - quality : int
              - gc (%)
      - trimming
        - total : int
        - surviving : int
        - drop : int
        - surviving_percent : int
      - bam
        - flagstat
          - reads : int
          - secondary : int
          - supp : int
          - pair_all : int
          - pair_good : int
          - read1 : int
          - read2 : int
          - sgltn : int
          - pair_map : int
          - mapped : int
          - dup : int
        - insert
          - size
        - <number> : int97 532     2 9
        - coverage
          - cutoff : list
          - on_target : int
          - roidepth
        - total : int
        - cutoff
          - <cutoff> : int
          - roibase
        - nb_base : int
        - mean : int
        - cutoff
          - <cutoff>
          - roigene
        - genes
          - <gene_name>
            - number : int
            - mean_cumulate : int
            - mean : int
      - vcf
        - total : int
        - snp : int
        - indel : int
        - ins : int
        - del : int
        - ti : int
        - tv : int
    """
    with open(foutput, "w") as fod:
        json.dump(data, fod, indent=4)
