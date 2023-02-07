import os
from collections import namedtuple
from typing import Dict, Tuple

import pytest

cur_dir = os.path.abspath(os.path.dirname(__file__))
data_dir = os.path.join(cur_dir, "dataset")

Sample = namedtuple(
    "Sample", ["fq_fwd_raw", "fq_rev_raw", "fq_fwd_trim", "fq_rev_trim", "bam", "vcf"]
)


@pytest.fixture(scope="session")
def data_directory():
    return data_dir


@pytest.fixture(scope="session")
def sample_one(data_directory) -> "Sample":
    return Sample(
        fq_fwd_raw=os.path.join(data_directory, "one.R1.fastq"),
        fq_rev_raw=os.path.join(data_directory, "one.R2.fastq"),
        fq_fwd_trim=os.path.join(data_directory, "one.R1.trim.fastq"),
        fq_rev_trim=os.path.join(data_directory, "one.R2.trim.fastq"),
        bam=os.path.join(data_directory, "one.bam"),
        vcf=os.path.join(data_directory, "one.vcf"),
    )


@pytest.fixture(scope="session")
def one_bed(data_directory) -> str:
    return os.path.join(data_directory, "one.bed")


@pytest.fixture(scope="session")
def reference_vcf(data_directory) -> str:
    return os.path.join(data_directory, "reference.vcf")


@pytest.fixture(scope="session")
def one_results(data_directory) -> Dict[str, str]:
    return dict(
        quality=os.path.join(data_directory, "one.quality.json"),
        depthmin=os.path.join(data_directory, "one.depthmin.xlsx"),
        depthtarget=os.path.join(data_directory, "one.depthtarget.xlsx"),
        infersexe=os.path.join(data_directory, "one.infersexe.xlsx"),
        extractvcf=os.path.join(data_directory, "one.extractvcf.xlsx"),
    )
