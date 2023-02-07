import os
from typing import Dict, Tuple

import pytest

cur_dir = os.path.abspath(os.path.dirname(__file__))
data_dir = os.path.join(cur_dir, "dataset")


@pytest.fixture(scope="session")
def data_directory():
    return data_dir


@pytest.fixture(scope="session")
def one_fq_raw(data_directory) -> Tuple[str, str]:
    return (
        os.path.join(data_directory, "one.R1.fastq"),
        os.path.join(data_directory, "one.R2.fastq"),
    )


@pytest.fixture(scope="session")
def one_fq_trim(data_directory) -> Tuple[str, str]:
    return (
        os.path.join(data_directory, "one.R1.trim.fastq"),
        os.path.join(data_directory, "one.R2.trim.fastq"),
    )


@pytest.fixture(scope="session")
def one_bam(data_directory) -> str:
    return os.path.join(data_directory, "one.bam")


@pytest.fixture(scope="session")
def one_vcf(data_directory) -> str:
    return os.path.join(data_directory, "one.vcf")


@pytest.fixture(scope="session")
def one_results(data_directory) -> Dict[str, str]:
    return dict(quality=os.path.join(data_directory, "one.quality.json"))
