import filecmp
import json
import sys
import tempfile
from typing import List

import pandas as pd
from hmnqc._version import __app_name__
from hmnqc.utils import run


class TestFunctional:
    @classmethod
    def compare_json(cls, one, two, del_keys: List[str] = None) -> bool:
        with open(one) as fod:
            one_data = json.load(fod)
        with open(two) as fod:
            two_data = json.load(fod)
        if del_keys is not None:
            for key in del_keys:
                del one_data[key]
                del two_data[key]

        return one_data == two_data

    @classmethod
    def compare_xlsx(cls, one, two) -> bool:
        df1 = pd.read_excel(one)
        df2 = pd.read_excel(two)
        return df1.equals(df2)

    def test_quality(self, sample_one, one_results):
        with tempfile.NamedTemporaryFile() as fd:
            args = [__app_name__, "quality"]
            args += ["--input-sample-name", "one"]
            args += ["--output-hmnqc-json", fd.name]
            args += ["-iffr", sample_one.fq_fwd_raw]
            args += ["-ifrr", sample_one.fq_rev_raw]
            args += ["-ifft", sample_one.fq_fwd_trim]
            args += ["-ifrt", sample_one.fq_rev_trim]
            args += ["--input-sample-bam", sample_one.bam]
            args += ["--input-sample-vcf", sample_one.vcf]

            ret = run(args)
            if ret.returncode > 0:
                print(ret.stderr)
                print(ret.stdout)
                sys.exit(1)

            assert TestFunctional.compare_json(
                fd.name, one_results["quality"], del_keys=["files"]
            )

    def test_depthmin(self, sample_one, one_bed, one_results):
        with tempfile.NamedTemporaryFile(suffix=".xlsx") as fd:
            args = [__app_name__, "depthmin"]
            args += ["--input-sample-bam", sample_one.bam]
            args += ["--output-hmnqc-xlsx", fd.name]
            args += ["--input-sample-bed", one_bed]
            # --parameter-cut-off <int>

            ret = run(args)
            if ret.returncode > 0:
                print(ret.stderr)
                print(ret.stdout)
                sys.exit(1)

            assert TestFunctional.compare_xlsx(fd.name, one_results["depthmin"])

    def test_depthtarget(self, sample_one, one_bed, one_results):
        with tempfile.NamedTemporaryFile(suffix=".xlsx") as fd:
            args = [__app_name__, "depthtarget"]
            args += ["--input-sample-bam", sample_one.bam]
            args += ["--output-hmnqc-xlsx", fd.name]
            args += ["--input-sample-bed", one_bed]
            # --parameter-mode target \

            ret = run(args)
            if ret.returncode > 0:
                print(ret.stderr)
                print(ret.stdout)
                sys.exit(1)

            assert TestFunctional.compare_xlsx(fd.name, one_results["depthtarget"])

    def test_infersexe(self, sample_one, one_bed, one_results):
        with tempfile.NamedTemporaryFile(suffix=".xlsx") as fd:
            args = [__app_name__, "depthtarget"]
            args += ["--input-sample-bam", sample_one.bam]
            args += ["--output-hmnqc-xlsx", fd.name]
            args += ["--input-sample-bed", one_bed]
            # --input-sample-bed <BED file> \

            ret = run(args)
            if ret.returncode > 0:
                print(ret.stderr)
                print(ret.stdout)
                sys.exit(1)

            assert TestFunctional.compare_xlsx(fd.name, one_results["infersexe"])

    def test_extractvcf(self, sample_one, reference_vcf, one_results):
        out = "/save/ggricourt/opt/guillaume-gricourt/HmnQc/tests/one.extractvcf.xlsx"
        with tempfile.NamedTemporaryFile(suffix=".xlsx") as fd:
            args = [__app_name__, "extractvcf"]
            args += ["--input-sample-vcf", sample_one.vcf]
            args += ["--input-reference-vcf", reference_vcf]
            args += ["--output-hmnqc-xlsx", fd.name]

            ret = run(args)
            if ret.returncode > 0:
                print(ret.stderr)
                print(ret.stdout)
                sys.exit(1)

            assert TestFunctional.compare_xlsx(fd.name, one_results["extractvcf"])
