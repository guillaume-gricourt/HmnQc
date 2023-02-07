import json
import sys
import tempfile
from typing import List

from hmnqc._version import __app_name__
from hmnqc.utils import run


class TestFunctional:
    @classmethod
    def compare_json(cls, one, two, del_keys: List[str] = None):
        with open(one) as fod:
            one_data = json.load(fod)
        with open(two) as fod:
            two_data = json.load(fod)
        if del_keys is not None:
            for key in del_keys:
                del one_data[key]
                del two_data[key]

        return one_data == two_data

    def test_quality(self, one_fq_raw, one_fq_trim, one_bam, one_vcf, one_results):
        with tempfile.NamedTemporaryFile() as fd:
            args = [__app_name__, "quality"]
            args += ["--input-sample-name", "one"]
            args += ["--output-hmnqc-json", fd.name]
            args += ["-iffr", one_fq_raw[0]]
            args += ["-ifrr", one_fq_raw[1]]
            args += ["-ifft", one_fq_trim[0]]
            args += ["-ifrt", one_fq_trim[1]]
            args += ["--input-sample-bam", one_bam]
            args += ["--input-sample-vcf", one_vcf]

            ret = run(args)
            if ret.returncode > 0:
                print(ret.stderr)
                print(ret.stdout)
                sys.exit(1)

            assert TestFunctional.compare_json(
                fd.name, one_results["quality"], del_keys=["files"]
            )
