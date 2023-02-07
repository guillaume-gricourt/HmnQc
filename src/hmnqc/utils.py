import json
import logging
import subprocess
from typing import List


def run(args: List[str], show_output: bool = True) -> subprocess.CompletedProcess:
    """Run a command line.
    Parameters
    ----------
    args: List[str]
    A list of argument
    show_output: bool (default: True)
    Output command line
    Return
    ------
    subprocess.CompletedProcess
    Return result obtained with subprocess
    """
    ret = subprocess.run(args, capture_output=True, encoding="utf8")
    if show_output and ret.stdout is not None:
        logging.info(ret.stdout)
    if show_output and ret.stderr is not None:
        logging.warning(ret.stderr)
    return ret


def abort(parser, msg=""):
    """Abort the program"""
    parser.error(msg)


def calculate_percentile(dico, target):
    total = 0
    ct = (target * (sum(dico.values()) + 1)) / 100
    for key, value in sorted(dico.items(), key=lambda x: x[0]):
        total += value
        if total >= ct:
            return key
    return -1


def read_json(path):
    with open(path) as fid:
        return json.load(fid)


def isNaN(num):
    return num != num
