import json


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
