import logging
import sys


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


def init_django(current_dir, static_dir):
    settings.configure(
        DEBUG=True,
        TEMPLATES=[
            {
                "BACKEND": "django.template.backends.django.DjangoTemplates",
                "DIRS": [current_dir],
            }
        ],
        INSTALLED_APPS=[
            "django.contrib.staticfiles",
        ],
        STATIC_URL="/static/",
        STATICFILES_DIRS=(static_dir,),
    )
    django.setup()


def read_json(path):
    with open(path) as fid:
        return json.load(fid)


def format_nb(nb, target=0, value2return=0):
    if target == 0 and value2return == 0:
        return "{:,}".format(round(float(nb))).replace(",", " ").split(".")[0]
    if value2return == 1:
        return "{:0.2f} %".format(float(target) * 100 / float(total))
    else:
        return "{} ({:0.2f} %)".format(
            formatNb(target), float(target) * 100 / float(total)
        )
