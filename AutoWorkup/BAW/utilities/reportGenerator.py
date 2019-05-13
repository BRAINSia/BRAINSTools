#!/bin/env python
"""
reportGenerator.py
============================
Description:
    The purpose of this is to...

Usage:
  reportGenerator.py -h | --help
  reportGenerator.py REPORT EXPERIMENT [--outdir D] [-f OUTFILE]

Arguments:
  REPORT       The report file csv with Python dictionary entry
  EXPERIMENT   The experiment directory to search for T1/T2 files

Options:
  -h, --help   Print this help and exit
  --outdir D   The directory to copy the T1/T2 tree and files, e.g. /Shared/sinapse/CACHE/report
  -f OUTFILE   The output file csv with average T1/T2 entries in Python dictionary [default: /tmp/autoworkup_report.csv]

"""


import csv
import glob
import os.path
import shutil


def only_t1t2(src, names):
    """
    This function...

    :param src:
    :param names:
    :return:
    """
    if src.endswith("TissueClassify"):
        # print "Keeping T1/T2!"
        try:
            names.remove("t1_average_BRAINSABC.nii.gz")
        except ValueError:
            pass
        try:
            names.remove("t2_average_BRAINSABC.nii.gz")
        except ValueError:
            pass
    else:
        names.remove("TissueClassify")
    # print "Ignoring these files..."
    # for name in names:
    #     print "\t" + name
    return names


def main(
    REPORT, EXPERIMENT, outdir=None, OUTFILE="/tmp/autoworkup_report.csv", **kwargs
):
    """
    This function...

    :param REPORT:
    :param EXPERIMENT:
    :param outdir:
    :param OUTFILE:
    :param kwargs:
    :return:
    """
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    if outdir is not None:
        EXPERIMENT = EXPERIMENT.rstrip(os.path.sep)
        outdir = os.path.join(outdir, os.path.basename(EXPERIMENT))
    with open(REPORT, "r") as iid, open(OUTFILE, "w") as oid:
        report_reader = csv.DictReader(iid, delimiter=",", quotechar='"')
        writer = None
        for row in report_reader:
            if writer is None:
                writer = csv.DictWriter(
                    oid, fieldnames=report_reader.fieldnames, quoting=csv.QUOTE_ALL
                )
                writer.writeheader()
            path = os.path.join(
                EXPERIMENT,
                row["project"],
                row["subject"],
                row["session"],
                "TissueClassify",
            )
            print(path)
            if outdir is not None:
                outpath = os.path.join(
                    outdir, row["project"], row["subject"], row["session"]
                )
                # assert os.path.isdir(path)
                shutil.copytree(path, outpath, ignore=only_t1t2)
                # HACK: copytree isn't copying the files, so do it manually
                try:
                    assert (
                        len(os.listdir(outpath)) > 2
                    ), "{0} doesn't have enough files".format(outpath)
                except AssertionError:
                    try:
                        fname = os.path.join(path, "t1_average_BRAINSABC.nii.gz")
                        newpath = os.path.join(outpath, "t1_average_BRAINSABC.nii.gz")
                        shutil.copyfile(fname, newpath)
                    except:
                        pass
                    try:
                        fname = os.path.join(path, "t2_average_BRAINSABC.nii.gz")
                        newpath = os.path.join(outpath, "t2_average_BRAINSABC.nii.gz")
                        shutil.copyfile(fname, newpath)
                    except:
                        pass
                # END HACK
                print(outpath)
            outdict = OrderedDict()
            olddict = eval(row["imagefiles"])
            for key in list(olddict.keys()):
                if key.startswith("T1"):
                    fname = os.path.join(path, "t1_average_BRAINSABC.nii.gz")
                    newpath = os.path.join(outpath, "t1_average_BRAINSABC.nii.gz")
                elif key.startswith("T2"):
                    fname = os.path.join(path, "t2_average_BRAINSABC.nii.gz")
                    newpath = os.path.join(outpath, "t2_average_BRAINSABC.nii.gz")
                else:
                    continue
                outdict[key] = [fname]
            row["imagefiles"] = outdict
            writer.writerow(row)


if __name__ == "__main__":
    from docopt import docopt

    args = docopt(__doc__)
    print(args)
    for key in list(args.keys()):
        if key.startswith("-"):
            value = args.pop(key)
            newkey = key.strip("-")
            args[newkey] = value
    args["OUTFILE"] = args.pop("f")
    main(**args)
