# \author Hans J. Johnson
# This little script computes various measure
# that help ideentify search reagions to be used for finding eyes
# relative to the center of head mass
import pandas as pd
import numpy as np


def read_fiducial_points(fcsv_file: str, fiducials: list) -> dict:
    col_names = "id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID"
    with open(fcsv_file) as fid:
        df = pd.read_csv(fid, sep=",", comment="#", names=col_names.split(","))

    fid_pts = {}
    for fid in fiducials:
        fid_pts[fid] = df[df["label"] == fid].iloc[0][1:4].values
        # flip RAS fiducials to be LPS
        fid_pts[fid][0] = -1 * fid_pts[fid][0]
        fid_pts[fid][1] = -1 * fid_pts[fid][1]
    return fid_pts


if __name__ == "__main__":
    import sys
    import math

    fiducials_fn = sys.argv[1]

    lmks = read_fiducial_points(fiducials_fn, ["LE", "RE", "CM"])

    print(f"{lmks}")
    diff_vector_re = lmks["RE"] - lmks["CM"]
    diff_ap_re = lmks["RE"][1] - lmks["CM"][1]
    distance_to_re = np.linalg.norm(diff_vector_re)
    diff_vector_le = lmks["LE"] - lmks["CM"]
    diff_ap_le = lmks["RE"][1] - lmks["CM"][1]
    distance_to_le = np.linalg.norm(diff_vector_le)

    angle_si_re = np.degrees(
        -1.0 * np.sign(diff_ap_re) * math.acos(-diff_ap_re / distance_to_re)
    )
    angle_si_le = np.degrees(
        -1.0 * np.sign(diff_ap_le) * math.acos(-diff_ap_le / distance_to_le)
    )

    print(
        f"measures,{distance_to_re},{distance_to_le},{diff_ap_re},{diff_ap_le},{angle_si_re},{angle_si_le}"
    )
