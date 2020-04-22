## -- https://research-git.uiowa.edu/SINAPSE/SINAPSE/issues/30
## partially complete fcsv landmark detector for slicer fiducials.
## placing here for future fine tuning.

import pandas as pd


def read_fiducial_points(fcsv_file: str, fiducials: list) -> dict:
    col_names = "id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID"
    with open(fcsv_file, "r") as fid:
        df = pd.read_csv(fid, sep=",", comment="#", names=col_names.split(","))

    fid_pts = {}
    for fid in fiducials:
        fid_pts[fid] = df[df["label"] == fid].iloc[0][1:4].values
        # flip RAS fiducials to be LPS
        fid_pts[fid][0] = -1 * fid_pts[fid][0]
        fid_pts[fid][1] = -1 * fid_pts[fid][1]
    return fid_pts
