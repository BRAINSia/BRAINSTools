"""
bidsUtils_test.py
============================
Description:
    This file is used to test the bidsUtils formatter

Author:

Usage: N/A not a command line script (no bang above)

"""
# I don't know if pytest is a dependency of this project, so I am avoiding adding it for the moments

from bidsUtils import BAWBIDSFormatter

if __name__ == "__main__":
    try:
        subject_info = {
            "subject": "0123456",
            "session": "43210",
            "space": "LPS",
            "run": "05",
            "favColor": "red",
            "spacing": "iso2mm",
        }
        bids_formatter = BAWBIDSFormatter()
        assert (
            "/Shared/sinapse/out/sub-0123456_ses-43210_favColor-red_run-05_space-LPS_spacing-iso2mm.nii.gz"
            == bids_formatter.get_bids_name(
                subject_data=subject_info,
                full_path="/Shared/sinapse/out/",
                ext="nii.gz",
            )
        )
        print("Success")
    except Exception:
        print("Failure")

    try:
        no_sub = {
            "session": "43210",
            "space": "LPS",
            "run": "05",
            "favColor": "red",
            "spacing": "iso2mm",
        }
        bids_formatter = BAWBIDSFormatter()
        bids_formatter.get_bids_name(subject_data=no_sub)
        print("Failure")
    except KeyError:
        print("Success")

    try:
        no_ses = {
            "subject": "0123456",
            "space": "LPS",
            "run": "05",
            "favColor": "red",
            "spacing": "iso2mm",
        }
        bids_formatter = BAWBIDSFormatter()
        bids_formatter.get_bids_name(subject_data=no_ses)
        print("Failure")
    except KeyError:
        print("Success")
