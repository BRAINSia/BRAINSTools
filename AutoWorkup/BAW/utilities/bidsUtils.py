"""
bidsUtils.py
============================
Description:
    This file is used for formatting pipline data into the BIDS format

Author:
    Alex Powers

Usage: N/A not a command line script (no bang above)

"""
from os.path import join


class BAWBIDSFormatter(object):
    def __init__(self):
        self.sub_kw = "subject"
        self.ses_kw = "session"
        self.REQUIRED_FIELDS = [self.ses_kw, self.sub_kw]

    def get_bids_name(
        self, subject_data: dict, full_path: str = None, ext: str = None
    ) -> str:
        """
        :param subject_data: a dictionary of information about the subject including subject id and session id
        :param full_path: a string representing the path to join the bids name to
        :param ext: an optional file extension parameter
        :return: a formatted string containing all of the subject information as a file name
        """
        # input validation (subject and session must be given)
        if self.sub_kw not in subject_data.keys():
            raise KeyError(
                "Subject must have a key of '{kw}' in order to generate a meaningful filename.".format(
                    kw=self.sub_kw
                )
            )
        if self.ses_kw not in subject_data.keys():
            raise KeyError(
                "Session must have a key of '{kw}' in order to generate a meaningful filename.".format(
                    kw=self.ses_kw
                )
            )

        # build the bids name
        bids_name = "sub-{subject}_ses-{session}".format(
            subject=subject_data[self.sub_kw], session=subject_data[self.ses_kw]
        )

        for k in sorted(subject_data.keys()):
            if k not in self.REQUIRED_FIELDS:
                bids_name = "{base}_{key}-{value}".format(
                    base=bids_name, key=k, value=subject_data[k]
                )

        # add the extension if defined
        if ext:
            bids_name = "{base}.{ext}".format(base=bids_name, ext=ext)

        # add full path if defined
        if full_path:
            bids_name = join(full_path, bids_name)

        return bids_name
