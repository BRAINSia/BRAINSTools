"""
autoworkup.py
=================
Description:

Author:

Usage:

"""
from future import standard_library

standard_library.install_aliases()
from . import _config

valid_schemes = ["BRAINS", "Nipype"]


def write_configuration(filename="example.config"):
    """
    This function...

    :param filename:
    """
    config = ConfigParser.SafeConfigParser()
    config.add_section("Results")
    config.set("Results", "directory", "/full/path/to/experiment/directory")
    config.set("Results", "segmentations", "SingleRFSegmentations")
    config.set("Results", "posteriors", "TissueClassify")
    config.set("Results", "scheme", "BRAINS")
    with open(filename, "wb") as configfile:
        config.write(configfile)


def load_configuration(configFile="/dev/null"):
    """
    This function...

    :param configFile:
    """
    _config.read(configFile)


# _config = load_configuration()
