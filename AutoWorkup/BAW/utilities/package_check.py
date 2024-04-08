"""
package_check.py
============================
Description:
    The purpose of this is to..

Usage:

"""

try:
    from nipype.utils.misc import package_check
except ImportError:
    import sys

    print(sys.path)
    raise ImportError(
        "Cannot import nipype.utils.misc.package_check(). \
        Verify that the sys.path includes the correct Nipype path.  \
        NOTE: this call should be made AFTER configuration file environment parameters are set!"
    )


def verify_packages(application="AutoWorkup"):
    """
    This function...

    :param application:
    :return:
    """
    package_version = [
        ("nipype", "0.9"),
        ("numpy", "1.8"),
        ("scipy", "0.13"),
        ("networkx", "1.8"),
        # ('IPython', '1.2'),
        # ('SimpleITK', '0.7')
    ]
    for item in package_version:
        package_check(*item, app=application)
