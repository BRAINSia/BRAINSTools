from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from configparser import SafeConfigParser as scp
_config = scp()

from .autoworkup import *
