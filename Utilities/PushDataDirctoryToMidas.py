#! /usr/bin/env python
"""
midas.py
========
This script will upload a directory to the Midas server for BRAINSTools

Usage:
  midas.py [--MD5CacheDir LOCALMD5CACHEDIR] [--email EMAIL] [--password PASS | --apikey KEY] DIR
  midas.py -h | --help

Arguments:
  DIR      The directory of files you wish to upload to the BRAINSTools community

Options:
  -h, --help            Show this help and exit
  --email EMAIL         Email account for Midas login (if not given, prompted)
  --password PASS       Password associated with login (if neither this or apikey given, prompted)
  --apikey KEY          Midas key associated with 'Default' application (see http://www.kitware.com/midaswiki/index.php/Documentation/Latest/Developer/Modules/Api for more info)
  --MD5CacheDir LOCALMD5CACHEDIR   Where the .md5 files should be written
"""
import os.path
import hashlib

import sys
import shutil

def md5_for_file(f, block_size=2 ** 20):
    """Generate a hash key from a file"""
    md5 = hashlib.md5()
    while True:
        data = f.read(block_size)
        if not data:
            break
        md5.update(data)
    return md5.hexdigest()

def MakeLocalKeyFile(fullpath,LocalMD5Dir):
    """
    fullpath - is the original input file to be uploaded to midas
    LocalMD5Dir - is the local key cache of md5 files
    """

    if not os.path.exists(LocalMD5Dir):
        raise "ERROR: Output path for key files does not exists"

    f = open(fullpath)
    md5HashValue = md5_for_file(f)
    f.close()
    baseFileName = os.path.basename(fullpath)
    localKeyFileFullPath = os.path.join(LocalMD5Dir,baseFileName+".md5")
    f = open(md5FileName, 'w')
    f.write(md5HashValue)
    f.close()

#===================================================================
#===================================================================
#===================================================================


try:
    import pydas
except ImportError:
    raise ImportError("""
Have you installed pydas?  If not, you will need to either system-wide or within a virtualenv:
  pip install pydas
or
  pip install --user pydas
""")

try:
    from docopt import docopt
except ImportError:
    raise ImportError("""
Have you installed docopt?  If not, you will need to either system-wide or within a virtualenv:
  pip install docopt
or
  pip install --user docopt
""")

argv = docopt(__doc__, version='0.1')

midasURL = "http://slicer.kitware.com/midas3"
if argv['--apikey'] is None:
    sessionToken = pydas.login(url=midasURL, email=argv['--email'], password=argv['--password'])
else:
    sessionToken = pydas.login(url=midasURL, email=argv['--email'], api_key=argv['--apikey'])
communicator = pydas.session.communicator
for driver in communicator.drivers:
    if isinstance(driver, pydas.drivers.CoreDriver):
        break
community = driver.get_community_by_name(name="BRAINSTools")
for folder in driver.folder_children(token=sessionToken, folder_id=community['folder_id'])['folders']:
    if folder['name'] == u'Public':
        break
publicFolderID = folder['folder_id']
dirpath = os.path.abspath(argv['DIR'])
assert os.path.isdir(dirpath), "Not a directory"
for root, dirs, files in os.walk(dirpath):
    for filename in files:
        fullpath = os.path.join(root, filename)
        # sessionToken = pydas.verify_credentials()
        response = driver.create_item(token=sessionToken, name=filename,
                                      parent_id=publicFolderID, privacy='Public')
        uploadToken = driver.generate_upload_token(token=sessionToken, item_id=response['item_id'],
                                                   filename=filename)
        upload_response = driver.perform_upload(uploadToken, filename=fullpath)

        if argv['--MD5CacheDir'] is not None:
            MakeLocalKeyFile(fullpath,argv['--MD5CacheDir'])
