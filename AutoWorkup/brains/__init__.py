__version__ = 'beta'
#import config
#import metrics

def check_file(filename):
    import os.path
    fullName = os.path.abspath(filename)
    assert os.path.exists(fullName), "File %s cannot be found!" % fullName
    return fullName
