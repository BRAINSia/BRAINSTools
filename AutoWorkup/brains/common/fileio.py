def check_file(path):
    import os.path
    full = os.path.abspath(path)
    if os.path.exists(full):
        return full
    return None

def parseLabelsFile():
    import os.path
    from ..config import _config
    build_tree = _config.get('Resources', 'build_directory')
    filename = _config.get('Resources', 'label_template')
    fullname = check_file(os.path.join(build_tree, filename))
    labelDict = {}
    with open(fullname, 'r') as fid:
        for line in fid.readlines():
            if line[0] != "#":
                parts = line.split(" ")
                labelDict[int(parts[0])] = parts[1]
    return labelDict

