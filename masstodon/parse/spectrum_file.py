from .mzml import read_mzml
from .txt import read_txt
from .npy import read_npy



readers = {
    'txt': read_txt,
    'mzml': read_mzml,
    'mzxml': read_mzml,
    'npy': read_npy
}

# TODO scan number need to be included in the mzml
def parse(path, scan_no=0):
    """Parse a path to spectrum."""
    p, ext = path.split('.')
    ext = ext.lower()
    try:
        mz, i = readers[ext](path)
    except KeyError as e:
        print("Supply a file with a proper extension: txt/mzml/mzxml")
        raise e
    return mz, i