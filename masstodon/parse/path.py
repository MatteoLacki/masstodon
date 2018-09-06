import os

def parse_path(path):
    """Parse a file path.

    Parameters
    ----------
    path : str
        Any path.
    Returns
    -------
    out : tuple
        Path of file, name of file, and file's extension.

    """
    file_path, file_ext = os.path.splitext(path)
    file_name = file_path.split('/')[-1]
    file_path = "/".join(file_path.split('/')[:-1]) + '/'
    return file_path, file_name, file_ext
