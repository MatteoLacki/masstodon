def parse(p):
    """Parse app-wise parameters.
    
    Arguments
    =========
    p : argparse.ArgumentParser
        An instance of the argument parser.

    Returns
    =======
    argparse.ArgumentParser
        The same parser with more details about what to parse.
    """
    p.add_argument("spectrum",
            help="Path to the spectrum file: parseable extensions include '.txt' for ASCII, '.mzxml', or '.mzml' (case insensitive).")
    return p
