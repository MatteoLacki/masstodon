import re


def parse_threshold(x = '0.0 Da'):
    """Parse a threshold string."""
    x_type = None
    floats = re.findall(r"[-+]?\d*[\.\,]\d+|\d+", x)
    if len(floats) > 1:
        raise ValueError("The threshold contained more then one numbers: {}".format(" ".join(floats)))
    v = float(floats[0].replace(",",'.'))
    x = x.lower()
    if "da" in x:
        x_type = 'abs'
    elif "th" in x:
        x_type = 'abs'
    elif "mmu" in x:
        x_type = 'abs'
        v /= 1000.0
    elif "ppm" in x:
        x_type = 'rel'
        v *= 1e-6
    else:
        raise ValueError("The threshold did not contain its unit.")
    return v, x_type


def test_parse_threshold():
    assert parse_threshold("0.05 Da") == (0.05, 'abs')
    assert parse_threshold("0.05 Th") == (0.05, 'abs')
    assert parse_threshold("0,05 Th") == (0.05, 'abs')
    assert parse_threshold("0,05 Da") == (0.05, 'abs')
    assert parse_threshold("50.0 mmu") == (0.05, 'abs')
    assert parse_threshold("50 mmu") == (0.05, 'abs')
    assert parse_threshold("5.0 ppm") == (5.0*1e-6,  'rel')
    assert parse_threshold("5,0 ppm") == (5.0*1e-6,  'rel')
