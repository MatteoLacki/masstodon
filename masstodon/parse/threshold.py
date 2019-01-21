import re

pattern = re.compile(r"(\d+([.,]\d*)?|([.,]\d*))([a-zA-Z]+)")


def parse(x = '0.0Da'):
    """Parse a threshold string."""
    g = re.match(pattern, x)
    unit = g[4].lower()
    assert unit in ('da', 'th', 'mmu', 'ppm'), "Wrong or missing unit."
    v = float(g[1].replace(',','.'))
    x_type = 'abs'
    if unit == 'mmu':
        v /= 1000.
    if unit == 'ppm':
        x_type = 'rel'
        v *= 1e-6
    return v, x_type


def test_parse():
    assert parse("0.05Da") == (0.05, 'abs')
    assert parse("0.05Th") == (0.05, 'abs')
    assert parse("0,05Th") == (0.05, 'abs')
    assert parse("0,05Da") == (0.05, 'abs')
    assert parse("50.0mmu") == (0.05, 'abs')
    assert parse("50mmu") == (0.05, 'abs')
    assert parse("5.0ppm") == (5.0*1e-6, 'rel')
    assert parse("5,0ppm") == (5.0*1e-6, 'rel')
    assert parse(",2ppm") == (.2*1e-6, 'rel')
    assert parse(".3ppm") == (.3*1e-6, 'rel')