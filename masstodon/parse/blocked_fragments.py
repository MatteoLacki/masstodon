"""Parser for blocked fragments."""
import re


pattern = re.compile('([a-z])([0-9])*')


def parse(blocked_fragments):
    return  set([e+c for e, c in re.findall(pattern, blocked_fragments)])
