import re

def parse_blocked_fragments(blocked_fragments):
    pattern = re.compile('([a-z])([0-9]*)')
    return  set([e+c for e, c in re.findall(pattern, blocked_fragments)])
