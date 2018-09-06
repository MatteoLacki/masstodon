def repr_long_list(L):
    if len(L) > 2:
        return list(L[0:1]).__repr__()[:-1] + ", ..., " + str(L[-1]) + "]"
    else:
        return list(L).__repr__()

def add_backslash(p):
    if p[-1] != '/':
        p += '/'
    return p
