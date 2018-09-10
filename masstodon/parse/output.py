from masstodon.parse.path import parse_path

def add_backslash(p):
    if p[-1] != '/':
        p += '/'
    return p

def parse_output(args):
    output_path, file_name, _ = parse_path(args.spectrum)
    output_path = add_backslash(output_path)
    if args.output_path:
        output_path = args.output_path
    else:
        output_path += 'output'
    output_path = add_backslash(output_path)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    return output_path
