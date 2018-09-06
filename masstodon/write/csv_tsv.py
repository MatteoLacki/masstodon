import csv
import os

from masstodon.parse.path import parse_path


def write_rows(flow_of_rows, path):
    """Write a stream of rows into a csv or tsv file, depending on file path."""
    file_path, file_name, file_ext = parse_path(path)
    assert file_ext in ('.csv', '.tsv'), "Writing only to csv or tsv."
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    delimiter = ',' if file_ext == '.csv' else '\t'
    with open(path, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=delimiter)
        for row in flow_of_rows:
            writer.writerow(row)
