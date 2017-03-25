from __future__ import print_function

import argparse
import difflib
import json
import sys

from nbformat import reads


def _parse_args():
    parser = argparse.ArgumentParser(
        description='Compare two very similar .ipynb files.')
    parser.add_argument('f1', metavar='FILE1', type=str,
                        help='The first file to compare.')
    parser.add_argument('f2', metavar='FILE2', type=str,
                        help='The second file to compare.')
    return parser.parse_args()


def parse_ipynb(fname):
    """Use nbformat.reads to parse a .ipynb file."""
    # Figure out the version number.
    with open(fname) as fh:
        contents = fh.read()
    version = json.loads(contents)['nbformat']

    # Parse and return.
    parsed = reads(contents, version)
    return parsed


if __name__ == '__main__':
    # Parse commandline arguments.
    args = _parse_args()

    # Parse the two .ipynb files.
    file1 = parse_ipynb(args.f1)
    file2 = parse_ipynb(args.f2)

    # Compare the number of cells.
    if len(file1.cells) != len(file2.cells):
        print('These notebooks have different numbers of cells.')
        print('{:s} has {:d} cells'.format(args.f1, len(file1.cells)))
        print('{:s} has {:d} cells'.format(args.f2, len(file2.cells)))

    # Compute cell diffs.
    diffs = difflib.SequenceMatcher(a=[c.source for c in file1.cells],
                                    b=[c.source for c in file2.cells])

    # Iterate over cell diffs.
    for tag, i1, i2, j1, j2 in diffs.get_opcodes():
        if tag == 'equal':
            pass
        elif tag == 'replace':
            assert i2-i1 == j2-j1
            # Iterate over range of cells in this opcode.
            for i in range(i2 - i1):
                # Compute and print diffs of the cell contents.
                print('\n')
                small_diff = difflib.unified_diff(
                    file1.cells[i1+i].source.split('\n'),
                    file2.cells[j1+i].source.split('\n'),
                    fromfile='Cell {:d} of {:s}'.format(i1+i, args.f1),
                    tofile='Cell {:d} of {:s}'.format(j1+i, args.f2),
                    lineterm='')
                for line in small_diff:
                    print(line)
        else:
            raise NotImplemented('IDK how to handle opcode "{:s}"'.format(tag))
