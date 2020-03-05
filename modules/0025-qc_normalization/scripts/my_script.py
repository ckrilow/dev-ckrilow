#!/usr/bin/env python

__author__ = "Forename Surname"
__date__ = "date"
__version__ = '0.1.0'

"""
my_script.py - skeleton for commandline script
"""

import argparse
import gzip
import sys

def run(f):
	"""
	This method will run on file handle (f)
	"""
	pass

def main():
    parser = argparse.ArgumentParser(
        description=
            """
            This is a generic skeleton for a Python script that is
            intended to be run from the command line, complete with
            rudimentary documentation and command line argument support.
            """)

    parser.add_argument('-v', '--version', action='version',
        version='%(prog)s {version}'.format(version=__version__))

    parser.add_argument('-f', '--file', action='store',
        dest='f',
        required=True,
        help=('tsv file with the following columns: '
            'file<tab>name. '
        )
    )

    options = parser.parse_args()

	# Example code that will read from stdin, gzip, or norm file
	if options.f == "stdin" or options.f == "-": # read from stdin
		run(sys.stdin)
	elif options.f.split(".")[-1] == "gz": # read as gzip file
		with gzip.open(options.f, 'rb') as f:
			run(f)
	else:
		with open(options.f, 'rb') as f: # read as normal file
			run(f)

if __name__ == '__main__':
	main()
