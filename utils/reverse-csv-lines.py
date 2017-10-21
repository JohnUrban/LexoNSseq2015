#!/usr/bin/env python2.7

import sys, argparse

parser = argparse.ArgumentParser(description="""

DESCRIPTION - format fasta file to have no line wrapping or desired line length.

    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--csv', '-f', '-i', '-c',
                   type= str, default=False, required=False,
                   help='''Path to input CSV file''')
args = parser.parse_args()



if args.csv == "-" or args.csv == "stdin" or not args.csv:
    args.csv = sys.stdin
else:
    args.csv = open(args.csv)


for line in args.csv:
    print (",").join( line.strip().split(",")[-1::-1] )
