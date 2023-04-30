import os
import sys

import argparse


# parse command line args
my_parser = argparse.ArgumentParser(description="RGBW Color Space Converter Playground")

# Add the arguments
my_parser.add_argument(
    "-i",
    "--input-file",
    action="store",
    help="The input csv or tsv file.  Type is detected from the file extension.  The first line must be column names, and the firt ror of the first column must be called 'Sample'.  Brief col names will render more nicely."
    default=False,
)


my_parser.add_argument(
    "-o",
    "--outfile_name",
    action="store",
    default="somefile_mqc.tsv",
    help="This is the output file, which needs to end in a '_mqc.tsv' or csv for it to be detected by MQC."
)


args = my_parser.parse_args()

suffix=''
if args.input_file.endswith('.tsv'):
    suffix = 'tsv'
elif args.input_file.endswith('csv'):
    suffix = 'csv'
else:
    raise Exception(f"Input file extensionn does not match tsv or csv")

print( "Adding headers to  {args.input_file} and creating new file {args.out_file}" )
