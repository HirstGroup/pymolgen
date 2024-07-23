import argparse
import sys

parser = argparse.ArgumentParser(description='Select structures from smiles or inchi file')
parser.add_argument('-i','--input', help='Input file',required=True)
parser.add_argument('--output1', help='Output file with structures in selection',required=True)
parser.add_argument('--output2', help='Output file with structures not in selection',required=True)
parser.add_argument('--sel', nargs='+', type=int, help='List of structures starting from zero, space-separated',required=True)

args = parser.parse_args()

if len(set([args.input, args.output1, args.output2])) != 3:
    sys.exit('Same input, output1, output2 files')

infile = open(args.input)

outfile1 = open(args.output1, 'w')
outfile2 = open(args.output2, 'w')

n = 0
for line in infile:
    n += 1
    
    if n in args.sel:
        outfile1.write(line)
    else:
        outfile2.write(line)

