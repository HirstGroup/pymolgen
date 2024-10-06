import argparse
import os

parser = argparse.ArgumentParser(description='Count number of fragments in inchi file by taking information from txt file')

# required arguments
parser.add_argument('-a','--aux', help='Auxiliary txt file', required=True)
parser.add_argument('-i','--input', help='Input inchi file', required=True)
parser.add_argument('-o','--output', help='Output inchi file', required=True)

# optional arguments
parser.add_argument('--delete', action='store_true', help='Delete input file after succesful run', required=False)

args = parser.parse_args()

with open(args.aux, "rb") as f:
    aux_lines = sum(1 for _ in f)

with open(args.input, "rb") as f:
    input_lines = sum(1 for _ in f)

assert aux_lines == input_lines

assert args.aux != args.output
assert args.input != args.output

auxfile = open(args.aux)
infile = open(args.input)
outfile = open(args.output, 'w')

def count_depth(line):
	depth = len(line.strip().split(':')[0].split('-')) - 1
	return depth


for auxline, inline in zip(auxfile, infile):
	depth = count_depth(auxline)
	inchi = inline.strip().split()[0]
	outfile.write(f'{inchi} 1 d{depth}\n')

auxfile.close()
infile.close()
outfile.close()

with open(args.input, "rb") as f:
    output_lines = sum(1 for _ in f)

assert output_lines == input_lines

if args.delete:
	os.remove(args.input)



