import argparse
import os
import sys

parser = argparse.ArgumentParser(description='This script submits calculator and rules scripts to HPC')
parser.add_argument('-i','--input', help='Input file name', required=True)
parser.add_argument('-l','--lines', type=int, help='Number of lines to split input file in', required=True)
parser.add_argument('--hpc', default='augusta', help='Name of HPC to look for sbatch header file', required=False)
args = parser.parse_args()

infile = open('%s' %args.input)

print( args.input, args.output)

lines = args.lines

n = 0
m = 0

for line in infile:
	if n%lines == 0: 
		outfile = open('%s_%s.txt' %(args.output, m), 'w' )
		m += 1
	outfile.write(line)
	n +=1

outfile.close()

outfile = open('array.sh', 'w')

last = n

outfile.write(f'''#!/bin/bash
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=16gb
#SBATCH --array=0-{last}

. ~/.bashrc

j=${{SLURM_ARRAY_TASK_ID}}
i={depth}

module load active_search/active_search
module load ruby-uon/gcc6.3.0/2.6.4

set -e

python ~/pymolgen/calculator.py -i whz-threshold8-depth${{i}}-unique_$j.txt -o whz-threshold8-depth${{i}}-unique_${{j}}-calculator.txt --mw_threshold 500.0 -bp

python ~/pymolgen/rules.py -i whz-threshold8-depth${{i}}-unique_${{j}}-calculator.txt -o whz-threshold8-depth${{i}}-unique_${{j}}-calculator-rules.txt -r whz-threshold8-depth${{i}}-unique_${{j}}-calculator-rules.smi --all

''')

outfile.close()

os.system('sbatch array.sh')
