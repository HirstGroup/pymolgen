import argparse

parser = argparse.ArgumentParser(description='Find missing inchis')

# Required arguments
parser.add_argument('-i','--input', help='Input file name with inchis to find in aux files',required=True)
parser.add_argument('-a','--aux', nargs='+', help='Auxiliary files with inchis',required=True)
parser.add_argument('-o','--output', help='Output file containing inchis from input not found in aux files',required=True)

args = parser.parse_args()

inchis_to_check  = set()

with open(args.input) as infile:

    for line in infile:

        inchi = line.strip().split()[0]

        inchis_to_check.add(inchi)

print('Number of inchis to check:', len(inchis_to_check))

for aux in args.aux:

    with open(aux) as auxfile:

        for line in auxfile:

            inchi = line.strip().split()[0]

            if inchi in inchis_to_check:

                inchis_to_check.remove(inchi)


print('Number of inchis not found:', len(inchis_to_check))

with open(args.output, 'w') as outfile:

    for inchi in inchis_to_check:

        outfile.write(f'{inchi}\n')