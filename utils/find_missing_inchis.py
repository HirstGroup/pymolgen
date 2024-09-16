import argparse

parser = argparse.ArgumentParser(description='Find missing inchis')

# Required arguments
parser.add_argument('-i','--input', help='Input file name with inchis to find in aux files',required=True)
parser.add_argument('-a','--aux', nargs='+', help='Auxiliary files with inchis',required=True)
parser.add_argument('-o','--output', help='Output file containing inchis from input not found in aux files',required=True)

args = parser.parse_args()

inchis_to_check  = set()

with open(args.input) as infile:

    for n, line in enumerate(infile):

        inchi = line.strip().split()[0]

        inchis_to_check.add(inchi)

        if n % 1000000 == 0:
            print('Number of inchis loaded', n)

print('Number of inchis to check:', len(inchis_to_check))

for aux in args.aux:

    print('Auxfile', aux)

    with open(aux) as auxfile:

        for n, line in enumerate(auxfile):

            inchi = line.strip().split()[0]

            if inchi in inchis_to_check:

                inchis_to_check.remove(inchi)

            if n % 1000000 == 0:
                print('Lines checked', n)

        print('Number of inchis not found:', len(inchis_to_check))


print('Number of inchis not found:', len(inchis_to_check))

with open(args.output, 'w') as outfile:

    for inchi in inchis_to_check:

        outfile.write(f'{inchi}\n')