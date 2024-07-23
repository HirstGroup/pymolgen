import argparse


parser = argparse.ArgumentParser(description='Find list of inchis in another list of inchis')

parser.add_argument('-i','--input', help='Input file containing inchis to look for', required=True)
parser.add_argument('-a','--aux', help='Aux file with large set of inchis to check if those in input are present', required=True)
parser.add_argument('-o1','--output1', help='Output file with found inchis', required=True)
parser.add_argument('-o2','--output2', help='Output file with not found inchis', required=True)

args = parser.parse_args()

inchi_list = []

if len(set([args.input, args.output1, args.output2, args.aux])) != 4:
    sys.exit('Same input, output1, output2, aux')

with open(args.input) as infile:
    for line in infile:
        inchi = line.strip().split()[0]
        
        inchi_list.append(inchi)

inchi_set = set(inchi_list)

n = len(inchi_set)
print('N =', n)

made = 0
made_list = []

with open(args.aux) as auxfile, open(args.output1, 'w') as outfile1:

    for line in auxfile:
    
        inchi = line.strip().split()[0]
        
        if inchi in inchi_set:
            made += 1
            made_list.append(inchi)
            outfile1.write(f'{inchi}\n')
           
        if made == n:
            sys.exit('All input molecules made')
            
not_made_list = [i for i in inchi_set if i not in made_list]

print('Not made n =', len(not_made_list) )

with open(args.output2, 'w') as outfile2:

    for i in not_made_list:
        outfile2.write(f'{i}\n')
            
        
        

