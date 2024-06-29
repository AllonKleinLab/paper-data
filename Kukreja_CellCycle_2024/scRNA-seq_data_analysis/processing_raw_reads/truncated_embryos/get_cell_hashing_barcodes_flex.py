import sys
import pickle
import itertools
from itertools import product, combinations
from collections import defaultdict

def string_hamming_distance(str1, str2):
    """
    Fast hamming distance over 2 strings known to be of same length.
    In information theory, the Hamming distance between two strings of equal 
    length is the number of positions at which the corresponding symbols 
    are different.
    eg "karolin" and "kathrin" is 3.
    """
    return sum(itertools.imap(operator.ne, str1, str2))

___tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
def rev_comp(seq):
    return ''.join(___tbl[s] for s in seq[::-1])

def seq_neighborhood(seq, n_subs=1):
    """
    Given a sequence, yield all sequences within n_subs substitutions of 
    that sequence by looping through each combination of base pairs within
    each combination of positions.
    """
    for positions in combinations(range(len(seq)), n_subs):
    # yields all unique combinations of indices for n_subs mutations
        for subs in product(*("ATGCN",)*n_subs):
        # yields all combinations of possible nucleotides for strings of length
        # n_subs
            seq_copy = list(seq)
            for p, s in zip(positions, subs):
                seq_copy[p] = s
            yield ''.join(seq_copy)

def build_barcode_neighborhoods(barcode_file, expect_reverse_complement=False, max_hamming_dist=2):
    """
    Given a set of barcodes, produce sequences which can unambiguously be
    mapped to these barcodes, within `max_hamming_dist` substitutions. If a sequence maps to 
    multiple barcodes, get rid of it. However, if a sequences maps to a bc1 with 
    1change and another with 2changes, keep the 1change mapping.
    """

    # contains all mutants that map uniquely to a barcode
    clean_mapping = dict()

    # contain single or double mutants 
    mappings = {i: defaultdict(set) for i in range(1, max_hamming_dist+1)}
    
    #Build the full neighborhood and iterate through barcodes
    with open(barcode_file, 'r') as f:
        # iterate through each barcode (rstrip cleans string of whitespace)
        for line in f:
            barcode = line.rstrip()
            if expect_reverse_complement:
                barcode = rev_comp(line.rstrip())

            # each barcode obviously maps to itself uniquely
            clean_mapping[barcode] = barcode

            # for each possible mutated form of a given barcode, either add
            # the origin barcode into the set corresponding to that mutant or 
            # create a new entry for a mutant not already in mapping1
            # eg: barcodes CATG and CCTG would be in the set for mutant CTTG
            # but only barcode CATG could generate mutant CANG
            for dist in mappings:
                for n in seq_neighborhood(barcode, dist):
                    mappings[dist][n].add(barcode)
            
    # take all single-mutants and find those that could only have come from one
    # specific barcode
    for dist in range(1, max_hamming_dist+1):
        for k, v in mappings[dist].items():
            if k not in clean_mapping:
                if len(v) == 1:
                    clean_mapping[k] = list(v)[0]
    
    del mappings
    return clean_mapping

barcode_file = sys.argv[1]
umi_filename = sys.argv[2]
read_filename = sys.argv[3]
pickle_filename = sys.argv[4]
if len(sys.argv) > 5:
    max_hamming_dist = int(sys.argv[5])
else:
    max_hamming_dist = 2

barcode_neighborhoods = build_barcode_neighborhoods(barcode_file, expect_reverse_complement=False, max_hamming_dist=max_hamming_dist)
hashtags = list(set(barcode_neighborhoods.values()))
barcode_length = list(set([len(bc) for bc in hashtags]))
if len(barcode_length) > 1:
    print('ERROR: multiple cell hashing barcode lengths detected')
    sys.exit()
barcode_length = barcode_length[0]

d = {}

for iL,l in enumerate(sys.stdin):
    if iL % 1e6 == 0:
        print(iL)
    if iL % 4 == 0:
        cols = l[1:].strip('\n').split(':')
        cb = cols[0]
        umi = cols[1]
        if cb not in d:
            d[cb] = {'no match': 0}
    elif iL % 4 == 1:
        l = l.strip('\n')
        found_barcode_match = False
        seq_start = 0
        while not found_barcode_match and seq_start < (len(l) - barcode_length):
            seq = l[seq_start:(seq_start + barcode_length)]
            if seq in barcode_neighborhoods:
                found_barcode_match = True
                barcode = barcode_neighborhoods[seq]
                if barcode not in d[cb]:
                    d[cb][barcode] = []
                d[cb][barcode].append((seq_start, umi))
            seq_start += 1
        if not found_barcode_match:
            d[cb]['no match'] += 1

# build UMI and Read count matrices
umi_file = open(umi_filename, 'w')
read_file = open(read_filename, 'w')
header = 'Cell_Barcode,' + ','.join(hashtags) + '\n'
read_file.write(header)
umi_file.write(header)
for cb,hits in d.items():
    umi_counts = {h:0 for h in hashtags}
    read_counts = {h:0 for h in hashtags}
    for htag,umis in hits.items():
        if htag != 'no match':
            umi_counts[htag] = len(set([u[1] for u in umis]))
            read_counts[htag] = len(umis)

    outline = cb + ','  + ','.join([str(umi_counts[htag]) for htag in hashtags]) + '\n'
    umi_file.write(outline)

    outline = cb + ',' + ','.join([str(read_counts[htag]) for htag in hashtags]) + '\n'
    read_file.write(outline)

umi_file.close()
read_file.close()

pickle.dump(d, open(pickle_filename, 'wb'), -1)