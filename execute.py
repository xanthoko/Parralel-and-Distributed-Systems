import sys
import os
from subprocess import call

d = 0.85


def reverse_mat(path, ind_f):
    """Read and reverse the dataset.
    
    Creates an new matrix which is the same as the original
    but for every row the columns are mirrored.
    (col1 col2 -->  col2 col1)
    If the dataset is 1-indexed performs a convertion to 
    0-indexed by substracting 1 from each element.
    """
    with open(path, 'r') as f:
        a = f.readlines()

    rev_sparse = []
    content = [x.strip() for x in a]
    for line in content:
        rc = line.split(' ')[0].split('\t')
        # rc = line.split(' ')
        if ind_f:
            new_rc = [int(rc[1]) - 1, int(rc[0]) - 1]
        else:
            new_rc = [int(rc[1]), int(rc[0])]
        rev_sparse.append(new_rc)

    return rev_sparse


def write_arrays(sparse):
    """Create the sparse array and the ingoing for pagerank.c.

    ingoing: for every site, the number of sites that have a link
        towards it.
    outgoing: for every site, the number of links it contains.
    sp_array: for every site, contains the id and the factor 
        of the sites that link to it.
    Saves the sp_array as sparse.txt and the ingoing as ingoing.txt
    """
    N = max(flatten(sparse)) + 1
    sp_array = [[] for i in range(N)]
    ingoing = [0] * N
    outgoing = [0] * N

    for rs in sparse:
        outgoing[rs[1]] += 1
        ingoing[rs[0]] += 1

    for rs in sparse:
        sp_array[rs[0]].append(rs[1])
        sp_array[rs[0]].append(-d/outgoing[rs[1]])

    with open('data/sparse.txt', 'w') as f:
        for row in sp_array:
            for col in row:
                f.write(str(col) + ' ')
            f.write('\n')

    with open('data/ingoing.txt', 'w') as f:
        for row in ingoing:
            f.write(str(row) + '\n')

    return N


def flatten(seq):
    for el in seq:
        if isinstance(el, list):
            for bar in flatten(el):
                yield bar
        else:
            yield el


def main():
    if len(sys.argv) != 4:
        print('You need to specify the dataset path, the number of max iterations' \
              ' and whether the dataset is 0 or 1 indexed (0,1)!')
        exit()

    # command line arguments
    dpath = sys.argv[1]
    iterations = sys.argv[2]
    ind = int(sys.argv[3])

    try:
        edited_data = reverse_mat(dpath, ind)
        print('Matrix reversed.')
    except IOError:
        print('Wrong file path.')
        exit()

    websites = write_arrays(edited_data)
    print('Data files created.')

    # run pagerank.c <number_of_websites> <number_of_max_iterations>
    call(["./pagerank", str(websites), iterations])


if __name__ == '__main__':
    main()
