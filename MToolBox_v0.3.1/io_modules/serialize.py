# encoding=utf8

from classifier import datatypes
# if '..' not in sys.path:
#     sys.path.append('..')

import cPickle

def save_tree_to_file(tree_class, fname):
    f = open(fname, 'wb')
    cPickle.dump(tree_class, f, cPickle.HIGHEST_PROTOCOL)
    
def load_tree_from_file(fname):
    f = open(fname, 'rb')
    return cPickle.load(f)

def load_fasta_file(fname):
    f = open(fname, 'r')
    name = f.readline()[1:-1].strip()
    seq = []
    for line in f: seq.append(line.strip())
    return datatypes.Sequence(name, ''.join(seq).upper())

def write_fasta_file(fname, seqs):
    f = open(fname, 'w')
    for name, seq in seqs:
        f.write('>'+name+'\n')
        f.write(seq+'\n')
    f.close()

