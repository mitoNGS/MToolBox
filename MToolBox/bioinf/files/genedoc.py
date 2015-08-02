# coding=utf8
"""
Formato genedoc (msf), supportata lettura

@note: info come checksum e altro viene perso, mantenuti solo il nome e la seq
"""
import os

FILE_GENEDOC = 'genedoc'
FILETYPE = 0
ERR_FILE_NOT_EXIST = "file %s doesn't exists"

def check_genedoc(fname):
    """
    Rudimentale controllo, vede se i primi 3 byte sono GDC
    
    @param fname:  il nome del file
    
    @raise ValueError: nel caso il file non esista
    """
    if not os.path.isfile(fname):
        raise ValueError(ERR_FILE_NOT_EXIST % fname)
    f = open(fname, 'U')
    if f.read(3) == 'GDC':
        f.close()
        return True
    else:
        f.close()
        return False

def load_genedoc(fname, addFunc):
    """
    Load Alignment from GeneDoc File
    
    @type fname: string
    @param fname: nome del file
    @type addFunc: function
    @param addFunc: puntatore alla funzione di aggiunta
    
    @raise ValueError: nel caso il file non esista
    """
    if not os.path.isfile(fname):
        raise ValueError(ERR_FILE_NOT_EXIST % fname)
    head = []
    seqs = {}
    f = open(fname,'U')
    head_gdc = head_info = head_seqs = False
    for line in f:
        # Set block being analysed
        if line[0:3] == 'GDC': # Don't know what this data represent, will add code if necessary
            head_gdc = True
        elif line.lstrip()[:5] == 'Name:': # don't know if 1st space is necessary or not, so...
            head_info = True; head_gdc = False
        elif line[0:2] == '//':
            head_seqs = True; head_info = False
        # Add info according to block
        if line == '\n': # skip blank lines
            pass
        elif head_gdc:
            head.append(line) # just add this data for later
        elif head_info:
            tmp = line.split()
            if tmp[0] == 'Name:': # may be first line after 'GDC'
                seqs[tmp[1].lower()] = []
        elif head_seqs:
            #  GeneDoc saves gaps as '.'
            tmp = line.split() # split eliminates whitespace
            if not tmp[0].isdigit() and len(tmp) > 1:
                seqs[tmp[0].lower()].extend([block.replace('.','-') for block in tmp[1:]])
    # now add sequences to input Alignment
    for name, seq in seqs.iteritems():
        addFunc(name, ''.join(seq), filetype=FILE_GENEDOC)
    f.close()
