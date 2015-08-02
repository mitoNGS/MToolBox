# coding=utf8
"""
Formato phylip, supportata lettura e scrittura

non ci sono controlli e si basa sul formato dato da CLC sequence viewer.
i file scritti non presenta nell'header il numero di sequenze contenute
nel file, per cui potrebbe presentare problemi nel caricamento in un altro
programma
"""
import os, re

FILE_PHYLIP = 'phylip'
FILETYPE = 0
ERR_FILE_NOT_EXIST = "file %s doesn't exists"
RE_HEADER=re.compile("(\d+) +(\d+)")

def check_phylip(fname):
    """
    Rudimentale controllo, vede se nella prima riga ci sono il numero di
    sequenze e la loro lunghezza
    
    @param fname:  il nome del file
    
    @raise ValueError: nel caso il file non esista
    """
    if not os.path.isfile(fname):
        raise ValueError(ERR_FILE_NOT_EXIST % fname)
    f = open(fname, 'U')
    header = f.readline()
    if RE_HEADER.match(header):
        return True
    else:
        return False
    

def load_phylip(fname, addFunc):
    """
    Loads Alignment from Phylip File

    Supporta il formato NOME[10]\wSEQUENZA[1-]
    il whitespace viene tolto.
    
    @type fname: string
    @param fname: nome del file
    @type addFunc: function
    @param addFunc: puntatore alla funzione di aggiunta
    
    @raise ValueError: nel caso il file non esista
    """
    if not os.path.isfile(fname):
        raise ValueError(ERR_FILE_NOT_EXIST % fname)
    extras = {"filetype":FILE_PHYLIP}
    nseq = 0 # così da controllare che il numero corrisponda
    f = open(fname, 'U')
    header = f.readline()
    tmp = RE_HEADER.match(header)
    #inutilizzati
    tot_seq = tmp.group(1)
    tot_len = tmp.group(2)
    # main loop to read file's sequences into an Alignment Object
    for line in f:
        name = line[:10].rstrip(' ')
        seq = line[10:].strip().upper()
        addFunc(name, seq, **extras)
        #inutilizzato
        nseq += 1
    f.close()

def write_phylip(fname, iterFunc):
    """
    Write a phylip file from an Alignment object

    date la limitazione di non sapere il numero di sequenze, il
    file sarà fuori standard
    
    @type fname: string
    @param fname: nome del file
    @type iterFunc: function
    @param iterFunc: funzione che restituisce una tupla (nome, seq)
    @type twidth: int
    @param twidth: numero massimo di caratteri per riga
    """
    f = open(fname, 'w')
    seqs = []
    for name, seq in iterFunc():
        seqs.append( (name, seq) )
    f.write("%d  %d\n" % (len(seqs), len(seqs[0][1])))
    for name, seq in seqs:
        f.write("% 10s  %s\n" % (name, seq))
    f.close()
