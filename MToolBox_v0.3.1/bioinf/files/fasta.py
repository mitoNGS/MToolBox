# coding=utf8
"""
Formato fasta, supportata lettura e scrittura

@note: il resto delle feature dell'extended fasta della ncbi vengono ignorati,
       solo il primo (l'id) viene usato come nome della sequenza 
"""
import os

TEXT_WRAP_DEFAULT = 80
"""
Il numero di caratteri per riga nello scrivere la sequenza, evitare la modifica
"""
FILE_FASTA = 'fasta'
FILETYPE = 0
ERR_FILE_NOT_EXIST = "file %s doesn't exists"

def check_fasta(fname):
    """
    Rudimentale controllo, vede se il primo byte è un >
    
    @param fname:  il nome del file
    
    @raise ValueError: nel caso il file non esista
    """
    if not os.path.isfile(fname):
        raise ValueError(ERR_FILE_NOT_EXIST % fname)
    f = open(fname, 'U')
    if f.read(1) == '>':
        f.close()
        return True
    else:
        f.close()
        return False

def load_fasta(fname, addFunc):
    """
    Loads Alignment from Fasta File
    
    @type fname: string
    @param fname: nome del file
    @type addFunc: function
    @param addFunc: puntatore alla funzione di aggiunta
    
    @raise ValueError: nel caso il file non esista
    """
    if not os.path.isfile(fname):
        raise ValueError(ERR_FILE_NOT_EXIST % fname)
    cur_name = ""
    last_name = ""
    extras = {"filetype":FILE_FASTA}
    nseq = 0
    # Better use ''.join(list) than seq+=seq, it's faster
    cur_seq = []
    f = open(fname, 'U')
    # main loop to read file's sequences into an Alignment Object
    for line in f:
        if line[0] == '>':
            if cur_seq != []:
                # start of next sequence
                addFunc(cur_name,''.join(cur_seq), **extras)
                nseq += 1
                cur_seq = []
                gi_index = ""
            # save previous name for loop's else clause
            last_name = cur_name
            # modified to handle ncbi names (features separated by |)
            # - split by | (if not found returns a list with on element
            # - takes 1st element (always present)
            # - start with 2nd char (exclude '>')
            # - strip whitespace (\n, \t, etc)
            cur_name = line.split('|')
            if len(cur_name) > 1:
                if 'gi_index' in extras:
                    del extras['gi_index']
                #caso ncbi/ebi fasta
                if cur_name[0][1:] == 'gi':
                    #per aggiungere l'indice ncbi
                    extras['gi_index'] = cur_name[1]
                    #ncbi extended fasta, prendere solo il quarto elemento per il nome
                    cur_name = cur_name[3].split('.')[0].upper()
                elif cur_name[0][1:] == 'embl':
                    cur_name = cur_name[1].upper()
                else:
                    cur_name = cur_name[0][1:]
            else:
                #fasta classico, prende tutti i caratteri tranne il primo ">"
                #e il whitespace a destra
                cur_name = cur_name[0][1:].rstrip()
        else:
            cur_seq.append(line.rstrip())
    else:
        # handles cases where there's no newline after last
        # portion of the sequence at EOF
        if nseq != 0:
            if cur_name != last_name:
                addFunc(cur_name,''.join(cur_seq), **extras)
        # case in which only one sequence is present
        else:
            addFunc(cur_name,''.join(cur_seq), **extras)
    f.close()

def write_fasta(fname, iterFunc, twidth=TEXT_WRAP_DEFAULT):
    """
    Write a fasta file from an Alignment object
    
    @type fname: string
    @param fname: nome del file
    @type iterFunc: function
    @param iterFunc: funzione che restituisce una tupla (nome, seq)
    @type twidth: int
    @param twidth: numero massimo di caratteri per riga
    """
    f = open(fname, 'w')
    for name, seq in iterFunc():
        seqLen = len(seq)
        # sequence name
        f.write(">" + name + "\n")
        cur_pos = 0
        # some kind of textwrap, default 80 chars
        while True:
            if cur_pos + twidth <= seqLen:
                f.write(seq[cur_pos:cur_pos + twidth] + '\n')
                cur_pos += twidth
            else:
                f.write(seq[cur_pos:] + '\n')
                break
    f.close()
