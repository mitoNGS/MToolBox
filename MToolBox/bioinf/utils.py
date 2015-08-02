# coding=utf8
from random import choice, seed, randint, sample

#: tabella usata per le transizioni
trs_tbl = {
           'A': 'G',
           'G': 'A',
           'C': 'T',
           'T': 'C'
          }

#: purine e pirimidine
pur = ('A','G')
pyr = ('T','C')

__deprecated__ = "genRandSeq"

class autoprop(type):
    """
    Classe scopiazzata che crea proprietà per una classe in cui si tovi un metodo
    che inizi con _get_ o _set_, imposta anche la docstring della proprietà se
    è presente _doc_ della proprietà
    
    Es. _get_value imposta per la proprietà value il getter _get_value
    """
    def __init__(cls, name, bases, dic):
        props = {}
        for name in dic:
            if name.startswith("_get_") or name.startswith("_set_"):
                props[name[5:]] = 1
        for name in props:
            fget = getattr(cls, '_get_%s' % name, None)
            fset = getattr(cls, '_set_%s' % name, None)
            doc  = getattr(cls, '_doc_%s' % name, None)
            setattr(cls, name, property(fget=fget, fset=fset, doc=doc))

def gen_rand_seq(length, comp={'A':0.25, 'G':0.25, 'C':0.25, 'T':0.25}):
    """
    Genera una sequenza nucleotidica casuale
    
    @type length: int
    @param length: la lunghezza della sequenza da generare
    @type A: float
    @param A: percentuale di A
    @type G: float
    @param G: percentuale di G
    @type C: float
    @param C: percentuale di C
    @type T: float
    @param T: percentuale di T
    
    @rtype: string
    @return: sequenza di nucleotidi
    
    @note: I valori di A, C, T, G sono in percentuale
    @note: Non viene fatto un controllo, ma la loro somma deve dare 1
    """

    nuc = ''.join([ str(x) * int(100 * comp[x]) for x in comp])
    #nuc = ['A'] * int(100 * A) + ['C'] * int(100 * C) +\
          #['T'] * int(100 * T) + ['G'] * int(100 * G)
    #inutile?
    #seed(tuple(nuc*randint(0,randint(3000, 15000))))
    seq = []
    for x in xrange(0, length):
        seq.append(choice(nuc))
    return ''.join(seq)

genRandSeq = gen_rand_seq

def rand_seq_nuc(seq, sost_rate=0.25, trs_rate=0.75):
    """
    Resituisce una sequenza a cui sono stati cambiate delle posizioni
    Si decide una percentuale di cambiamenti (B{sost_rate}, viene creata una lista
    di posizioni da cambiare casuale e di queste si prende una lista di transizioni
    rappresentata in percentuale da B{trs_rate}; i restanti sono le trasversioni
    che sono saranno applicate
    
    @type seq: string
    @param seq: sequenza da cambiare
    @type sost_rate: float
    @param sost_rate: percentuale di sostituzioni (M{0 < n_sost S{<=} 1})
    @type trs_rate: float
    @param trs_rate: percentuale di transizioni, il resto saranno trasversioni
                     (M{0 S{<=} trs_rate S{<=} 1})
    
    @rtype: string
    @return: sequenza di nucleotidi
    
    @note: le ambiguità non vengono cambiate
    """
    
    #trasformarla in lista ci permette di sostituire velocemente i valori
    #nel caso seq sia una stringa 
    seq = list(seq)
    
    #il totale delle sostituzioni da eseguire
    sost = set(sample(xrange(len(seq)), int(len(seq)*sost_rate)))
    trs = set(sample(sost, int(len(sost)*trs_rate)))
    tsv = sost - trs
    
    for idx in trs:
        nuc = seq[idx]
        if nuc in trs_tbl:
            seq[idx] = trs_tbl[nuc]
    
    for idx in tsv:
        nuc = seq[idx]
        if nuc in pur:
            seq[idx] = choice(pyr)
        elif nuc in pyr:
            seq[idx] = choice(pur)
    
    return ''.join(seq)

def rand_seq_del(seq, del_rate=0.05, del_char='-'):
    """
    Inserisce delezioni nella sequenza, di default il numero di delezioni e' il
    5% della lunghezza della sequenza e alla posizione cancellata viene
    sostituito un gap
    
    @type seq: string
    @param seq: la sequenza da cambiare
    @type del_rate: float
    @param del_rate: il tasso di delezione (M{0 < del_rate < 1})
    @type del_char: string
    @param del_char:  il carattere con cui viene sostituito la posizione, in alternativa
                      è possibile usare None o un qualunque valore negativo
                      ("", 0, ecc.) per cancellare la posizione
    
    @rtype: string
    @return: sequenza di nucleotidi
    """
    
    del_pos = sample(xrange(len(seq)), int(len(seq)*del_rate))
    
    seq = list(seq)
    
    for idx in del_pos:
        if del_char:
            seq[idx] = del_char
        else:
            #nel caso in cui le posizioni vadano cancellate e' piu' veloce usare
            #un replace sulla stringa finale contro un carattere che non puo'
            #esserci nella stringa iniziale/finale
            seq[idx] = '#'
    
    return ''.join(seq).replace('#','')

def rand_alg(alg_len, seq_len, del_rate=None, gen_rand_seq_args=None, rand_seq_nuc_args=None):
    """
    Crea un allineamento casuale, di cui la prima è quella generata da L{gen_rand_seq},
    mentre le altre hanno due passaggi: prima rand_seq_nuc e poi L{rand_seq_del}
    (quest'ultimo SOLO se del_rate diverso da 0 o None) per evitare che possano
    esserci siti "vuoti" (solo gap)
    
    @type alg_len: int
    @param alg_len: numero di sequenze nell'allineamento
    @type seq_len: int
    @param seq_len: lunghezza delle sequenze nell'allineamento
    
    @type del_rate: float
    @param del_rate: il tasso di delezione (M{0 S{<=} del_rate < 1})
    @type gen_rand_seq_args: dict
    @param gen_rand_seq_args: dizionario contenente i valori da passare a L{gen_rand_seq} (tranne length)
    @type rand_seq_nuc_args: dict
    @param rand_seq_nuc_args: dizionario contenente i valori da passare a L{rand_seq_nuc} (tranne seq)
    
    @rtype: generator
    @return: sequenze nucleotidiche random   
    """
    
    if gen_rand_seq_args:
        seq = gen_rand_seq(seq_len, gen_rand_seq_args)
    else:
        seq = gen_rand_seq(seq_len)
    
    #restituiamo la prima come riferimento
    yield seq
    
    for x in xrange(alg_len - 1):
        if rand_seq_nuc_args:
            seq_new = rand_seq_nuc(seq, **rand_seq_nuc_args)
        else:
            seq_new = rand_seq_nuc(seq)
        if del_rate:
            seq_new = rand_seq_del(seq_new, del_rate)
        yield seq_new

def rand_list(list_len, seq_max_len, seq_min_len=None, gen_rand_seq_args=None):
    """
    Generazione di una lista di sequenze casuali
    numero di sequenze list_len, lunghezza massima seq_max_len.

    se specificato seq_min_len, verrà estratto una porzione casuale della
    sequenza di lunghezza compresa tra seq_max_len e seq_min_len
    l'altro parametro è un dizionario che passa la composizione nucletidica
    alla funzione gen_rand_seq
    """
    for x in xrange(list_len):
        
        if gen_rand_seq_args:
            seq = gen_rand_seq(seq_max_len, gen_rand_seq_args)
        else:
            seq = gen_rand_seq(seq_max_len)

        if seq_min_len:
            seq = ''.join(sample(seq, randint(seq_min_len, seq_max_len)))

        yield seq


SNP_TRS = 0
SNP_TSV = 1
SNP_DEL = 2
SNP_INS = 3
TRANSITIONS = {'A':'G','G':'A','C':'T','T':'C'}

class BadNucleotideError(Exception):
    """
    Sollevata nel caso non sia applicabile la tabella delle transizioni
    
    In pratica quando il nucleotide da controllare è un'ambiguità
    """
    pass

def snp_to_seq(ref, snps):
    """
    Partendo da una sequenza di riferimento, genera una sequenza modificata
    sulla base degli SNPs dati
    
    elementi di snps:
    (posizione, tipo, dati)
    
        - posizione: numero che indica la posizione dello SNP (pos S{>=} 0)
        - tipo: codice indicante il tipo:
            - 0: transizione
            - 1: trasversione
            - 2: delezione
            - 3: inserzione
        - dati varia per tipo:
            - transizione: automatizzato, se presente utilizzato (in questo caso
                           equivale ad una trasversione)
            - trasversione: cambio nucleotidico (C se A S{->} C, ecc.)
            - delezione: prima posizione non interessata dalla delezione
            - inserzione: sequenza inserita prima della posizione 
    
    @raise BadNucleotideError: nel caso non sia possibile discernere la transizione
           o non sia specificato un nucleotide per la trasversione
    @note: nell'indicare le posizioni, usare sempre le regole degli slice
    
    @type ref: string
    @param ref: sequenza di riferimento
    @type snps: list
    @param snps: lista di SNPs (deve avere il metodo sort con reverse come argomento) 
    @rtype: string
    @return: sequenza modificata
    """
    ref = list(ref.upper())
    
    snps.sort(reverse=True)
    
    for pos, ctype, data in snps:
        if ctype == SNP_TRS:
            if data:
                #nel caso sia stato specificato un cambio
                ref[pos] = data
            else:
                #cambio automatico
                try:
                    ref[pos] = TRANSITIONS[ref[pos]]
                except KeyError:
                    raise BadNucleotideError("unknown: %s" % ref[pos])
        elif ctype == SNP_TSV:
            if not data:
                raise BadNucleotideError("no nucleotide specified")
            ref[pos] = data
        elif ctype == SNP_DEL:
            del ref[pos:data]
        elif ctype == SNP_INS:
            ref.insert(pos, data.upper())
    return ''.join(ref)
        
