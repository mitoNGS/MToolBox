#!/usr/bin/env python
# encoding=utf8

#tenere presente:
#- i numeri alla cui sinistra non c'è un aplogruppo, vanno considerati come
#  comuni a tutti gli aplogruppi al di sotto (purché abbiano un nome valido di aplogruppo/metaaplogruppo)
#- una posizione è revertita se è si ritrova a monte nell'albero (segnare come permanenza?)
#- la notazione N1'5 indica che le posizioni elencate sono da considerarsi per gli aplogruppi N1 ed N5
#- L1-6 indica che le posizioni elencate vanno considerate per L1, L2, ..., L6 
#- un ! alla destra di una posizione indica una retromutazione, per cui il nucleotide sarà uguale ad rCRS

# import sys
# if '..' not in sys.path:
#     sys.path.append('..')

from classifier import datatypes

def common_filter(element):
    "Resituisce None se l'elemento è da scartare"
    #i numeri fa parentesi non vanno usati
    #if element.find('(') > -1 or element.find(')') > -1: return ''
    #li usiamo ma togliamo le parentesi, e i frequenti doppi spazi
    element = element.replace('(','').replace(')','').replace('  ','')
    #sembra che questi caratteri siano presenti nella conversione
    #excel->csv, perlomeno con OO o Gnumeric
    return element.replace('\xa0', '').replace('\xc2', '').strip()

def detect_haploname(element):
    "Verifica che l'elemento sia il nome di un aplogruppo"
    return element[0].isalpha()

def find_level(line):
    "Restituisce il livello della linea (conta le caselle vuote a sinistra)"
    count = 0
    for x in line:
        if not x.strip():
            count += 1
        else:
            break
    return count

def parse_line(line, level, parent=None):
    name = common_filter(line[level])
    #controllo il tipo di class di aplogruppo da usare
    haplo_class = datatypes.detect_haplotype(name)
    positions = []
    for pos in line[level+1:]:
        filtered = common_filter(pos)
        #se l'elemento non è vuoto oppure da scarta passa all'iterazione succesiva
        if filtered:
            positions.append(datatypes.detect_feature(filtered))
    return haplo_class(name, parent=parent, pos_list=positions)

def parse_csv(file_handle, aplo_list=None, parent=None):
    if aplo_list is None: aplo_list = []
    
    parent_level = [parent]
    
    for n, line in enumerate(file_handle):
        try:
			# trova il "livello" della riga (conta le caselle vuote a sinistra)" 
            level = find_level(line)
        except IndexError:
            print "WARNING: La riga:", n+1, "è vuota"
        #print n, line[:3], level
        if len(parent_level) - 1 < level:
            parent_level.append(aplo_list[-1])
        elif len(parent_level) -1 > level:
            parent_level = parent_level[:level+1]
        try:
            value = parse_line(line, level, parent_level[-1])
        except ValueError, e:
            print "parent level after, ", parent_level
            print "level is ", level
            print value
            print "ERROR: Più di un valore per campo alla riga:", n+1, "oppure un aprogruppo contiene sia * che ' o -", e.args, line
            #sys.exit(2)
            raise ValueError
        except IndexError, e:
            print "WARNING: La riga:", n+1, "è vuota oppure contiene una posizione maggiore della lunghezza di rCRS - ", "(", str(e), ") - ", line
        aplo_list.append(value)  
    return aplo_list
