# encoding=utf8

import datatypes, itertools, pickle
from copy import copy

def traverse_haplogroup(haplogroup, aplo_list=None):
    """
    Resituisce l'elenco completo degli aplogruppi di cui fa parte quello
    passato, compreso esso stesso
    """
    if not aplo_list:
        aplo_list.append(haplogroup)
    if not haplogroup.parent is None:
        aplo_list.append(haplogroup.parent)
        traverse_haplogroup(haplogroup.parent, aplo_list)

def check_branch_slice(pos, pos_list, deleted):
    
    for x in pos_list:
        if pos.start == x.start:
            print pos.start, x.start, pos.is_reverted(x), pos.change, x.change

def group_equal(pos_list):
    """"Elimina le posizioni uguali"""
    found = set()
    if pos_list is None:
        return []
    for idx, pos in enumerate(pos_list[:-1]):
        for dup in itertools.ifilter(pos.__eq__, pos_list[idx+1:]):
            #print "Uguale:", pos, "->", dup
            found.add(hash(dup))
    return [x for x in pos_list if not hash(x) in found]

def group_reverted(pos_list):
    """Se due SNP sono retrmoutati li raggruppa in un altro tipo di classe"""
    retromutated = []
    deleted = []
    
    for idx in reversed(xrange(len(pos_list)-1)):
        retrom = [x for x in itertools.ifilter(pos_list[idx].is_reverted, pos_list[idx+1:])]
        if retrom:
            if pos_list[idx] not in deleted and retrom[0] not in deleted:
                #print "Retromutata:", pos_list[idx], retrom[0]
                retromutated.append(datatypes.Retromutated(pos_list[idx], retrom[0]))
                deleted.append(pos_list[idx]); deleted.append(retrom[0])
    #print "cancellati:", deleted
    #print "retromutati:", retromutated
    #nlist = []
    #for x in pos_list:
    #    if x in deleted:
    #        print x, "deleted"
    #    else:
    #        print x, "added"
    #        nlist.append(x)
    #return nlist, retromutated
#    return [x for x in pos_list if x not in deleted] + retromutated        
    return [x for x in pos_list if x not in deleted]

def group_equal_alt_(pos_list):
    """Quelli duplicati li trasforma in Retromutazioni
    FUNZIONA
    """
    found = set()
    
    for idx, pos in enumerate(pos_list[:-1]):
        for dup in itertools.ifilter(pos.__eq__, pos_list[idx+1:]):
            #print "Uguale, conversione in retromutazione:", pos, "->", dup
            found.add(hash(dup))
    #controllare quali sono le posizioni uguali e vedere se vale le seguenti righe (solo SNPs?)
    new_list = []
    for mut in pos_list:
        if hash(mut) in found:
            new_list.append(datatypes.Retromutation("%d!" % mut.start))
            #print "Retromutata", mut, "->", new_list[-1]
        else:
            #print "Mantenuta", mut
            new_list.append(mut)
    return new_list

def filter_positions(haplogroup, reverse=True):
    """
    
    """
    
    #genera la lista degli aplogruppi del ramo
    aplo_list = []
    traverse_haplogroup(haplogroup, aplo_list)
    if reverse:
        aplo_list.reverse()
    
    #genera la lista delle posizioni di questo ramo
    pos_list = []
    for haplo in aplo_list:
        pos_list.extend(haplo.pos_list)

    dict_set_start_list = dict(pos_list)
    new_list, weird_guys = choose_terminal_mutation(dict_set_start_list)
    return new_list
    """
    #mantiene tra una lista di posizioni uguali solo la prima
    #prima si raggruppano le posizioni uguali, facendo diventare la seconda una retromutazione
    #poi si passa ad eliminarle
    pos_list = group_equal(group_equal_alt_(pos_list))
    
    #per ultimo si raggruppano le posizioni retromutate
    tmp = group_reverted(pos_list)
    #print group_reverted(pos_list)
    #----------acrocchio
    #al momento si vuole togliere le posizioni uguali dal conto
    tmp = [x for x in pos_list if not isinstance(x, (datatypes.Retromutation, datatypes.Retromutated))]
    #-----------
    #if haplogroup.name in ['F1a1','M40']:
    #    print len(tmp), tmp
    return tmp
    """

# ==============================================================
# nuove funzioni
#
# dizionario {start : [eventi associati]} (per ogni aplogruppo)

def dict(pos_list):
    rev_pos_list = copy(pos_list)
    rev_pos_list.reverse()
    
    set_start_list = set([i.start for i in rev_pos_list])
    dict_set_start_list = {}
    
    for i in set_start_list:
        dict_set_start_list[i] = []
    
    for pos in rev_pos_list:
        if pos.start in dict_set_start_list.keys():
            dict_set_start_list[pos.start].append(pos)
    
    return dict_set_start_list
            

# questa funzione dovrebbe scegliere l'ultimo evento relativo ad un sito

def choose_terminal_mutation(dict_set_start_list):

    new_list = []
    weird_guys = {}
    
    # Analizza l'insieme di eventi associati ad una posizione.
    # Se l'ultimo evento è Retromutation, il sito viene scartato.
    
    for pos_num in dict_set_start_list.keys():
        pos_event_dict = dict_set_start_list[pos_num]
        # se l'ultimo evento (che è il primo, nella lista revertita) è una
        # retromutazione, viene scartato
        if [event.mutation_type() for event in pos_event_dict][0] == 'Retromutation':
            pass
        # se non c'è retromutazione ma comunque più di due eventi,
        # escludendo delezioni, possono essere Transversion e Transition
        elif len(pos_event_dict) > 1:
            # stampa la voce del dizionario		
            weird_guys[pos_num] = pos_event_dict
            # il prodotto di una transizione e di una trasversione è una trasversione
            # il cui change è quello dell'ultima trovata lungo l'albero
            if sorted([event.mutation_type() for event in pos_event_dict]) == ['Transition', 'Transversion']:
                new_list.append(datatypes.Transversion("%s" % str(pos_num)+pos_event_dict[0].change))
            # il prodotto di due trasversioni è una transizione
            # il cui change è quello dell'ultima trovata lungo l'albero
            elif sorted([event.mutation_type() for event in pos_event_dict]) == ['Transversion', 'Transversion']:
                new_list.append(datatypes.Transition("%d" % pos_num))
        else:
            new_list.append(dict_set_start_list[pos_num][0])
        """
		A -> G -> T
        Transition -> Transversion = Transversion

        A -> T -> C
        Transversion -> Transition = Transversion

        A -> T -> G
        Transversion -> Transversion = Transition

        A -> T -> A
        Transition -> Transition = Retromutation
		"""
    return new_list, weird_guys

# =================================================================

class HaplogroupTree(object):
    _aplo_dict = {}
    _filt_dict = {}
    _weird_dict = {}
    def __init__(self, aplo_list=None, pickle_data=None):
        if aplo_list:
            for haplogroup in aplo_list:
                self.add_haplogroup(haplogroup)
        if pickle_data:
            self.deserialize(pickle_data)
    def add_haplogroup(self, haplogroup):
        if not isinstance(haplogroup, (datatypes.Haplogroup, datatypes.MetaGroup)):
            raise TypeError("Non è stato passato un Aplogruppo")
        self._aplo_dict[haplogroup.name] = haplogroup
    def __getitem__(self, key):
        if not isinstance(key, str):
            raise TypeError
        return self._aplo_dict[key]
    def get_haplogroup_branch(self, haplo_name):
        try:
            aplo_list = []
            traverse_haplogroup(self[haplo_name], aplo_list)
            return aplo_list
        except KeyError:
            raise ValueError("Aplogruppo non trovato")
    def get_branch_positions(self, haplo_name):
        #genera la lista degli aplogruppi del ramo
        aplo_list = self.get_haplogroup_branch(haplo_name)
        aplo_list.reverse()
        
        #genera la lista delle posizioni di questo ramo
        pos_list = []
        for haplo in aplo_list:
            pos_list.extend(haplo.pos_list)
        return pos_list
    def get_filtered_positions(self, haplo_name):
        #print "Filtrazione ramo:", haplo_name
        pos_list = self.get_branch_positions(haplo_name)
        
        ## PROVA
        dict_set_start_list = dict(pos_list)
        new_list, weird_guys = choose_terminal_mutation(dict_set_start_list)
#        if weird_guys != {}:
#        for guy in weird_guys.keys():
#            print haplo_name, guy, ','.join(weird_guys[guy])
        ##
        """
        #mantiene tra una lista di posizioni uguali solo la prima
        #prima si raggruppano le posizioni uguali, facendo diventare la seconda una retromutazione
        #poi si passa ad eliminarle
        new_list = group_equal(group_equal_alt_(pos_list))
        
        #per ultimo si raggruppano le posizioni retromutate
        new_list = group_reverted(new_list)
        """
        if len(pos_list) > len(new_list):
            self._filt_dict[haplo_name] = (len(pos_list), len(new_list), new_list)
        #if weird_guys != {}:
        self._weird_dict[haplo_name] = weird_guys
        return new_list 
    def serialize(self):
        return pickle.dumps(self._aplo_dict, 0)
    def deserialize(self, data):
        self._aplo_dict = pickle.loads(data)
    def __len__(self):
        return len(self._aplo_dict)
    def __iter__(self):
        for haplo_name in self._aplo_dict:
            yield haplo_name

