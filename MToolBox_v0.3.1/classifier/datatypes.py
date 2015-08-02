# encoding=utf8

import consts

class BaseSNP(object):
    start = 0
    def __init__(self, start=None):
        if start:
            self.start = abs(int(start))
    def __eq__(self, other):
        return self.start == other.start
    def __str__(self):
        return "%s(%d)" % (self.mutation_type(), self.start)
    def __repr__(self):
        return str(self)
    #Per valutare se la mutazione è revertita
    def is_reverted(self, other):
        """Deve essere innanzitutto uguale la posizione della mutazione"""
        return self.start == other.start
    def mutation_type(self):
        return self.__class__.__name__
    def print_snp(self):
        return ''
    def print_table(self):
        return (self.start,)

class Insertion(BaseSNP):
    seq = ''
    def __init__(self, element=None):
        if element:
            pos, seq = element.split('.')
            self.seq = seq.upper()
            super(Insertion, self).__init__(pos)
    def __eq__(self, other):
        if not isinstance(other, Insertion): return False
        if super(Insertion, self).__eq__(other):
            return self.seq == other.seq
        else:
            return False
    def is_reverted(self, other):
        """
        Una inserzione revertita equivale ad una Delezione che ha come inizio
        la posizione dell'inserzione e come fine della delezione la posizione
        d'inserzione a cui va sommata la lunghezza della sequenza inserita
        """
        if not isinstance(other, Deletion): return False
        if super(Insertion, self).is_reverted(other):
            if (self.start + len(self.seq)) == other.end:
                return True
            else:
                return False
        else:
            return False
    #def pprint(self):
    #    return str(self) + " -> " + self.seq
    def pprint(self):
        return str(self.start) + "." + self.seq
    def print_snp(self):
        return super(Insertion, self).print_snp() + self.seq
    def print_table(self):
        return super(Insertion, self).print_table() + ('I', '', self.seq, 0, '')
    def __hash__(self):
        return hash((self.mutation_type(), self.start, self.seq))


class Deletion(BaseSNP):
    end = 0
    def __init__(self, element=None):
        if element:
            if element.find('-') > -1:
                start, end = element.split('-')
            else:
                start = element[:-1]
                end = element
            self.end = abs(int(end[:-1]))
            super(Deletion, self).__init__(start)
    def __eq__(self, other):
        if not isinstance(other, Deletion): return False
        if super(Deletion, self).__eq__(other):
            return self.end == other.end
        else:
            return False
    def is_reverted(self, other):
        """
        Una delezione revertita equivale ad una Inserzione che ha come
        posizione d'inserzione l'inizio della Delezione e come sequenza 
        d'inserzione quella che si trova nel range della Delezione in
        rCRS

        Nota: forse è inutile il check della lunghezza. Basterebbe fare
              un check sulla sequenza.
        """
        if not isinstance(other, Insertion): return False
        if super(Deletion, self).is_reverted(other):
            if (other.start + len(other.seq)) == self.end:
                return other.seq == consts.RCRS[self.start-1:self.end]
            else:
                return False
        else:
            return False
    #def pprint(self):
    #    return str(self) + " -> " + consts.RCRS[self.start-1:self.end]
    def pprint(self):
        if self.start == self.end:
            return str(self.start) + "d"
        else:
            return str(self.start) + "-" + str(self.end) + "d"
    def print_snp(self):
        return super(Deletion, self).print_snp() + 'end: ' + str(self.end)
    def print_table(self):
        return super(Deletion, self).print_table() + ('D', '', '', self.end, '')
    def __hash__(self):
        return hash((self.mutation_type(), self.start, self.end))

class SNP_MixIn(object):
    """
    Contiene metodi condivisi dalle classi che dichiarano uno SNP
    
    Nota: porre come prima classe nella dichiarazione
    """
    def __eq__(self, other):
        if self.__class__ == other.__class__:
            if super(SNP_MixIn, self).__eq__(other):
                try:
                    return self.change == other.change
                except AttributeError:
                    pass
        return False
	# NEWLY ADDED on Oct18 2013
    def pprint(self):
		if hasattr(self, 'ambiguity'):
			return str(self.start) + self.change + '(' + self.ambiguity + ')'
		else:
			return str(self.start) + self.change
    #def pprint(self):
    #    return str(self) + " " + consts.RCRS[self.start-1] + " -> " + self.change
    def is_reverted(self, other):
        """
        La direzione delle due mutazioni è importante. la prima (quella
        che richiama questo metodo è quella ancestrale (self), quella
        invece che viene passata (other) è quella più a valle nel ramo.
        
        Premesso che la posizione della mutazione deve essere la stessa,
        il cambio di quella a valle deve essere necessariamente uguale ad
        anderson. Entrambe devono essere una tra i seguenti tipi di mutazione:
        Transizione, Trasversione o Retromutazione.
        """
        if super(SNP_MixIn, self).is_reverted(other):
            try:
                #è revertita solo sel il cambio nell'altra mutazione è uguale
                #ad anderson
                return consts.RCRS[self.start-1] == other.change
            except AttributeError:
                pass
        return False
    def print_snp(self):
        return super(SNP_MixIn, self).print_snp() + consts.RCRS[self.start-1] + ' -> ' + self.change
    def print_table(self):
        return super(SNP_MixIn, self).print_table() + ('S', self.change, '', 0, '')
    def __hash__(self):
        return hash((self.mutation_type(), self.start, self.change))

class Transition(SNP_MixIn, BaseSNP):
    def __init__(self, element=None, start=0, change=''):
        try:
            if element:
                super(Transition, self).__init__(element)
                self.change = consts.TRS_TBL[consts.RCRS[self.start-1]]
            else:
                #nel caso si voglia modificare dopo i cambi
                super(Transition, self).__init__(start)
                self.change = change
        except KeyError:
            print self.start-1
            pass
    #il valore revertito dell'istanza Trs(100, 'A')->Trs(100,'G')
    def revert(self):
        return Transition(start=self.start, change=consts.TRS_TBL[self.change])

class Transversion(SNP_MixIn, BaseSNP):
    def __init__(self, element=None):
        try:
            if element:
                #changed since it returned an error. changed according to similar definition in Transition class.
                #super(Transversion, self).__init__(element[:-1])
                super(Transversion, self).__init__(element[:-1])
                self.change = element[-1].upper()
            else:
                super(Transversion, self).__init__(0)
                self.change = ''
        except KeyError:
            print self.start-1
            pass

class Unknown(SNP_MixIn, BaseSNP):
    pass

class Retromutation(Transition):
    def __init__(self, element=None):
        if element:
            super(Retromutation, self).__init__(element[:-1])
            self.change = consts.RCRS[self.start-1]
        else:
            super(Retromutation, self).__init__()

class Haplogroup(object):
    def __init__(self, name, parent=None, pos_list=None):
        self.parent = parent
        self.pos_list = [x for x in pos_list]
        self.name = name
    def __contains__(self, item):
        return item in self.pos_list
    def __eq__(self, other):
        if not issubclass(other.__class__, Haplogroup): return False
        if self.parent == other.parent and self.pos_list == other.pos_list: return True
        return False
    def __ne__(self, other):
        return not self.__eq__(other)
    def __str__(self):
        return "%s(%d)" % (self.name, len(self.pos_list))
    def __repr__(self):
        return str(self)
    def __iter__(self):
        for x in self.pos_list:
            yield x

class MetaGroup(Haplogroup):
    def __init__(self, name, parent=None, pos_list=None):
        super(MetaGroup, self).__init__(name, parent, pos_list)
        if "'" in name:
            self.groups = self._apix_sep(name)
        elif '-' in name:
            self.groups = self._line_sep(name)
        else:
            raise TypeError('Tipo di metagruppo non riconosciuto')
        self.pos_list = pos_list
    def _sep(self, groups, sep="'"):
        suffix = groups.split(sep)
        prefix = suffix[0][:-1]
        suffix[0] = suffix[0][-1]
        return prefix, suffix
    def _apix_sep(self, groups):
        prefix, suffix = self._sep(groups, "'")
        return tuple(prefix+x for x in suffix)
    def _line_sep(self, groups):
        prefix, suffix = self._sep(groups, "-")
        start = consts.CHARS.index(suffix[0])
        end = consts.CHARS.index(suffix[1]) + 1
        return tuple(prefix+consts.CHARS[x] for x in range(start, end))
    def __contains__(self, item):
        #se ci si chiede se un aplogruppo fa parte di questo metagruppo
        if isinstance(item, (str, unicode)):
            return item in self.groups
        else:
            #gli altri casi al momento cotemplano esclusivamente gli interi
            return super(MetaGroup, self).__contains__(item)
    def __eq__(self, other):
        if not issubclass(other.__class__, MetaGroup): return False
        if not super(MetaGroup, self).__eq__(other): return False
        if self.groups == other.groups: return True
        return False
    def __ne__(self, other):
        return not self.__eq__(other)

class Retromutated(object):
    """
    Contiene una coppia di mutazioni che si possono escludere (retromutata)
    """
    mutations = []
    start = None
    def __init__(self, anc, desc):
        self.mutations = [anc, desc]
        self.start = anc.start
    def __eq__(self, other):
        """
        Funziona solo se il termine a sinistra è Retomutated
        """
        return other in self.mutations
    def __ne__(self, other):
        return not self.__eq__(other)
    def __str__(self):
        return "%s(%s)" % (self.__class__.__name__, ','.join( str(x) for x in self.mutations))
    def __repr__(self):
        return str(self)

class Sequence(object):
    name = None
    seq = None
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

def detect_feature(element):
    value = None
    #delezione
    if element[-1] == 'd':
        value = Deletion(element)
    #inserzione
    elif element.find('.') > -1:
        value = Insertion(element)
    #Trasversione
    elif element[-1].upper() in consts.DNA:
        value = Transversion(element)
    #Retromutazione
    elif element[-1] == '!':
        value = Retromutation(element) 
    else:
        #gli altri casi sono Transizioni
        value = Transition(element)
    return value

def detect_haplotype(element):
    #i "meta" aplogruppi contengono i caratteri "-" o "'"
    if "'" in element or "-" in element:
        return MetaGroup
    else:
        return Haplogroup

