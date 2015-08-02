# coding=utf8
"""
Done:
    - Creato oggetto composition e modificato il comportamento di chi ne fa uso:
    RegEx e SiteVar (lo script)
ToDo:
    - Eliminare il calcolo della composizione dalla classe Alignment
    - Cambiare i nomi delle funzioni per mantenere una notazione some_funtion invece
      di someFunction, i nomi delle classi vanno bene cosi'
    - Col tempo togliere i metodi/attributi nascosti __var e sostituirli coi meno
      problematici _var (che si sottintende non vadano toccati)
"""

__deprecated__ = ""

#python 2.6 raises a DeprecationWarning if sha is not imported from
#hashlib module
try:
    from hashlib.sha1 import sha
except ImportError:
    import sha

from array import array
import os, imp
from files import write_funcs, load_funcs, check_funcs, file_types, plugins
from utils import autoprop
import comp

ERR_SEQ_LEN   = 'Sequence in not of correct length, must be %d'

#this are some kind of gaps found in certain sequences:
#? was used for an old version of sitev and maybe other software as well
#. is used in genedoc files
gap_types = ['?','.']

class BadSequenceType(Exception):
    """Sollevata quando i tipi di sequenze non coincidono"""
    pass

def contained(seq1, seq2):
    """Check if seq1 contains seq2"""
    if seq1.find(seq2) != -1:
        return True
    else:
        return False

def check_item(item):
    """
    These checks and conversions are used for some
    comparison functions, less than, more than, etc.
    """
    if not isinstance(item, (str, array, BaseSequence)):
        raise TypeError
    if isinstance(item, str):
        item = item.upper()
    if isinstance(item, BaseSequence):
        item = item.seq
    if isinstance(item, array):
        item = item.tostring().upper()
    return item

def calc_hash(seq):
    """
    Returns an SHA1 hash of seq
    
    It's used in comparison functions and when seq is assigned
    to a Sequence Class
    
    @type seq: sequence
    @param seq: a string of chars
    """
    return array('c', sha.new(seq).digest())

class BaseSequence(object):
    """
    Base Class for all Sequences types
    
    The only (by now) method overridden by Classes
    derived by BaseSequence is __init__ cooperatively
    called by all BaseSequence derived classes
    
    @note: You shouldn't use this class but L{Sequence} or one of
           its derived classes
    """
    _contained = staticmethod(contained)
    _check_item = staticmethod(check_item)
    _calc_hash = staticmethod(calc_hash)
    __metaclass__ = autoprop
    def __init__(self, name='', seq=''):
        self._set_seq(seq)
        self._set_name(name)
    def _set_seq(self, seq):
        if not isinstance(seq, str):
            raise TypeError
        for gap in gap_types:
            seq = seq.replace(gap,'-')
        self._seq = array('c', seq.upper())
        self._hash = self._calc_hash(self._seq)
    def _get_seq(self):
        return self._seq
    _doc_seq = 'Set/Get sequence'
    def _set_name(self, name):
        if not isinstance(name, str):
            raise TypeError
        self._name = name
    def _get_name(self):
        return self._name
    _doc_name = 'Set/Get name of sequence'
    def _get_hash(self):
        return self._hash
    _doc_hash='Get SHA1 Hash for sequence'
    def __len__(self):
        """Returns length of sequence"""
        return len(self._seq)
    def __contains__(self, item):
        """
        Checks if item is a substring of sequence
        
        item can be a string, an array of caracters, of another
        istance of BaseSequence and its derivative
        """
        item = self._check_item(item)
        return self._contained(self._seq.tostring(), item)
    def __repr__(self):
        """Used by print statement"""
        return "%s: %s" % (self._name, self._seq.tostring())
    def __str__(self):
        """Convert to string"""
        return self._seq.tostring()
    def __lt__(self, item):
        item = self._check_item(item)
        if len(self._seq) < len(item):
            return self._contained(item, self._seq.tostring())
        else:
            return False
    def __le__(self, item):
        item = self._check_item(item)
        if len(self._seq) <= len(item):
            return self._contained(item, self._seq.tostring())
        else:
            return False
    def __gt__(self, item):
        item = self._check_item(item)
        if len(self._seq) > len(item):
            return self._contained(self._seq.tostring(), item)
        else:
            return False
    def __ge__(self, item):
        item = self._check_item(item)
        if len(self._seq) >= len(item):
            return self._contained(self._seq.tostring(), item)
        else:
            return False
    def __eq__(self, item):
        if not isinstance(item, (str, array, BaseSequence)):
            raise TypeError
        if isinstance(item, array):
            item = item.tostring()
        if isinstance(item, str):
            hash = self._calc_hash(item.upper())
        if isinstance(item, BaseSequence):
            hash = item.hash
        if self._hash == hash:
            return True
        else:
            return False
    def __ne__(self, item):
        if not isinstance(item, (str, array, BaseSequence)):
            raise TypeError
        if isinstance(item, array):
            item = item.tostring()
        if isinstance(item, str):
            hash = self._calc_hash(item.upper())
        if isinstance(item, BaseSequence):
            hash = item.hash
        if self._hash != hash:
            return True
        else:
            return False
    def __nonzero__(self):
        return bool(self._seq)
    # indexing and slicing
    def __getitem__(self, item):
        if isinstance(item, (int, long)):
            return self._seq[item]
        elif isinstance(item, slice):
            return self._seq[item.start:item.stop:item.step]
        else:
            raise TypeError
    def __delitem__(self, item):
        if isinstance(item, (int, long)):
            del self._seq[item]
            self._calc_hash(self._seq)
        elif isinstance(item, slice):
            del self._seq[item.start:item.stop:item.step]
            self._hash = self._calc_hash(self._seq)
        else:
            raise TypeError
    def __setitem__(self, item, value):
        """
        Value can be an array or a string of another instance
        
        @attention: la lunghezza di value deve essere uguale a quella della stringa
                    sostituita
        """
        #uppercase if string
        if isinstance(value, str):
                value = value.upper()
        
        if isinstance(item, (int, long)) and len(value) == 1:
            if not isinstance(value, str):
                raise TypeError
            self._seq[item] = value
            self._calc_hash(self._seq)
        elif isinstance(item, slice):
            if isinstance(value, str):
                value = array('c', value)
            elif isinstance(value, array):
                value = array('c', value.tostring().upper())
            elif isinstance(value, BaseSequence):
                value = value.seq
            self._seq[item.start:item.stop:item.step] = value
            self._hash = self._calc_hash(self._seq)
        else:
            raise TypeError
    # iterators
    def __iter__(self):
        for item in self._seq:
            yield item
    # emulates numeric types
    def __add__(self, item):
        """Accept another instance or a string"""
        if isinstance(item, BaseSequence):
            self._seq = self._seq + item.seq
        elif isinstance(item, str):
            self._seq.fromstring(item.upper())

class Sequence(BaseSequence):
    """
    This class can associates a dictionary indicating its composition
    
    It serves the purpose of identifying an erroneus sequence
    via checkSeq method, you have to choose one from L{comp.comp_types}
    """
    def __init__(self, name='', seq='', comp_type=None):
        try:
            self._comp_type = comp.comp_types[comp_type]
        except KeyError:
            self._comp_type = {}
        super(Sequence, self).__init__(name, seq)
    def calc_comp(self):
        comp = self._comp_type.copy()
        seq = self.seq
        for key in comp:
            comp[key] = seq.count(key)
        self._comp = comp
    def _get_comp(self):
        if not hasattr(self, '_comp'):
            self.calc_comp()
        return self._comp
    _doc_comp = 'Get composition of sequence'
    def check_seq(self):
        for x in self.seq:
            if x not in self._comp_type:
                return False
        return True

class NucSequence(Sequence):
    """
    Gives quick access to a nucleotide sequence without
    having to specify a compType in Sequence __init__
    moreover NucSequence has other methods for translating
    its sequence in an amminoacid sequence etc.
    """
    def __init__(self, name='', seq=''):
        super(NucSequence, self).__init__(name, seq, 'dna')
    def translate(self, tbl=None):
        """To be implemented"""
        pass
class AaSequence(Sequence):
    def __init__(self, name='', seq=''):
        super(NucSequence, self).__init__(name, seq, 'aa')
    

class NCBI_Entry(Sequence):
    """Non usare"""
    def __init__(self, info=None):
        if not info:
            name = seq = ''
        else:
            name = info['locus']['name']
            self.locus = info['locus'].copy()
            self.definition = info['definition']
            self.version = info['version'].copy()
            self.keywords = info['keywords']
            self.source = []
            self._add_values(self.source, info, 'source')
            self.reference = []
            self._add_values(self.reference, info, 'reference')
            self.features = []
            self._add_values(self.features, info, 'features')
            seq = ''.join(info['origin'])
            if self.locus['seqtype'].upper() == 'DNA':
                comp = 'dna'
            else:
                comp = None
        super(NCBI_Entry, self).__init__(name, seq, comp)
    def _add_values(self, var, info, c_key):
        for item in info[c_key]:
                var.append( {} )
                c_item = var[-1]
                for key, value in item.iteritems():
                    if isinstance(value, list):
                        if key == 'pos':
                            tmp = []
                            for pos in value:
                                tmp.append( tuple([ int(x) for x in pos.split('..')]) )
                            value = tmp
                        else:
                            value = ''.join(value)
                    c_item[key] = value

class FileIO(object):
    __metaclass__ = autoprop
    def __init__(self):
        self._file_formats = plugins
        self._file_types = file_types
        self._check_funcs = check_funcs
        self._load_funcs = load_funcs
        self._write_funcs = write_funcs
        self._default_format = 'fasta'
    def _get_write_formats(self):
        return list(self._write_funcs.keys())
    #writeFormats = property(_get_writeFormats)
    def _get_default_format(self):
        return self._default_format
    def _set_default_format(self, format='fasta'):
        if format not in self._write_funcs:
            format = 'fasta'
        self._default_format = format
    #defaultFormat = property(_get_defaultFormat, _set_defaultFormat)
    def load_file(self, fname, add_func=None):
        if not add_func:
            raise ValueError('nessuna funzione di aggiunta dichiarata')
        l_func = None
        for key, func in self._check_funcs.iteritems():
            if func(fname):
                l_func = self._load_funcs[key]
        if not l_func:
            raise TypeError("non e' stato possibile determinare il formato del file")
        l_func(fname, add_func)
    def write_file(self, fname, format=None, get_func=None):
        if not format:
            format = self.default_format
        try:
            w_func = self._write_funcs[format]
        except KeyError:
            raise ValueError('formato di file errato')
        if not get_func:
            raise ValueError('nessuna funzione di aggiunta dichiarata')
        w_func(fname, get_func)

class SeqList(FileIO):
    """
    It gives access to a list of sequences
    
    FileFormats declares functions to load/write
    sequence from/to a file
    
    Notes:
        - getSeqByIndex is manteined for compatibility
          with old package, but will probably removed
    """
    __metaclass__ = autoprop
    def __init__(self, seq_cls=BaseSequence):
        super(SeqList, self).__init__()
        self._seqs = []
        self._seq_cls = seq_cls
    def add_seq(self, name='', seq='', **kwds):
        """Adds a sequence, keywords arguments are ignored for now"""
        if 'filetype' in kwds:
            if kwds['filetype'] == 'ncbi_flat':
                #caso in cui sia una flat entry NCBI
                self._seqs.append(NCBI_Entry(kwds['info']))
                return
        #se la sequenza è una stringa
        if isinstance(seq, str):
            self._seqs.append(self._seq_cls(name, seq))
        #se è un'istanza della classe usata
        elif isinstance(seq, self._seq_cls):
            self._seqs.append(seq)
        elif isinstance(seq, array):
            self.seqs.append(self._seq_cls(name, seq.tostring()))
        #l'indice che va usato per prendere dati all'ncbi
        if 'gi_index' in kwds:
            self._seqs[-1].gi_index = kwds['gi_index']
    def remove_seqs(self):
        """Removes ALL Sequences"""
        self._seqs = []
    def get_seq_by_index(self, index):
        """
        Returns Sequence by Position in the Alignment (deprecated)
        use a an index SeqList[index], iteration or slice instead
        """
        return self[index]
    def get_seq_by_name(self, name):
        """Returns Sequence by Name in the Alignment
        
        First sequence which have that name will be returned"""
        for seq in self:
            if seq.name == name:
                return seq
        return None
    def get_seq_by_hash(self, hash):
        """Returns Sequence by Hash in the Alignment
        
        First sequence which have that hash will be returned"""
        for seq in self:
            if seq.hash == hash:
                return seq
        return None
    # Getters for informations (iterators)
    def get_index_by_name(self, name):
        """Return index of sequence of given name"""
        for idx, seq in enumerate(self):
            if seq.name == name:
                return idx
    def get_index_by_hash(self, hash):
        """Return index of sequence of given name"""
        for idx, seq in enumerate(self):
            if seq.hash == hash:
                return idx
    def get_pos(self):
        """Returns (position, name)"""
        for pos, seq in enumerate(self):
            yield (pos, seq.name)
    def get_items(self):
        """Returns (name, seq), used by writeFile of FileFormats"""
        for seq_cls in self:
            yield (seq_cls.name, str(seq_cls))
    def get_duplicates(self):
        dupl = {}
        for idx, seq in enumerate(self):
            try:
                dupl[seq.hash.tostring()].append(idx)
            except KeyError:
                dupl[seq.hash.tostring()] = [idx]
        for dup in dupl.itervalues():
            yield dup
    def get_hash(self):
        "ritorna un hash sha1 calcolato usando gli hash di tutte le sequenze"
        hsh = array('c')
        for seq in self:
            hsh += seq.hash
        return calc_hash(hsh)
    def load_file(self, fname):
        super(SeqList, self).load_file(fname, self.add_seq)
    def write_file(self, fname, format=None):
        super(SeqList, self).write_file(fname, format, self.get_items)
    # --------------------------------------------
    # Special Methods
    # --------------------------------------------
    # used by len() builtin
    def __len__(self):
        """Length of internal list"""
        return len(self._seqs)
    # used by "in" keyword
    def __contains__(self, item):
        if item in self._seqs:
            return True
        else:
            return False
    def __repr__(self):
        data = []
        data.append( "Contains %d sequence(s)" % len(self) )
        for seq in self:
            data.append( '%-12s: %s' % (seq.name[:12], seq.seq.tostring()) )
        return '\n'.join(data)
    # rich comparison methods
    def __lt__(self, other):
        if len(other) > len(self):
            return True
        else:
            return False
    def __le__(self, other):
        if len(other) >= len(self):
            return True
        else:
            return False
    def __eq__(self, other):
        if not isinstance(other, SeqList):
            raise TypeError
        for seq1, seq2 in zip(self, other):
            if seq1 != seq2:
                return False
        return True
    def __ne__(self, other):
        if not isinstance(other, SeqList):
            raise TypeError
        for seq1, seq2 in zip(self, other):
            if seq1 != seq2:
                return True
        return False
    def __gt__(self, other):
        if len(other) < len(self):
            return True
        else:
            return False
    def __ge__(self, other):
        if len(other) <= len(self):
            return True
        else:
            return False
    def __nonzero__(self):
        return bool(self._seqs)
    # indexing and slicing
    def __getitem__(self, item):
        """
        E' possibile richiedere la/e sequenza/e tramite indice, slice, nome o hash
        """
        if isinstance(item, (int, long)):
            return self._seqs[item]
        elif isinstance(item, slice):
            return self._seqs[item.start:item.stop:item.step]
        elif isinstance(item, str):
            return self.get_seq_by_name(item)
        elif isinstance(item, array):
            return self.get_seq_by_hash(item) 
        else:
            raise TypeError
    def __delitem__(self, item):
        if isinstance(item, (int, long)):
            del self._seqs[item]
        elif isinstance(item, slice):
            del self._seqs[item.start:item.stop:item.step]
        else:
            raise TypeError
    def __setitem__(self, item, value):
        if isinstance(item, (int, long)):
            if value.__class__.__name__ != self._seq_cls.__name__:
                raise BadSequenceType
            self._seqs[item] = value
        elif isinstance(item, slice):
            for x in value:
                if x.__class__.__name__ != self._seq_cls.__name__:
                    raise BadSequenceType
            self._seqs[item.start:item.stop:item.step] = value
        else:
            raise TypeError
    # iterators
    def __iter__(self):
        for item in self._seqs:
            yield item

class Alignment(SeqList):
    """
    Defines a class for an Alignment
    
    length of alignment when a sequence is added is
    enforced by adding gaps to the shortest sequence
    between those contained in the alignment and that
    being added.
    """
    def __init__(self, seq_cls=BaseSequence, seq_len=None, track=False):
        """
        track indica se cambiare la lunghezza di tutte le sequenze in base
        a quella appena aggiunta (ciuccia memoria che è una bellezza)
        """
        super(Alignment, self).__init__(seq_cls)
        if seq_len:
            self._seq_len = seq_len
        else:
            self._seq_len = 0
        self.track = track
    def _get_seq_len(self):
        """Gets current length of the alignment"""
        return self._seq_len
    def _set_seq_len(self, seq_len):
        """
        You can cut/expand the alignment changing these
        value
        """
        if self.track:
            if seq_len < self._seq_len:
                for seq in self:
                    del seq[seq_len:]
            elif seq_len > self._seq_len:
                for seq_cls in self:
                    seq_cls + '-'*(seq_len-self._seq_len)
        self._seq_len = seq_len
    _doc_seq_len = 'See doc for get/set length'
    def add_seq(self, name='', seq='', **kwds):
        """Additional checks for sequence length"""
        if len(seq) < self.seq_len:
            seq += '-'*(self.seq_len-len(seq))
        self.seq_len = len(seq)
        super(Alignment, self).add_seq(name, seq, **kwds)
    def remove_seqs(self):
        """Additional checks"""
        self._seq_len = 0
        super(Alignment, self).remove_seqs()
    #these applies to all sequences and accept slices objects
    def del_site(self, site):
        """Deletes a position of the alignment"""
        if len(self) == 0:
            return
        for seq in self:
            if isinstance(site, (int, long)):
                del seq[site]
            elif isinstance(site, slice):
                del seq[site.start:site.stop:site.step]
        self._seq_len = len(self[0])
    def chg_site(self, site, value='-'):
        """Change a site of the alignment to value"""
        if len(self) == 0:
            return
        for seq in self:
            if isinstance(site, (int, long)):
                seq[site] = value
            elif isinstance(site, slice):
                seq[site.start:site.stop:site.step] = value
    def add_site(self, site, value='-'):
        """Adds a site (defaults to '-') after specified position"""
        if len(self) == 0:
            return
        seq_len = self._seq_len
        for seq_cls in self:
            if isinstance(site, (int, long)):
                if site < seq_len:
                        seq_cls.seq.insert(site+1, value)
                elif site == seq_len:
                    seq_cls.seq.append('-')
            elif isinstance(site, slice):
                for pos in reversed(xrange(site.start, site.stop, site.step)):
                    seq_cls.seq.insert(pos+1, value)
        self._seq_len = len(self[0])

