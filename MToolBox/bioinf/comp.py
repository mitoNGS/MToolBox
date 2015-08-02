# coding=utf8
"""
Idealmente dovremmo usare una classe che dato un allineamento nel calcola
la composizione (o usarla come mixin?)

Il perche' deriva da situazioni come il sitevar in cui va calcolata su regioni
dell'allineamento e solo su alcune di queste sequenze. L'alternativa sarebbe di
utilizzare differenti allineamenti, ma non mi piace.

La modifica da fare alla classe dell'allineamento sarebbe di settare direttamete
la classe della composizione opportuna (aa o dna) che si occupa della
composizione e quindi renderla piu' neutrale

Magari fare delle classi che gia' inizializzano la composizione che serve (fatta
una) e mettere dei check sull'allineamento che viene dato

Done:
    - L'oggetto di base (o quello che serve e basta?)
    - Riadattata la roba per calcolare la composizione
    - Amminoacidi
    - Controllare il tutto e confrontare il calcolo con quanto fatto dalla
      classe Alignment (Fatte prove e confrontati i risultati
    - La proprieta' comp restituisce una lista con la composizione
    - E' possibile iterare la composizione

ToDo:
    - Commentare e sistemare un po' tutto
"""

__deprecated__ = "Composition.getCompositionRange Composition.calcComp"

from utils import autoprop
import seqs
from math import sqrt

#: dizionario contenente i tipi di composizione usati
comp_types = {'dna':
                {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'U': 0, '-': 0,
                 'R': 0, 'Y': 0, 'K': 0, 'M': 0, 'S': 0, 'W': 0,
                 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0
                },
              'aa':
                {
                 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0,
                 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'O': 0, 'P': 0, 'Q': 0, 'R': 0,
                 'S': 0, 'T': 0, 'U': 0, 'V': 0, 'W': 0, 'Y': 0, '-': 0
                }
             }

class Composition(object):
    """
    Classe che calcola/contiene la composizione di un allineamento.
    
    Gli oggetti supportati sono L{Alignment} e L{SeqList}.
    
    @attention: L{SeqList} DEVE avere sequenze tuute della stessa lunghezza  
    """
    __metaclass__ = autoprop #: la mia bella metaclasse (scopiazzata)
    
    #: il tipo di composizione da usare
    _comp_type = comp_types['dna']
    def _set_comp_type(self, type='dna'):
        """
        Imposta il tipo di composizione, usare una di quelle definite da
        L{comp_types} oppure utilizzarne un'altra conforme a quelle specificate
        
        @type type: string
        @keyword type: il tipo di composizione, B{dna} o B{aa} sono i due tipi
                       supportati
        """
        self._comp_type = comp_types[type]
    _alg = None #: l'allineamento di cui calcolare la composizione
    def _set_alg(self, alg):
        """
        Imposta l'allineamento da utilizzare
        
        @param alg: L{Alignment} o L{SeqList}
        
        @warning: non ci sono controlli sulla validita' di alg al momento
        """
        self._alg = alg
    def _get_alg(self):
        """
        Restituisce l'allineamento utilizzato
        """
        return self._alg
    
    _starts = None #: lista/tupla delle posizioni d'inizio
    _ends = None #: lista/tupla delle posizioni di fine
    _seq_lists = None #: lista/tupla contenente gli elenchi di sequenze
    def set_ranges(self, starts=None, ends=None, seq_lists=None):
        """
        Imposta i la/le regione/i su cui calcolare la composizione e le relative
        sequenze.
        
        starts, ends, seq_lists B{must be iterables} as must be each element of
        seq_lists (a list sequence which will be used in the calculation each
        tuple (start, end) B{must be the same as used in a slice of a sequence}
        
        @attention: i valori di starts ed end devono essere in ordine crescente!!!
        
        Esempio:
        C{starts = ( 0,  40, 80 )}
        C{ends   = ( 40, 80, 120 )}
        C{seq_lists = ( (0, 1, 2, 3), (2, 3), (0, 1, 2, 3) )}
        
        In tabella::
            start |  end |  seq_list
            -------------------------
            0     |  40  |  (0, 1, 2, 3)
            40    |  80  |  (2, 3)
            80    | 120  |  (0, 1, 2, 3)
        
        Nel loop di calcolo succedera' questo nell'ordine:
            1. per lo slice [0:40]:
                - verra' calcolata la composizione per le sequenze di indice
                  0, 1, 2 e 3
            2. per lo slice [40:80]:
                - verra' calcolata la composizione per le sequenze di indice
                  2 e 3
            3. per lo slice [80:120]:
                - verra' calcolata la composizione per le sequenze di indice
                  0, 1, 2 e 3
        
        In generale si puo' evitare di usarla (o usarla senza parametri) dato che:
            1. in generale interessa la composizione dell'intero allineamento
            2. nel richiamare L{calc_comp}, viene controllato che i parametri necessari
               per il calcolo siano impostati ed imposta automaticamente i valori
            3. anche se si richiamasse L{get_comp_range} o la proprieta' L{comp}, la
               la composizione verrebbe comunque calcolata prima
        
        Possibili default "misti":
            - Se starts o ends sono vuoti vengo settati ad una solo regione
              (allineamento intero)
            - Se seq_lists e' vuoto si assume che si voglia fare il calcolo su tutto
              l'allineamento
        Quindi e' possibile impostare solo un "gruppo" di parametri se si vuole
        
        @type starts: iterable
        @keyword starts: sequenza che resituisca i valori in di inizio
        @type ends: iterable
        @keyword ends: sequenza che restituisca i valori di fine
        @type seq_lists: iterable
        @keyword seq_lists: sequenza che restituisca una sequenza per ogni suo
                            elemento, quest'ultima contenente gli indici delle
                            sequenze su cui calcolare la composizione
        """
        #se entrambi sono None (o vuoti) allora setta dei valori totali
        if not (starts and ends):
            self._starts = (0,)
            self._ends = ( self._alg.seq_len, )
        else:
            self._starts = starts
            self._ends = ends
        #se non sono state specificate le sequenze su cui agire, si sottintende
        #tutte quelle dell'allineamento
        if not seq_lists:
            if starts:
                self._seq_lists = ( range(len(self._alg)), ) * len(starts)
            else:
                self._seq_lists = ( range(len(self._alg)), )
        else:
            self._seq_lists = seq_lists
            
    def calc_comp(self):
        """
        Calcola la composizione utilizzando i parametri impostati, viene chiamata
        L{set_ranges} per l'intero allineamento se almeno un parametro non e' stato
        specificato
        """
        if not self._alg:
            return
        if not (self._starts and self._ends and self._seq_lists):
            self.set_ranges()
        #mandiamo il calcolo
        self._calc_comp()
    calcComp = calc_comp
    
    def _calc_comp(self):
        """
        Il vero e proprio calcolo, nel caso mi venisse il piccio di cambiare
        l'algoritmo, mi basterebbe ridefinire il puntatore a questa funzione
        
        Per il momento va bene cosi', ma la si potrebbe sparare in C se
        necessario, bisogna vedere che punto (forse solo i due loop interni)
        """
        starts = self._starts
        ends = self._ends
        seq_lists = self._seq_lists
        comp_type = self._comp_type
        self._comp = []
        for start, end, seq_list in zip(starts, ends, seq_lists):
            #costruisce una lista di valori di composizione vuoti
            comp = [comp_type.copy() for idx in range(end-start)]
            list_len = len(seq_list)
            for idx in seq_list:
                seq = self._alg[idx][start:end]
                for pos, item in enumerate(seq):
                    if item in comp_type:
                        comp[pos][item] += 1.0 / list_len
            self._comp += comp
        
    def _check_comp(self):
        """Controlla che sia stata calcolata la composizione"""
        if not hasattr(self, '_comp'):
            self.calc_comp()
        
    def _get_comp(self):
        """
        Ritorna una lista contenente l'allineamento
        
        @rtype: list
        @return: la composizione dell'allineamento in una lista
        """
        self._check_comp()
        return [site for site in self]
    _doc_comp = 'Get Alignment composition'
    
    def get_comp_range(self, spos=None, epos=None):
        """
        Ritorna solo la composizione dell'allineamento tra le posizioni specificate
        
        @attention: come per la maggior parte delle funzioni che implementano
            il ritorno di un range, viene utilizzata la notazione degli slice
            in python ovvero [spos:epos]
        
        @type spos: int
        @param spos: estremo inferiore
        @type epos: int
        @param epos: estremo superiore
        
        @rtype: generator
        @return: la composizione dell'allineamento
        """
        #la vecchia regex prevedeva dei controlli spos ed epos se fossero uguali a zero
        #la situazione cambia leggermente se si usano slice ed i valori direttamente
        #per cui se sono uguali a 0 si vuole l'intero range che viene dato con [None:None]
        self._check_comp()
        if spos == epos == 0:
            spos = epos = None
        for site in self._comp[spos:epos]:
            yield site
    getCompositionRange = get_comp_range

    def __iter__(self):
        """
        Differisce dalla proprieta' L{comp} per l'utilizzo di un generatore, invece
        di una intera lista.
        
        @rtype: generator
        @return: la composizione dell'allineamento
        """
        self._check_comp()
        for site in self._comp:
            yield site
    
    def __len__(self):
        """
        Ritorna la lunghezza della composizione
        @rtype: int
        @return: lunghezza della composizione
        """
        if getattr(self, '_comp', None):
            return len(self._comp)
        else:
            return 0
    # indexing and slicing
    def __getitem__(self, item):
        if isinstance(item, (int, long)):
            return self._comp[item]
        elif isinstance(item, slice):
            return self._comp[item.start:item.stop:item.step]
        else:
            raise TypeError

class NucComposition(Composition):
    """
    Classe che contiene l'inizializzazione del tipo di composizione a B{dna}
    """
    def __init__(self):
        """L'inizializzazione viene lasciata che non si sa mai"""
        self._set_comp_type('dna')
    
class AaComposition(Composition):
    """
    Classe che contiene l'inizializzazione del tipo di composizione a B{aa}
    """
    def __init__(self):
        self._set_comp_type('aa')

class Statistics(object):
    """
    Accetta un oggetto che supporta il protocollo di iterazione, e a ogni
    iterazione restituirà un oggetto L{Sequence} oppure una stringa (che verrà
    comunque convertita in un oggetto L{Sequence} per calcolare la composizione)
    """
    __metaclass__ = autoprop #: la mia bella metaclasse (scopiazzata)
    _stats = {}
    def __init__(self, seqs=None, seq_type='dna'):
        self._seqs = seqs
        self._seq_type = seq_type
    def _set_comp_type(self, comp_type):
        self._comp_type = comp_type
    def _get_comp_type(self):
        return self._comp_type
    def _set_seqs(self, seqs):
        self._seqs = seqs
    def _get_seqs(self):
        return self._seqs
    def calc_len(self):
        self.calc_avg_len()
        self.calc_dev_len()
    def calc_avg_len(self):
        """
        Lunghezza media delle sequenze
        """
        average = 0
        for seq in self._seqs:
            average += len(seq)
        average = average / float(len(self._seqs))
        self._stats['avg_len'] = average
    def calc_dev_len(self):
        """
        Calcolo della deviazione standard dalla lunghezza media
        """
        avg = self._stats.get('avg_len', None)
        if not avg:
            self.calc_avg_len()
            avg = self._stats['avg_len']
        avg = float(avg)
        #la sommatoria
        values = [ (len(x)-avg)**2 for x in self._seqs ]
        #il resto dei calcoli
        std_dev = sqrt(sum(values) / len(self._seqs))
        self._stats['dev_len'] = std_dev
    def calc_comp(self):
        self.calc_avg_comp()
        self.calc_dev_comp()
    def calc_avg_comp(self):
        """
        Composizione media dei nucleotidi
        """
        comp = comp_types[self._seq_type].copy()
        for seq in self._seqs:
            if not hasattr(seq, 'calc_comp'):
                seq = seqs.Sequence('x', seq, self._seq_type)
            tmp = seq.comp
            length = float(len(seq))
            for c in tmp:
                if tmp[c] > 0:
                    comp[c] += tmp[c] / length
        length = len(self._seqs)
        for c in comp:
            if comp[c] > 0:
                comp[c] = comp[c] / length
        self._stats['avg_comp'] = comp
    def calc_dev_comp(self):
        """
        calcolo della deviazione standard della composizione nucleotidica
        """
        comp = comp_types[self._seq_type].copy()
        avg = self._stats.get('avg_comp', None)
        if not avg:
            self.calc_avg_comp()
        for seq in self._seqs:
            if not hasattr(seq, 'calc_comp'):
                seq = seqs.Sequence('x', seq, self._seq_type)
            tmp = seq.comp
            length = float(len(seq))
            for c in tmp:
                if tmp[c] > 0:
                    comp[c] += ((tmp[c] / length) - avg[c]) ** 2
        for c in comp:
            if comp[c] > 0:
                comp[c] = sqrt(comp[c] / len(self._seqs))
        self._stats['dev_comp'] = comp
    def calc_min_len(self):
        self._stats['min_len'] = min([ len(x) for x in self._seqs ])
    def calc_max_len(self):
        self._stats['max_len'] = max([ len(x) for x in self._seqs ])
    def _get_stats(self):
        return self._stats
