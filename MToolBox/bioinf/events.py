# coding=utf8
pur = ('A','G','R')
pyr = ('T','C','Y')

amb = {
        'M':('A','C'), 'W':('A','T'), 'S':('C','G'),
        'K':('G','T'), 'V':('A','C','G'), 'H':('A','C','T'),
        'D':('A','G','T'), 'B':('C','G','T')
      }

def isTransition(ch1, ch2):
    """
    Checks if the event is a Transition
    
    ATTENZIONE: ritorna true anche nel caso che i due nucleotidi siano uguali
    si assume che sia stato gia' fatto un controllo
    """
    if ch1 in pur and ch2 in pur:
        return True
    elif ch1 in pyr and ch2 in pyr:
        return True
    return False

def isTransversion(ch1, ch2):
    """Checks if the event is a Transversion"""
    if ch1 in pyr and ch2 in pur:
        return True
    if ch1 in pur and ch2 in pyr:
        return True
    return False

def isAmbiguity(ch1, ch2):
    """Controlla se l'evento puo' essere ricondotto ad una ambiguita'"""
    if ch1 in amb:
        if ch2 in amb[ch1]: return True
    elif ch2 in amb:
        if ch1 in amb[ch2]: return True
    return False
