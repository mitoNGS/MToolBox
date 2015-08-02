# coding=utf8
"""
Non mi ricordo neanche QUANTO funzioni evitare in ogni caso
"""

import os, time, datetime
import re #amiche mie... salvatemi voi!

FILE_NCBI_FLAT = 'ncbi_flat'
FILETYPE = 1 #e' una vera e propria entry che contiene molte informazioni

def check_ncbi_flat(fname):
    """
    Rudimentale controllo, vede se i primi 5 byte sono LOCUS
    
    @param fname:  il nome del file
    
    @raise ValueError: nel caso il file non esista
    """
    if not os.path.isfile(fname):
        raise ValueError(ERR_FILE_NOT_EXIST % fname)
    f = open(fname, 'U')
    if f.read(5) == 'LOCUS':
        f.close()
        return True
    else:
        f.close()
        return False

locus = re.compile("""LOCUS\s+
        ([\w\d]+)\s+         #locus name
        ([\d]+)\s            #length of sequence (not sure if a single space is enough)
        (bp|aa)\s+           #bp or aa (aa is my guess)
        (DNA|AA)\s+          #sequence type (AA is my guess)
        (circular|linear)\s+ #type of molecule (linear is my guess)
        (PRI|PUB|PAT)\s+         #my guess: private sequence (if so PUB should be public)
        (\d{2}-\w{3}-\d{4})  #date of submission?
            """, re.IGNORECASE | re.VERBOSE)
definition = re.compile("DEFINITION\s+(.+)", re.VERBOSE|re.IGNORECASE)
accession = re.compile("ACCESSION\s+(.+)", re.VERBOSE|re.IGNORECASE)
version = re.compile("""VERSION\s+
        [\w\d]+\.(\d{1,2})\s+ #accession number including version number
        (GI:\d+)              #maybe an index?
        """, re.VERBOSE|re.IGNORECASE)
keywords = re.compile("KEYWORDS\s+(.+)", re.VERBOSE|re.IGNORECASE)
#facenti parte del campo SOURCE
source = re.compile("SOURCE\s+(.+)", re.VERBOSE)
organism = re.compile("\s+ORGANISM\s+(.+)", re.VERBOSE|re.IGNORECASE)
#facenti parte del campo REFERENCE
reference = re.compile("""REFERENCE\s+
        (\d{1,2})\s+       #index of reference?
        \((bases|aa)\s     #bases (my guess aa)
        (\d+)\sto\s(\d+)\) #range of sequence covered by reference
        """, re.VERBOSE|re.IGNORECASE)
authors = re.compile("\s+AUTHORS\s+(.+)", re.VERBOSE|re.IGNORECASE)
title = re.compile("\s+TITLE\s+(.+)", re.VERBOSE|re.IGNORECASE)
journal = re.compile("\s+JOURNAL\s+(.+)", re.VERBOSE|re.IGNORECASE)
pubmed = re.compile("\s+PUBMED\s+(\d+)", re.VERBOSE|re.IGNORECASE)
#facenti parte del campo FEATURES
features = re.compile("FEATURES\s+Location/Qualifiers", re.VERBOSE|re.IGNORECASE)
f_genes = re.compile("""
        (source|D-loop|rRNA|tRNA|gene|CDS|variation)\s+  #features
        (complement)?\(?                       #case of gene on complementary strand
        (join)?\(?                             #case in which a gene is compromised of more than one region
        ([.\d,]+)\)?\)?                        #takes the whole series of ranges (to be processed)
        """, re.VERBOSE|re.IGNORECASE)
f_info = re.compile("""
        /(.+)     # name of field
        ="?       #quotes are used for a string value (ignored)
        ([^"\n]+) #argument 
        "?
        """, re.VERBOSE|re.IGNORECASE)
#campo ORIGIN
origin = re.compile("ORIGIN\s+", re.VERBOSE|re.IGNORECASE)

root_rex   = {'locus': locus, 'definition': definition, 'accession': accession, 'version': version,
              'keywords': keywords, 'source': source, 'reference': reference, 'features': features,
              'origin': origin}
source_rex = {'organism': organism}
ref_rex    = {'authors': authors, 'title': title, 'journal': journal, 'pubmed': pubmed}
feat_rex   = {'genes': f_genes, 'info': f_info}

def load_ncbi_flat(fname, addFunc):
    """
    Carica un file ncbi flat entry
    
    @type fname: string
    @param fname: nome del file
    @type addFunc: function
    @param addFunc: puntatore alla funzione di aggiunta
    
    @raise ValueError: nel caso il file non esista
    """
    if not os.path.isfile(fname):
        raise ValueError(ERR_FILE_NOT_EXIST % fname)
    f = open(fname, 'U')
    in_source = in_ref = in_feat = in_origin = False
    last_field = ''
    last_feat = None
    ref_info = []
    feat_info = []
    source_info = []
    origin_info = []
    seq_info = {'source': source_info, 'reference': ref_info, 'features': feat_info, 'origin': origin_info}
    for line in f:
        match = None
        field = None
        for key, rex in root_rex.iteritems():
            match = rex.search(line)
            if match != None:
                field = key
                in_source = in_ref = in_feat = in_origin = False
                last_field = ''
                break
        if field == 'source':
            in_source = True
            source_info.append( { field: match.group(1)} )
            last_feat = source_info[-1]
        elif field == 'reference':
            in_ref = True
            ref_info.append( {'index': int(match.group(1)), 'covers': (int(match.group(3)), int(match.group(4)))} )
            last_feat = ref_info[-1]
        elif field == 'features':
            in_feat = True
            last_feat = last_field = None
        elif field == 'origin':
            in_origin = True
        #non vanno su altre righe
        elif field == 'locus':
            #convertire la data in oggetto datetime e controllare che fare di 'bp'
            seq_info[field] = {'name': match.group(1), 'length': int(match.group(2)),
                               'seqtype': match.group(4), 'molecule': match.group(5),
                               'date': datetime.date(*(time.strptime(match.group(7),'%d-%b-%Y')[:3]))}
        elif field == 'definition':
            seq_info[field] = match.group(1)
        elif field == 'accession':
            seq_info[field] = match.group(1)
        elif field == 'version':
            seq_info[field] = {'version': int(match.group(1)), 'other': match.group(2)}
        elif field == 'keywords':
            seq_info[field] = match.group(1)
        else:
            if line[:2].strip() != '//':
                if in_source:
                    field = match = None
                    for key, rex in source_rex.iteritems():
                        match = rex.search(line)
                        if match != None:
                            field = key
                            break
                    if field == 'organism':
                        last_feat[field] = match.group(1)
                    else:
                        last_field = 'taxonomy'
                        try:
                            last_feat[last_field].append(line.strip())
                        except KeyError:
                            last_feat[last_field] = [line.strip()]
                elif in_ref:
                    field = match = None
                    for key, rex in ref_rex.iteritems():
                        match = rex.search(line)
                        if match != None:
                            field = key
                            break
                    if field == 'authors':
                        last_field = 'authors'
                        last_feat[last_field] = [match.group(1)]
                    elif field == 'title':
                        last_field = 'title'
                        last_feat[last_field] = [match.group(1)]
                    elif field == 'journal':
                        last_field = 'journal'
                        last_feat[last_field] = [match.group(1)]
                    elif field == 'pubmed':
                        last_field = 'pubmed'
                        last_feat[last_field] = int(match.group(1))
                    else:
                        last_feat[last_field].append(line.strip())
                elif in_feat:
                    field = match = None
                    for key, rex in feat_rex.iteritems():
                        match = rex.search(line)
                        if match != None:
                            field = key
                            break
                    if field == 'genes':
                        pos = match.group(4).split(',')
                        feat_info.append( {'pos':pos, 'type':match.group(1), 'complement': True if match.group(2) else False})
                        last_feat = feat_info[-1]
                    elif field == 'info':
                        if match.group(2).isdigit():
                            value = int(match.group(2))
                        elif match.group(1) in ['translation', 'note']:
                            value = [match.group(2)]
                        else:
                            value = match.group(2)
                        last_feat[match.group(1)] =  value
                        last_field = match.group(1)
                    else:
                        last_feat[last_field].append(line.strip(' \n\t"'))
                elif in_origin:
                    origin_info.append(''.join(line.split()[1:]))
            else:
                addFunc(filetype=FILE_NCBI_FLAT, info=seq_info)
                in_source = in_ref = in_feat = in_origin = False
                last_field = ''
                last_feat = None
                ref_info = []
                feat_info = []
                source_info = []
                origin_info = []
                seq_info = {'source': source_info, 'reference': ref_info, 'features': feat_info, 'origin': origin_info}
