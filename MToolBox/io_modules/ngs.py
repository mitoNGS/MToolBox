def get_mapping_regions(contigs_file):
    """
    input: contig file with headings like
    >Contig.1|1-10000
    >Contig.1|11000-15000
    
    output: list of interval tuples
    [(1,10000), (11000,15000)]
    """
    contigs = open(contigs_file, 'r')
    j = []
    l = contigs.readline()
    while l:
        if l.startswith('>'):
            y = l.strip().split('|')[1].split('-')
            j.append((int(y[0]), int(y[1])))
        l = contigs.readline()
    return j