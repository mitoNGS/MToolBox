"""
Utilities for VCF files.
"""


def walk_together(*readers):
    """ Simultaneously iteratate two or more VCF readers and return
        lists of concurrent records from each
        reader, with None if no record present.  Caller must check the
        inputs are sorted in the same way and use the same reference
        otherwise behaviour is undefined.
    """
    # if one of the VCFs has no records, StopIteration is
    # raised immediately, so we need to check for that and
    # deal appropriately
    nexts = []
    for reader in readers:
        try:
            nexts.append(reader.next())
        except StopIteration:
            nexts.append(None)

    while True:
        min_next = min([x for x in nexts if x is not None])

        # this line uses equality on Records, which checks the ALTs
        # not sure what to do with records that have overlapping but different
        # variation
        yield [x if x is None or x == min_next else None for x in nexts]

        # update nexts that we just yielded
        for i, n in enumerate(nexts):

            if n is not None and n == min_next:
                try:
                    nexts[i] = readers[i].next()
                except StopIteration:
                    nexts[i] = None

        if all([x is None for x in nexts]):
            break


def trim_common_suffix(*sequences):
    """
    Trim a list of sequences by removing the longest common suffix while
    leaving all of them at least one character in length.

    Standard convention with VCF is to place an indel at the left-most
    position, but some tools add additional context to the right of the
    sequences (e.g. samtools). These common suffixes are undesirable when
    comparing variants, for example in variant databases.

        >>> trim_common_suffix('TATATATA', 'TATATA')
        ['TAT', 'T']

        >>> trim_common_suffix('ACCCCC', 'ACCCCCCCC', 'ACCCCCCC', 'ACCCCCCCCC')
        ['A', 'ACCC', 'ACC', 'ACCCC']

    """
    if not sequences:
        return []
    reverses = [seq[::-1] for seq in sequences]
    rev_min = min(reverses)
    rev_max = max(reverses)
    if len(rev_min) < 2:
        return sequences
    for i, c in enumerate(rev_min[:-1]):
        if c != rev_max[i]:
            if i == 0:
                return sequences
            return [seq[:-i] for seq in sequences]
    return [seq[:-(i + 1)] for seq in sequences]
