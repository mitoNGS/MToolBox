from model import _Call

cdef _map(func, iterable, bad='.'):
    '''``map``, but make bad values None.'''
    return [func(x) if x != bad else None
            for x in iterable]

INTEGER = 'Integer'
FLOAT = 'Float'
NUMERIC = 'Numeric'

def parse_samples(
        list names, list samples, samp_fmt,
        list samp_fmt_types, list samp_fmt_nums, site):

    cdef char *name, *fmt, *entry_type, *sample
    cdef int i, j
    cdef list samp_data = []
    cdef dict sampdict
    cdef list sampvals
    n_samples = len(samples)
    n_formats = len(samp_fmt._fields)

    for i in range(n_samples):
        name = names[i]
        sample = samples[i]

        # parse the data for this sample
        sampdat = [None] * n_formats

        sampvals = sample.split(':')

        for j in range(n_formats):
            if j >= len(sampvals):
                break
            vals = sampvals[j]

            # short circuit the most common
            if vals == '.' or vals == './.':
                sampdat[j] = None
                continue

            entry_type = samp_fmt_types[j]
            # TODO: entry_num is None for unbounded lists
            entry_num = samp_fmt_nums[j]

            # we don't need to split single entries
            if entry_num == 1 or ',' not in vals:

                if entry_type == INTEGER:
                    sampdat[j] = int(vals)
                elif entry_type == FLOAT or entry_type == NUMERIC:
                    sampdat[j] = float(vals)
                else:
                    sampdat[j] = vals

                if entry_num != 1:
                    sampdat[j] = (sampdat[j])

                continue

            vals = vals.split(',')

            if entry_type == INTEGER:
                sampdat[j] = _map(int, vals)
            elif entry_type == FLOAT or entry_type == NUMERIC:
                sampdat[j] = _map(float, vals)
            else:
                sampdat[j] = vals

        # create a call object
        call = _Call(site, name, samp_fmt(*sampdat))
        samp_data.append(call)

    return samp_data
