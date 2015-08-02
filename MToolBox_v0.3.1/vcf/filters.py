try:
    from rpy2 import robjects
except:
    robjects = None


class Base(object):
    """ Base class for vcf_filter.py filters.

        Use the class docstring to provide the filter description
        as it appears in vcf_filter.py
    """

    name = 'f'
    """ name used to activate filter and in VCF headers """

    @classmethod
    def customize_parser(self, parser):
        """ hook to extend argparse parser with custom arguments """
        pass

    def __init__(self, args):
        """ create the filter using argparse ``args`` """
        self.threshold = 0

    def __call__(self):
        """ filter a site, return not None if the site should be filtered """
        raise NotImplementedError('Filters must implement this method')

    def filter_name(self):
        """ return the name to put in the VCF header, default is ``name`` + ``threshold`` """
        return '%s%s' % (self.name, self.threshold)


class SiteQuality(Base):
    """ Filter low quailty sites """

    name = 'sq'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--site-quality', type=int, default=30,
                help='Filter sites below this quality')

    def __init__(self, args):
        self.threshold = args.site_quality

    def __call__(self, record):
        if record.QUAL < self.threshold:
            return record.QUAL


class VariantGenotypeQuality(Base):
    """ Filters sites with only low quality variants.

        It is possible to have a high site quality with many low quality calls.  This
        filter demands at least one call be above a threshold quality.
    """

    name = 'mgq'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--genotype-quality', type=int, default=50,
                help='Filter sites with no genotypes above this quality')

        def __init__(self, args):
            self.threshold = args.genotype_quality

    def __call__(self, record):
        if not record.is_monomorphic:
            vgq = max([x['GQ'] for x in record if x.is_variant])
            if vgq < self.threshold:
                return vgq


class ErrorBiasFilter(Base):
    """ Filter sites that look like correlated sequencing errors.

        Some sequencing technologies, notably pyrosequencing, produce mutation
        hotspots where there is a constant level of noise, producing some reference
        and some heterozygote calls.

        This filter computes a Bayes Factor for each site by comparing
        the binomial likelihood of the observed allelic depths under:

        * A model with constant error equal to the MAF.
        * A model where each sample is the ploidy reported by the caller.

        The test value is the log of the bayes factor.  Higher values
        are more likely to be errors.

        Note: this filter requires rpy2
    """

    name = 'eb'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--eblr', type=int, default=-10,
                help='Filter sites above this error log odds ratio')

    def __init__(self, args):
        self.threshold = args.eblr
        if robjects is None:
            raise Exception('Please install rpy2')
        self.ll_test = robjects.r('''
            function(ra, aa, gt, diag=F) {
                ra_sum = sum(ra)
                aa_sum = sum(aa)
                ab = aa_sum / (ra_sum + aa_sum)
                gtp = 0.5 + (0.48*(gt-1))

                error_likelihood = log(dbinom(aa, ra+aa, ab))
                gt_likelihood = log(dbinom(aa, ra+aa, gtp))

                if (diag) {
                    print(ra)
                    print(aa)
                    print(gtp)
                    print(ab)
                    print(error_likelihood)
                    print(gt_likelihood)
                }
                error_likelihood = sum(error_likelihood)
                gt_likelihood = sum(gt_likelihood)
                c(error_likelihood - gt_likelihood, ab)
            }
            ''')

    def __call__(self, record):
        if record.is_monomorphic:
            return None
        passed, tv, ab = self.bias_test(record.samples)
        if tv > self.threshold:
            return tv

    def bias_test(self, calls):
        calls = [x for x in calls if x.called]
        #TODO: single genotype assumption
  
        try:
            # freebayes
            ra = robjects.IntVector([x['RO'][0] for x in calls])
            aa = robjects.IntVector([x['AO'][0] for x in calls])
        except AttributeError:
            # GATK
            ra = robjects.IntVector([x['AD'][0] for x in calls])
            aa = robjects.IntVector([x['AD'][1] for x in calls])

        gt = robjects.IntVector([x.gt_type for x in calls])
        test_val, ab = self.ll_test(ra, aa, gt)

        return test_val < 0, test_val, ab


class DepthPerSample(Base):
    'Threshold read depth per sample'

    name = 'dps'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--depth-per-sample', type=int, default=5,
                help='Minimum required coverage in each sample')

    def __init__(self, args):
        self.threshold = args.depth_per_sample

    def __call__(self, record):
        # do not test depth for indels
        if record.is_indel:
            return

        mindepth = min([sam['DP'] for sam in record.samples])
        if mindepth < self.threshold:
            return mindepth


class AvgDepthPerSample(Base):
    'Threshold average read depth per sample (read_depth / sample_count)'

    name = 'avg-dps'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--avg-depth-per-sample', type=int, default=3,
              help='Minimum required average coverage per sample')

    def __init__(self, args):
        self.threshold = args.avg_depth_per_sample

    def __call__(self, record):
        avgcov = float(record.INFO['DP']) / len(record.samples)
        if avgcov < self.threshold:
            return avgcov


class SnpOnly(Base):
    'Choose only SNP variants'

    name = 'snp-only'

    def __call__(self, record):
        if not record.is_snp:
            return True

    def filter_name(self):
        return self.name
