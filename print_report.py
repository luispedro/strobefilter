def print_report(final, size_fasta_sub, ofile):

    frac = final.dropna(axis=1)
    frac /= size_fasta_sub
    frac *= 100
    fra = frac.copy()
    frac['study'] = frac.index.str.split('/').str[0]

    table = frac.groupby('study').describe().round(2).T.loc['strict']

    study_counts =final.index.str.split('/').str[0].value_counts()
    study_counts.index.name = 'study'

    with open(ofile, 'w') as f:
        f.write(f'''
# Results of filtering

We estimate the fraction of unigenes (sequences) kept after filtering the
[GMGCv1](https://gmgc.embl.de).

Using two datasets

1. Human gut dataset (from [Zeller et al., 2014](https://doi.org/10.15252/msb.20145645))
1. Dog gut dataset (from [Coelho et al., 2018](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0450-3))

Dataset size (number of samples):

{study_counts.to_markdown()}

## Fraction (%) of elements kept after filtering

### Using a minimum of 1 matched hash per unigene (for it to be kept)
{table.loc['min1'].drop('count').to_markdown()}

### Using a minimum of 2 matched hashes per unigene (for it to be kept)
{table.loc['min2'].drop('count').to_markdown()}


Note that the _fraction_ of kept unigenes is counted, which will not correspond
to the fraction of basepairs as longer genes are more likely to be kept.

''')


