# Prefilter for mapping

**Problem:** Mapping metagenomes to large databases, such as the
[GMGCv1](https://gmgc.embl.de), takes too much memory. Partitioning the dataset
is a common solution (and [supported by
NGLess](https://ngless.embl.de/Mapping.html#low-memory-mode), but it has
drawbacks (it's slow).

This repository explores the possibility of prefiltering the database by
removing sequences that are extremely unlikely to be matches.

## Approach

1. Parse all the reads and collect all randstrobes (or rather their hashes)
2. Parse the database and select only unigenes that are expected to be present in the reads
3. Map as usual to the pre-filtered database

For `2`, different strategies are possible. The simplest is to keep any unigene
that shares any hash with the set of hashes from the reads. Currently being
considered

- `min1`: keep all references that match at least one hash
- `min2`: keep all references that match at least two hashes

We also tested counting the exact value or using a hacky Bloom filter structure
that uses a single fixed size array, but the hacky version gave bad estimates.

### Requirements

- Python, including NumPy and Pandas
- [Jug](https://jug.rtfd.io/)
- [NGLess](https://ngless.embl.de/)
- [Strobealign](https://github.com/ksahlin/strobealign) ([Sahlin, 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02831-7)), including the Python bindings
- [tabulate](https://pypi.org/project/tabulate/) is used to print the final table

To install most dependencies (assuming you have conda-forge & bioconda set up):

```
conda install python=3.11 numpy pandas requests tabulate jug ngless
```

To install `strobealign`'s Python bindings (which will **not** be installed by default with `conda`):

```
# To ensure you have a recent C++ compiler (not always needed)
conda install gxx_linux-64 gcc_linux-64
export CC CXX

git clone https://github.com/ksahlin/strobealign
cd strobealign
pip install .
```

### Data

1. _Database_ GMGCv1 (from ([Coelho et al., 2022](https://www.nature.com/articles/s41586-021-04233-4)). This can be is downloaded by `jugfile.py`
2. _Metagenomes_: Dog dataset (from [Coelho et al., 2018](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0450-3)) and human gut dataset (from [Zeller et al., 2014](https://doi.org/10.15252/msb.20145645). These can be downloaded with [ena-mirror](https://github.com/BigDataBiology/ena-mirror). More guidance will be provided on how to do it soon, but [get in touch](https://github.com/luispedro/strobefilter/issues) if you have questions.

Note that running this benchmark will use a lot of disk storage!


### Author

- [Luis Pedro Coelho](https://luispedro.org) (Queensland University of Technology). [luis@luispedro.org](mailto:luis@luispedro.org)
