# Prefilter for mapping

**Problem:** Mapping metagenomes to large databases, such as the
[GMGCv1](https://gmgc.embl.de), takes too much memory. Partioning the dataset
is a common solution, but it has drawbacks.

This repository explores prefiltering the database by removing sequences that
are extremely unlikely to be matches.

## Approach

1. Parse all the reads and collect all randstrobes (or rather their hashes)
2. Parse the database and select only unigenes that are expected to be present in the reads
3. Map as usual to the pre-filtered database

For `2`, different strategies are possible. The simplest is to keep any unigene
that shares any hash with the set of hashes from the reads. Currently being considered

- `strict`: keep all references that match at least one hash
- `min2`: keep all references that match at least two hashes
- `packed`: a hacky Bloom filter structure that uses a single fixed size array
  (break up the 64 bit hash into two 32 bit hashes and set both to `1`)


### Requirements

- Python
- [Jug](https://jug.rtfd.io/)
- [NGLess](https://ngless.embl.de/)
- [Strobealign](https://github.com/ksahlin/strobealign) ([Sahlin, 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02831-7)), including the Python bindings

### Data

1. GMGCv1 (from ([Coelho et al., 2022](https://www.nature.com/articles/s41586-021-04233-4)). This can be is downloaded by `jugfile.py`

### Author

- [Luis Pedro Coelho](https://luispedro.org) (Queensland University of Technology). [luis@luispedro.org](mailto:luis@luispedro.org)
