# Prefilter for mapping

## prefiltered mapping

1. Parse all the reads and collect all randstrobes (or rather their hashes)
2. Parse the database and select only unigenes that are expected to be present in the reads
3. Map as usual to the pre-filtered database

For `2`, different strategies are possible. The simplest is to keep any unigene that shares any hash with the set of hashes from the reads.

### Requirements

- Python
- [Jug](https://jug.rtfd.io/)
- [NGLess](https://ngless.embl.de/)

### Data

1. GMGCv1 (from ([Coelho et al., 2022](https://www.nature.com/articles/s41586-021-04233-4)). This can be is downloaded by `jugfile.py`
