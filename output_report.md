
# Results of filtering

We estimate the fraction of unigenes (sequences) kept after filtering the
[GMGCv1](https://gmgc.embl.de).

Using two datasets

1. Human gut dataset (from [Zeller et al., 2014](https://doi.org/10.15252/msb.20145645))
1. Dog gut dataset (from [Coelho et al., 2018](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0450-3))

Dataset size (number of samples):

| study            |   count |
|:-----------------|--------:|
| Coelho_dogs_2018 |      20 |
| Zeller_CRC_2014  |      20 |

## Fraction (%) of elements kept after filtering

### Using a minimum of 1 matched hash per unigene (for it to be kept)
|      |   Coelho_dogs_2018 |   Zeller_CRC_2014 |
|:-----|-------------------:|------------------:|
| mean |               4.21 |              5.83 |
| std  |               0.41 |              1.03 |
| min  |               3.61 |              2.69 |
| 25%  |               3.91 |              5.3  |
| 50%  |               4.21 |              5.94 |
| 75%  |               4.48 |              6.51 |
| max  |               4.98 |              7.13 |

### Using a minimum of 2 matched hashes per unigene (for it to be kept)
|      |   Coelho_dogs_2018 |   Zeller_CRC_2014 |
|:-----|-------------------:|------------------:|
| mean |               3.3  |              4.97 |
| std  |               0.33 |              0.92 |
| min  |               2.83 |              2.16 |
| 25%  |               3.06 |              4.53 |
| 50%  |               3.28 |              5.08 |
| 75%  |               3.48 |              5.58 |
| max  |               3.91 |              6.15 |


Note that the _fraction_ of kept unigenes is counted, which will not correspond
to the fraction of basepairs as longer genes are more likely to be kept.

