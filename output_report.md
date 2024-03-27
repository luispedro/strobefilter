
# Results of filtering

We estimate the fraction of unigenes (sequences) kept after filtering the
[GMGCv1](https://gmgc.embl.de).

Using two datasets

1. Human gut dataset (from [Zeller et al., 2014](https://doi.org/10.15252/msb.20145645))
1. Dog gut dataset (from [Coelho et al., 2018](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0450-3))

Dataset size (number of samples):

| study                     |   count |
|:--------------------------|--------:|
| Coelho_dogs_2018          |      20 |
| Coelho_dogs_2018:dog-gut  |      20 |
| TaraOceans                |      20 |
| TaraOceans:marine         |      20 |
| Zeller_CRC_2014           |      20 |
| Zeller_CRC_2014:human-gut |      20 |

## Fraction (%) of elements kept after filtering

### Using a minimum of 1 matched hash per unigene (for it to be kept)
|      |   Coelho_dogs_2018 |   Coelho_dogs_2018:dog-gut |   TaraOceans |   TaraOceans:marine |   Zeller_CRC_2014 |   Zeller_CRC_2014:human-gut |
|:-----|-------------------:|---------------------------:|-------------:|--------------------:|------------------:|----------------------------:|
| mean |               4.21 |                      63.18 |         8.48 |               25.21 |              5.83 |                       28.04 |
| std  |               0.41 |                       7.31 |         4.5  |               14.39 |              1.03 |                        4.78 |
| min  |               3.61 |                      50.8  |         2.04 |                5.68 |              2.69 |                       12.97 |
| 25%  |               3.91 |                      57.88 |         3.13 |                7.71 |              5.3  |                       25.86 |
| 50%  |               4.21 |                      63.33 |        10    |               29.35 |              5.94 |                       28.73 |
| 75%  |               4.48 |                      68.59 |        12.07 |               35.48 |              6.51 |                       30.72 |
| max  |               4.98 |                      75.8  |        13.99 |               43.67 |              7.13 |                       34.06 |

### Using a minimum of 2 matched hashes per unigene (for it to be kept)
|      |   Coelho_dogs_2018 |   Coelho_dogs_2018:dog-gut |   TaraOceans |   TaraOceans:marine |   Zeller_CRC_2014 |   Zeller_CRC_2014:human-gut |
|:-----|-------------------:|---------------------------:|-------------:|--------------------:|------------------:|----------------------------:|
| mean |               3.3  |                      60.78 |         6.84 |               21.56 |              4.97 |                       25.6  |
| std  |               0.33 |                       7.54 |         3.89 |               12.89 |              0.92 |                        4.63 |
| min  |               2.83 |                      48.08 |         1.63 |                4.47 |              2.16 |                       11.12 |
| 25%  |               3.06 |                      55.35 |         2.1  |                5.81 |              4.53 |                       23.49 |
| 50%  |               3.28 |                      61.05 |         7.95 |               24.88 |              5.08 |                       26.28 |
| 75%  |               3.48 |                      66.3  |         9.94 |               30.84 |              5.58 |                       28.26 |
| max  |               3.91 |                      73.91 |        11.87 |               38.78 |              6.15 |                       31.51 |


Note that the _fraction_ of kept unigenes is counted, which will not correspond
to the fraction of basepairs as longer genes are more likely to be kept.

