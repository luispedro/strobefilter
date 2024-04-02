
# Results of filtering

We estimate the fraction of unigenes (sequences) kept after filtering the
[GMGCv1](https://gmgc.embl.de).

Using two datasets

1. Human gut dataset (from [Zeller et al., 2014](https://doi.org/10.15252/msb.20145645))
1. Dog gut dataset (from [Coelho et al., 2018](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0450-3))

Dataset size (number of samples):

| study                              |   count |
|:-----------------------------------|--------:|
| Coelho_dogs_2018                   |      20 |
| Coelho_dogs_2018:dog-gut           |      20 |
| TaraOceans                         |      20 |
| TaraOceans:marine                  |      20 |
| Zeller_CRC_2014                    |      20 |
| Zeller_CRC_2014:human-gut          |      20 |
| Coelho_dogs_2018-combined          |       1 |
| Coelho_dogs_2018-combined:dog-gut  |       1 |
| TaraOceans-combined                |       1 |
| TaraOceans-combined:marine         |       1 |
| Zeller_CRC_2014-combined           |       1 |
| Zeller_CRC_2014-combined:human-gut |       1 |

## Fraction (%) of elements kept after filtering (with [NGLess](https://ngless.embl.de) preprocessing)

### Using a minimum of 1 matched hash per unigene (for it to be kept)
|      |   Coelho_dogs_2018 |   Coelho_dogs_2018-combined |   Coelho_dogs_2018-combined:dog-gut |   Coelho_dogs_2018:dog-gut |   TaraOceans |   TaraOceans-combined |   TaraOceans-combined:marine |   TaraOceans:marine |   Zeller_CRC_2014 |   Zeller_CRC_2014-combined |   Zeller_CRC_2014-combined:human-gut |   Zeller_CRC_2014:human-gut |
|:-----|-------------------:|----------------------------:|------------------------------------:|---------------------------:|-------------:|----------------------:|-----------------------------:|--------------------:|------------------:|---------------------------:|-------------------------------------:|----------------------------:|
| mean |               4.21 |                        7.36 |                               95.68 |                      63.18 |         8.48 |                 23.28 |                        62.13 |               25.21 |              5.83 |                      13.76 |                                59.24 |                       28.04 |
| std  |               0.41 |                      nan    |                              nan    |                       7.31 |         4.5  |                nan    |                       nan    |               14.39 |              1.03 |                     nan    |                               nan    |                        4.78 |
| min  |               3.61 |                        7.36 |                               95.68 |                      50.8  |         2.04 |                 23.28 |                        62.13 |                5.68 |              2.69 |                      13.76 |                                59.24 |                       12.97 |
| 25%  |               3.91 |                        7.36 |                               95.68 |                      57.88 |         3.13 |                 23.28 |                        62.13 |                7.71 |              5.3  |                      13.76 |                                59.24 |                       25.86 |
| 50%  |               4.21 |                        7.36 |                               95.68 |                      63.33 |        10    |                 23.28 |                        62.13 |               29.35 |              5.94 |                      13.76 |                                59.24 |                       28.73 |
| 75%  |               4.48 |                        7.36 |                               95.68 |                      68.59 |        12.07 |                 23.28 |                        62.13 |               35.48 |              6.51 |                      13.76 |                                59.24 |                       30.72 |
| max  |               4.98 |                        7.36 |                               95.68 |                      75.8  |        13.99 |                 23.28 |                        62.13 |               43.67 |              7.13 |                      13.76 |                                59.24 |                       34.06 |

### Using a minimum of 2 matched hashes per unigene (for it to be kept)
|      |   Coelho_dogs_2018 |   Coelho_dogs_2018-combined |   Coelho_dogs_2018-combined:dog-gut |   Coelho_dogs_2018:dog-gut |   TaraOceans |   TaraOceans-combined |   TaraOceans-combined:marine |   TaraOceans:marine |   Zeller_CRC_2014 |   Zeller_CRC_2014-combined |   Zeller_CRC_2014-combined:human-gut |   Zeller_CRC_2014:human-gut |
|:-----|-------------------:|----------------------------:|------------------------------------:|---------------------------:|-------------:|----------------------:|-----------------------------:|--------------------:|------------------:|---------------------------:|-------------------------------------:|----------------------------:|
| mean |               3.3  |                        5.88 |                               95.26 |                      60.78 |         6.84 |                 19.96 |                        57.34 |               21.56 |              4.97 |                      12.11 |                                56.56 |                       25.6  |
| std  |               0.33 |                      nan    |                              nan    |                       7.54 |         3.89 |                nan    |                       nan    |               12.89 |              0.92 |                     nan    |                               nan    |                        4.63 |
| min  |               2.83 |                        5.88 |                               95.26 |                      48.08 |         1.63 |                 19.96 |                        57.34 |                4.47 |              2.16 |                      12.11 |                                56.56 |                       11.12 |
| 25%  |               3.06 |                        5.88 |                               95.26 |                      55.35 |         2.1  |                 19.96 |                        57.34 |                5.81 |              4.53 |                      12.11 |                                56.56 |                       23.49 |
| 50%  |               3.28 |                        5.88 |                               95.26 |                      61.05 |         7.95 |                 19.96 |                        57.34 |               24.88 |              5.08 |                      12.11 |                                56.56 |                       26.28 |
| 75%  |               3.48 |                        5.88 |                               95.26 |                      66.3  |         9.94 |                 19.96 |                        57.34 |               30.84 |              5.58 |                      12.11 |                                56.56 |                       28.26 |
| max  |               3.91 |                        5.88 |                               95.26 |                      73.91 |        11.87 |                 19.96 |                        57.34 |               38.78 |              6.15 |                      12.11 |                                56.56 |                       31.51 |


## Fraction (%) of elements kept after filtering (passthru, i.e., no FQ preprocessing)

### Using a minimum of 1 matched hash per unigene (for it to be kept)
|      |   Coelho_dogs_2018 |   Coelho_dogs_2018-combined |   Coelho_dogs_2018-combined:dog-gut |   Coelho_dogs_2018:dog-gut |   TaraOceans |   TaraOceans-combined |   TaraOceans-combined:marine |   TaraOceans:marine |   Zeller_CRC_2014 |   Zeller_CRC_2014-combined |   Zeller_CRC_2014-combined:human-gut |   Zeller_CRC_2014:human-gut |
|:-----|-------------------:|----------------------------:|------------------------------------:|---------------------------:|-------------:|----------------------:|-----------------------------:|--------------------:|------------------:|---------------------------:|-------------------------------------:|----------------------------:|
| mean |               4.31 |                        7.63 |                               95.87 |                      63.76 |         9.08 |                 24.69 |                        63.68 |               26.29 |              5.89 |                      13.84 |                                59.44 |                       28.29 |
| std  |               0.42 |                      nan    |                              nan    |                       7.36 |         4.58 |                nan    |                       nan    |               14.58 |              1.03 |                     nan    |                               nan    |                        4.77 |
| min  |               3.73 |                        7.63 |                               95.87 |                      51.47 |         2.34 |                 24.69 |                        63.68 |                6.49 |              2.74 |                      13.84 |                                59.44 |                       13.17 |
| 25%  |               3.99 |                        7.63 |                               95.87 |                      58.21 |         3.44 |                 24.69 |                        63.68 |                8.7  |              5.39 |                      13.84 |                                59.44 |                       26.21 |
| 50%  |               4.27 |                        7.63 |                               95.87 |                      63.64 |        10.75 |                 24.69 |                        63.68 |               30.7  |              6.01 |                      13.84 |                                59.44 |                       29.09 |
| 75%  |               4.55 |                        7.63 |                               95.87 |                      69.01 |        12.75 |                 24.69 |                        63.68 |               36.4  |              6.53 |                      13.84 |                                59.44 |                       30.83 |
| max  |               5.06 |                        7.63 |                               95.87 |                      76.26 |        14.84 |                 24.69 |                        63.68 |               45.48 |              7.22 |                      13.84 |                                59.44 |                       34.37 |

### Using a minimum of 2 matched hashes per unigene (for it to be kept)
|      |   Coelho_dogs_2018 |   Coelho_dogs_2018-combined |   Coelho_dogs_2018-combined:dog-gut |   Coelho_dogs_2018:dog-gut |   TaraOceans |   TaraOceans-combined |   TaraOceans-combined:marine |   TaraOceans:marine |   Zeller_CRC_2014 |   Zeller_CRC_2014-combined |   Zeller_CRC_2014-combined:human-gut |   Zeller_CRC_2014:human-gut |
|:-----|-------------------:|----------------------------:|------------------------------------:|---------------------------:|-------------:|----------------------:|-----------------------------:|--------------------:|------------------:|---------------------------:|-------------------------------------:|----------------------------:|
| mean |               3.39 |                        6.05 |                               95.45 |                      61.47 |         7.28 |                 21.09 |                        58.97 |               22.54 |              5.03 |                      12.18 |                                56.77 |                       25.85 |
| std  |               0.34 |                      nan    |                              nan    |                       7.58 |         3.97 |                nan    |                       nan    |               13.15 |              0.92 |                     nan    |                               nan    |                        4.63 |
| min  |               2.94 |                        6.05 |                               95.45 |                      48.93 |         1.84 |                 21.09 |                        58.97 |                5.12 |              2.2  |                      12.18 |                                56.77 |                       11.3  |
| 25%  |               3.12 |                        6.05 |                               95.45 |                      55.77 |         2.32 |                 21.09 |                        58.97 |                6.24 |              4.6  |                      12.18 |                                56.77 |                       23.84 |
| 50%  |               3.35 |                        6.05 |                               95.45 |                      61.45 |         8.52 |                 21.09 |                        58.97 |               26.2  |              5.16 |                      12.18 |                                56.77 |                       26.63 |
| 75%  |               3.55 |                        6.05 |                               95.45 |                      66.82 |        10.45 |                 21.09 |                        58.97 |               31.73 |              5.59 |                      12.18 |                                56.77 |                       28.34 |
| max  |               4.01 |                        6.05 |                               95.45 |                      74.44 |        12.55 |                 21.09 |                        58.97 |               40.56 |              6.23 |                      12.18 |                                56.77 |                       31.82 |


Note that the _fraction_ of kept unigenes is counted, which will not correspond
to the fraction of basepairs as longer genes are more likely to be kept.

