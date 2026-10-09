# EDA

Train: 2,539,608 spectra, 275,810 structures. Test sample: 1,213 spectra, 400 molecules.

## Libraries

| ingest_lib         |    spectra |   structures |   med_peaks |   pos_frac |   med_ppm |
|:-------------------|-----------:|-------------:|------------:|-----------:|----------:|
| enveda-180         | 1153785.00 |    182941.00 |      138.00 |       0.77 |      0.98 |
| pluskal_ms2        |  527581.00 |     46821.00 |       21.00 |       0.61 |      1.66 |
| riken              |  347171.00 |     15892.00 |        9.00 |       0.90 |      3.75 |
| gnps               |  220849.00 |     45750.00 |       46.00 |       0.81 |      1.67 |
| massbank           |  101727.00 |      9180.00 |       11.00 |       0.67 |      1.95 |
| mona               |   92416.00 |     11681.00 |       57.00 |       0.75 |      1.80 |
| spectraverse       |   50933.00 |      9631.00 |       47.00 |       0.66 |      0.65 |
| msdial             |   40765.00 |      9127.00 |       11.00 |       0.58 |      1.47 |
| drug_plus          |    2545.00 |      2539.00 |      108.00 |       1.00 |      1.73 |
| enveda-np-examples |    1184.00 |       250.00 |      141.50 |       0.67 |      1.65 |
| masaryk            |     652.00 |       416.00 |        9.00 |       1.00 |      1.48 |


## Top train adducts

| adduct         |   train | in_test_set   |
|:---------------|--------:|:--------------|
| [M+H]+         | 1324585 | True          |
| [M-H]-         |  465688 | True          |
| [2M+Na]+       |  190033 | False         |
| [M+Na]+        |  178092 | True          |
| [2M+H]+        |   92041 | False         |
| [M+CH2O2-H]-   |   83188 | True          |
| [2M-H]-        |   56172 | False         |
| [M+NH4]+       |   53452 | True          |
| [M+Cl]-        |   27511 | True          |
| [M-H2O+H]+     |   21682 | False         |
| [M+K]+         |   18955 | True          |
| [M]+           |    8818 | False         |
| [M+C2H4O2-H]-  |    7436 | False         |
| [M-2H2O+H]+    |    3439 | False         |
| [M+2H]2+       |    1717 | False         |
| [2M+Na-2H]-    |    1431 | False         |
| [2M+CH2O2-H]-  |     665 | False         |
| [M-H2O]+       |     624 | False         |
| [M-CH3]-       |     619 | False         |
| [M-2H]-        |     454 | False         |
| [M]-           |     446 | False         |
| [2M+NH4]+      |     378 | False         |
| [M-H2O-H]-     |     285 | False         |
| [M-H5O3]+      |     230 | False         |
| [2M+C2H4O2-H]- |     187 | False         |


Share of train spectra with a test-set adduct: 84.7%

## Test adducts / mode

| adduct       |   spectra |
|:-------------|----------:|
| [M+H]+       |       959 |
| [M-H]-       |       193 |
| [M+CH2O2-H]- |        31 |
| [M+Na]+      |        22 |
| [M+NH4]+     |         4 |
| [M+K]+       |         2 |
| [M+Cl]-      |         2 |


| ionization_mode   |   spectra |
|:------------------|----------:|
| positive          |       987 |
| negative          |       226 |


## Spectra per structure (train)

|       |    n_spectra |
|:------|-------------:|
| count | 275810       |
| mean  |      9.20782 |
| std   |     21.0696  |
| min   |      1       |
| 25%   |      4       |
| 50%   |      6       |
| 75%   |      9       |
| max   |   1468       |


## Spectra per molecule (test sample)

|       |   n_spectra |
|:------|------------:|
| count |   400       |
| mean  |     3.0325  |
| std   |     1.38247 |
| min   |     1       |
| 25%   |     2       |
| 50%   |     3       |
| 75%   |     4       |
| max   |     9       |


## Collision-energy list length

|   n_ce |   train |   test |
|-------:|--------:|-------:|
|      0 |  337443 |      0 |
|      1 | 1706663 |    915 |
|      2 |   12562 |      0 |
|      3 |  394599 |    298 |
|      4 |   83712 |      0 |
|      5 |    4613 |      0 |
|      6 |      16 |      0 |


## Precursor error (ppm), train

|       |      abs_ppm |
|:------|-------------:|
| count |   2535322.00 |
| mean  |      1845.73 |
| std   |    511752.51 |
| min   |         0.00 |
| 50%   |         1.51 |
| 90%   |         5.04 |
| 99%   |      2884.42 |
| max   | 469771545.80 |


## Precursor m/z, train vs test

|       |     train |   test |
|:------|----------:|-------:|
| count | 2539608.0 | 1213.0 |
| mean  |     419.8 |  327.9 |
| std   |    1536.0 |   38.1 |
| min   |       2.0 |  245.1 |
| 25%   |     308.1 |  300.1 |
| 50%   |     350.1 |  329.2 |
| 75%   |     489.2 |  347.2 |
| max   | 2430333.0 |  460.2 |


## Instrument-matched rows (enveda-*)

1,154,969 spectra, 183,191 structures.
