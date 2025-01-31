# Double-Sided Greedy Median Algorithm for Fast and Accurate Trace Reconstruction in DNA Data Storage
Double-Sided Greedy Median (DSGM) Algorithm is a Trace Reconstruction algorithm that outputs reconstructed DNA sequences from clusters of noisy strands.

In this repository, we present the implementation of the DSGM proposed [[1]](#1). 

Python 3 is required (we recommend 3.9 or newer)


## Table of Contents
- [Installation](#installation)
- [Usage](#usage)


## Installation
To try out DSGM, clone this repository to a directory of your choice with the command:

```shell
$ git clone https://github.com/prongs1996/DSGMRecon.git
```
Then install the Levenshtein module
```shell
$ cd DSGMRecon/Levenshtein

$ pip install .

$ cd ..
```

## Usage

The following three algorithms are implemented:

- **Double-Sided Greedy Median (DSGM)**
- **Double-Sided Greedy Median - Beam (DSGM-Beam)**
- **Double-Sided Greedy Median - Refined (DSGM-Refined)**

To generate the reconstructed strand, you will need to pass the cluster of sequences as a list of strings to the respective function. we have provided a demo script ([**demo.py**](./demo.py)) that uses our nanopore sequencing dataset found in */data/our_nanopore_UnderlyingClusters.txt*. 

A reference file with the ground truth strands can also be specified to compute the accuracy and average edit distance of the reconstructed strands. The name of the output file which contains all the reconstructed strands can also be specified. Otherwise, the default output file will be *ReconstructedStrands.txt* in the */data* folder.


To run the demo on DSGM or DSGM-Refined, use the command:
```shell
$ python3 demo.py --i $input_file --ALG $algorithm_number [--o $output_file] [--r $reference_file]
```
For DSGM-Beam, an additional parameter *beamsize* needs to be specified. Use the command:
```shell
$ python3 demo.py --i $input_file --ALG $algorithm_number [--o $output_file] [--r $reference_file] --b $beamsize
```
**Algortihm Number and Implemented Algorithms**: <br>
        0 - *DSGM* <br>
        1 - *DSGM-Beam* <br>
        2 - *DSGM-Refined* <br>
        

For comparisons with prior work, use [DNA Storage Toolkit](https://github.com/prongs1996/DNAStorageToolkit) [[1]](#1) for BMA, DBMA and NWA trace reconstruction. Use the [Reconstruction Algorithms for DNA Storage Systems
](https://github.com/omersabary/Reconstruction/) repo [[2]](#2) for Iterative and DivBMA trace reconstruction. Use the [TrellisBMA](https://github.com/microsoft/TrellisBMA) repo[[3]](#3) for Trellis BMA reconstruction.

## References
<a id="1">[1] </a> Sharma, P., Yipeng, G. G., Gao, B., Ou, L., Lin, D., Sharma, D., & Jevdjic, D. (2024, May). **DNA Storage Toolkit: A Modular End-to-End DNA Data Storage Codec and Simulator**. In 2024 IEEE International Symposium on Performance Analysis of Systems and Software (ISPASS) (pp. 144-155). IEEE.

<a id="2">[2] </a> Sabary, O., Yucovich, A., Shapira, G., & Yaakobi, E. (2024). **Reconstruction algorithms for DNA-storage systems**. Scientific Reports, 14(1), 1951.

<a id="3">[3] </a> Srinivasavaradhan, S. R., Gopi, S., Pfister, H. D., & Yekhanin, S. (2021, July). **Trellis BMA: Coded trace reconstruction on IDS channels for DNA storage**. In 2021 IEEE International Symposium on Information Theory (ISIT) (pp. 2453-2458). IEEE.

