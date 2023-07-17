# NanoRMS2

NanoRMS2 is a pipeline for detection of DNA/RNA modifications.
It requires two samples that differ in modification status, for example KO and WT.

The pipeline consists of three steps / modules:
- get_features: retrieve & store features in BAM file
- predict_mods: train ML models for modified bases, estimate per-position stoichiometry
  and much more...
- encode_mods: predict & encode per-read modification status in BAM file
  using ML models generated by predict step. 

All these steps can be run separately or as a pipeline by executing `nanoRMS2/run`.  

## Requirements
You can find the software versions and requirements [here](./requirements.txt).

## Additional documentation
For more information, please visit 
[full documentation](https://public-docs.crg.es/enovoa/public/lpryszcz/src/nanoRMS2/readthedocs).  

## Base-callin models 
The basecalling models are stored in [/data](/data) directory. 

## Citation

If you find this work useful, please cite:

Sonia Cruciani \*, 
Anna Delgado-Tejedor \*,
Leszek P. Pryszcz <sup>#,</sup>\*,
Rebeca Medina, Laia Llovera and Eva Maria Novoa <sup>#</sup>
(2023)
*De novo* basecalling of m6A RNA modifications at single-molecule and single-nucleotide resolution using direct RNA nanopore sequencing
In preparation.

Note, the order of co-first authorship (\* above) was determined
by flipping the coin on 15 Jun 2023 at 14:53 CET
[@41.3848953,2.1939749](https://www.google.com/maps?q=loc:41.3848953,2.1939749).  
In our opinion, such process summarises pretty well
the randomness of success in science (and life in general). 
