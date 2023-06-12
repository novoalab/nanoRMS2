# NanoRMS2

NanoRMS2 is a pipeline for detection of DNA/RNA modifications.
It requires two samples that differ in modification status, for example KO and WT.

The pipeline consists of three steps / modules:
- get_features: retrieve & store features in BAM file
- predict_mods: train ML models for modified bases, estimate per-position stoichiometry and [much more](LINK)
- encode_mods: predict & encode per-read modification status in BAM file using ML models generated by predict step. 

All these steps can be run separately or as a pipeline by executing `nanoRMS2/run`.  

For more information, please visit 
[full documentation](https://public-docs.crg.es/enovoa/public/lpryszcz/src/nanoRMS2/readthedocs).  

The basecalling models are stored in [/data](/data) directory. 

## Citation

If you find this work useful, please cite:
Pryszcz L, etc... bioRxiv


