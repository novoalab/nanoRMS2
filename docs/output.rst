Program output
==============
If the program was executed with command such as 

.. code-block:: bash

   ~/src/nanoRMS2/run -e -f reference.fasta \
    -o nanoRMS2/projectName -i fast5_folder_KO fast5_folder_WT


it'll store the features in BAM files named as 
``fast5_folder_KO.bam`` and ``fast5_folder_WT.bam``. 
Additionally, it'll generate the output directory ``nanoRMS2/projectName`` containing:

- KS.tsv.gz - statistic and P-value for two features (SI_0 and TR_0) for all 7-mers
  reported by `two-sample Kolmogorov-Smirnov test <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html>`_
- denovo.tsv.gz - internal format with all predicted positions with likely modifications
- models.lzma - joblib dump of trained models classifying reads as un- or modified
- feature_importances.tsv.gz - importance of features for every model (7-mer)
- meme/meme.html - motif enrichment analysis report
- denovo.motifs.tsv.gz - similar to denovo.tsv.gz with added motif information
- denovo.bed & denovo.motifs.bed - bedMethyl formatted files representing modified positions
  reported in denovo.tsv.gz & denovo.motifs.bed, respectively
- encode/\*.bam - BAM files with encoded modification probabilities

  
All \*.tsv.gz files are pandas DataFrames saved in tab-delimited format.


Data formats
------------
While using nanoRMS2 it'll be good to familiarise yourself with couple of data formats.

bedMethyl
^^^^^^^^^
modPhred reports information about modified positions in bedMethyl format.
This is tab-delimited files, compatible with BED (positions are 0-based, half-open)
with several additional fields:

#. Reference chromosome or scaffold
#. Start position in chromosome
#. End position in chromosome
#. Name of item - short modification name
#. Score - modification frequency difference between two samples.
#. Strandedness, plus (+), minus (-), or unknown (.)
#. Start of where display should be thick (start codon) - same as 2.
#. End of where display should be thick (stop codon) - same as 3.
#. Color value (RGB) - different color for modifications each canonical base. 
   Color intensity depends on modification frequency
   (stronger means modification is more frequent).
#. Coverage - number of reads at this position
#. Percentage of modified reads - similar to field #5, but expressed as rounded percentage

For example, the output for ``NC_000913.3:1,061-1,253`` looks like this:

.. code-block:: bash

   NC_000913.3     1091    1092    RCCWGGY 0.667   -       1091    1092    170,0,0 19      67
   NC_000913.3     1166    1167    NGATCNN 0.2     +       1166    1167    51,0,0  10      20
   NC_000913.3     1167    1168    NGATCNN 0.11    -       1167    1168    28,0,0  19      11
   NC_000913.3     1169    1170    NGATCNN 0.333   -       1169    1170    85,0,0  19      33
   NC_000913.3     1207    1208    RCCWGGY 0.2     +       1207    1208    51,0,0  10      20
   NC_000913.3     1208    1209    RCCWGGY 0.2     +       1208    1209    51,0,0  10      20
   NC_000913.3     1208    1209    RCCWGGY 0.875   -       1208    1209    223,0,0 17      88
   NC_000913.3     1209    1210    RCCWGGY 0.2     +       1209    1210    51,0,0  10      20

bedMethyl files can be visualised in many genome browsers ie IGV.


BAM
^^^
A binary format for storing raw genomic data. Reads in BAM files are
typically aligned to reference, compressed and sorted by reference position.

nanoRMS2 reports two types of BAM files:

- ``get_features`` store read features as SAM tags
- ``encode_mods`` store base modification probability for every base
  instead of base qualities. 


Tags in get_features BAM
~~~~~~~~~~~~~~~~~~~~~~~~
We store features using several custom BAM tags:

- per base normalised signal intensity mean [tag si:B:f]
- reference base probability [tag tr:B:C] and probabilities of A [tag tA:B:C], C [tag tC:B:C], 
  G [tag tG:B:C] and T [tag tT:B:C] retrieved from guppy (trace scaled 0-255)
- dwell time [tag dt:B:C] in signal step capped at 255
- probability of modification calculated by tombo [tag mp:B:f]

All features are matched versus padded reference sequnce blocks 
excluding introns and large deletions from reference. 
Those blocks (2-D array of start & ends) are stored as flattened 1-D array [tag bs:B:i]
ie. exons [(8114, 8244), (8645, 8797)] will be stored as array('I', [8114, 8244, 8645, 8797]). 

More information regarding features can be found in
Supplementary Methods in nanoRMS2 `paper <index.html#paper>`__.


