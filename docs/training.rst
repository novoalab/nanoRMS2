Training of guppy models
========================


Note, you'll need to get and install customised
`taiyaki <https://github.com/novoalab/taiyaki>`_
and/or `bonito <https://github.com/novoalab/bonito>`_. 


Data preparation
----------------

Typically, ~16,000 reads should be more than enough to train decent model.
At least for taiyaki. For bonito, having more reads is not a luxury.
Especially for ``sup`` models.

Balancing RNA dataset
^^^^^^^^^^^^^^^^^^^^^

Remember, for RNA, you'll need to balanced your dataset,
because highly expressed transcripts will dominate your data.

.. code-block:: bash

    s=Mettl3-KO; i=1 # WT
    ref=ref/HUMAN.transcripts.fa;
    # you can use BAM or FastQ files
    f=~/cluster/rna_mods/human/modPhred/xPore/minimap2/HEK293T-${s}-rep$i.bam
    # just make sure you align on transcripts
    minimap2 -ax map-ont -k13 $ref <(samtools fastq $f) 2>/dev/null|samtools sort --write-index -o $f.transcripts.bam
    # and sample up to 5 reads per each transcript
    ~/src/taiyaki/bin/bam2sample.py -f $ref -i $f.transcripts.bam > $f.transcripts.bam.ids; wc -l $f.transcripts.bam.ids
    # finall subset your Fast5 files
    fast5_subset -t6 -i /no_backup_isis/enovoa/data/ont/xPore/HEK293T-${s}-rep$i \
      -s reads.balanced/HEK293T-${s}-rep$i -l $f.transcripts.bam.ids



Modification-unaware models
---------------------------

taiyaki
^^^^^^^

We'll train a modification-unaware model, but we'll prepare it such,
that we would be able to train it further if we wish to add
up to 1-modification for each canonical base. 


DNA model
.........

.. code-block:: bash

    cd ~/cluster/rna_mods/curlcake/taiyaki
    source ~/src/taiyaki/venv/bin/activate
    size=256; stride=2; winlen=19;
    opts="--device cuda:0 --step 25 --size $size --stride $stride --winlen $winlen \
      --save_every 500 -r ../../ref/yeast.fa --min_seq_len 50 --max_seq_len 200 \
      --chunk_len 1000 --batch_size 100 --max_chunks 10000000" #DNA 10/base RNA 40/base
    s=DNA_yeast;
    ~/src/taiyaki/bin/train_flipflop_tombo.py $opts --mod_factor 0.01 \
      --alphabet ACGT --mod 1 A 3mA --mod 2 C 5mC --mod 3 G 7mG --mod 4 T 4sU
      --outdir models.tombo/$s.$size.$winlen.$stride \
      -i ../../../dna_mods/yeast/_archives/raw/WO_S1_1_FAL12282/treated_untreated/batch_1?.fast5 \
      -m ~/src/taiyaki/models/mGru_cat_mod_flipflop.py \
      -c dna_r9.4.1_450bps_pcr_fast.cfg \
      --host ~/src/ont-guppy_3.6.1/bin/guppy_basecall_server #--overwrite --max_reads 10



RNA model
.........


.. code-block:: bash

    cd ~/cluster/rna_mods/curlcake/taiyaki
    source ~/src/taiyaki/venv/bin/activate
    size=256; stride=10; winlen=31;
    opts="--rna --device cuda:0 --step 10 --size $size --stride $stride --winlen $winlen \
      --save_every 500 -r ../../ref/curlcake_human_transcripts.fa --chunk_len 4000 \
      --max_per_ref 100 --max_chunks 5000000 --batch_size 50" #DNA 10/base RNA 40/base
    s=IVT_Hs_mRNA_fast5;
    ~/src/taiyaki/bin/train_flipflop_tombo.py $opts --mod_factor 0.01 \
      --alphabet ACGU --mod 1 A m6A --mod 2 C m5C --mod 3 G m5G --mod 4 U 4sU \
      -i reads/$s/*.fast5 --outdir models.tombo/$s.$size.$winlen.$stride \
      -m ~/src/taiyaki/models/mGru_cat_mod_flipflop.py \
      --host ~/src/ont-guppy_3.6.1/bin/guppy_basecall_server #--overwrite --max_reads 10


bonito
^^^^^^


RNA model
.........

.. code-block:: bash

    cd ~/cluster/bonito
    
    # prepare dataset
    conda activate nanoRMS2
    ~/src/taiyaki/bin/get_bonito_dataset.py --rna -o data/rna_xpore50_balanced_rep123123 \
      -r ref/HUMAN.transcripts.fa -i reads.balanced/HEK293T-*-rep?/*.fast5 \
      --min_seq_len 30 --max_seq_len 270 --step 50 --max_chunks 8000000 \
      -c rna_r9.4.1_70bps_ivt_fast.cfg
      
    # run training
    source ~/src/bonito/venv3/bin/activate
    d=rna_xpore50_balanced_rep123123; m=sup; # or hac
    bonito train --config models/configs/rna_r9.4.1_${m}@v3.3.toml \
      --directory data/$d models/$d.$m # --epochs 10 # --batch 128


Modification-aware models
-------------------------


taiyaki
^^^^^^^

RNA model
.........

Here, it's crucial that you already have modification-unaware model trained (``-m``).
This model has to have all modifications of interest already defined in the model.

On top of that, you'll need to provide :doc:`nanoRMS2 models <output>` using ``--models``.
The script will perform on-the-fly feature retrieval and prediction of modified bases.
This will be subsequently used for training. 

.. code-block:: bash

    cd ~/cluster/rna_mods/curlcake/taiyaki
    source ~/src/taiyaki/venv/bin/activate
    size=256; stride=10; winlen=31;
    opts="--rna --device cuda:0 --size $size --stride $stride --winlen $winlen \
      --save_every 500 -r ../../ref/curlcake_human_transcripts.fa \
      --chunk_len 4000 --step 25 --batch_size 50 --max_per_ref 100"
    ko=HEK293T-Mettl3-KO-rep1; sample=HEK293T-WT-rep1;
    ~/src/taiyaki/bin/train_flipflop_tombo_mods.py $opts --mod_factor 1.0 \
      --alphabet ACGU --mod Y A m6A --motifs AGT,AG,A,C,ACT ACGT,C,A,G,ACGT \
      --outdir models.tombo/${sample}_$ko.$size.$winlen.$stride \
      -m models.tombo/IVT_human.$size.$winlen.$stride/model_final.checkpoint \
      --models ~/cluster/rna_mods/human/modPhred_de_novo/guppy3.6.1.ivt_fast/xPore/Mettl3-KO-rep1/models.lzma \
      -k ~/cluster/rna_mods/human/taiyaki/reads/$ko/*.fast5 \
      -i ~/cluster/rna_mods/human/taiyaki/reads/$sample/*.fast5 \
      --max_chunks 1000000 --max_chunks_per_file 15000 \
      -c rna_r9.4.1_70bps_ivt_fast.cfg \
      --host ~/src/ont-guppy_3.6.1/bin/guppy_basecall_server 



bonito
^^^^^^

Currently, only taiyaki can be used to train modification-aware models.


Post-processing
---------------

So you have trained the model. Now what?

Evaluate model
^^^^^^^^^^^^^^

Below will basecall the reads using custom model,
align the reads and report stats from BAM files. 

.. code-block:: bash

    cd ~/cluster/rna_mods/curlcake/taiyaki
    source ~/src/taiyaki/venv/bin/activate
    ver=3.6.1;
    ref=../../ref/curlcake_human_transcripts.fa
    s=HEK293T-WT-rep1_HEK293T-Mettl3-KO-rep1.256.31.10; i=00046; mode=hac;
    mkdir -p guppy$ver/$s.$i;
    # get model json
    m=models.tombo/$s/model_checkpoint_$i.checkpoint;
    if [ ! -s $m.json ]; then dump_json.py $m > $m.json; fi;
    # basecall using new model
    for d in ~/cluster/rna_mods/human/taiyaki/reads/*/; do
        n=`echo $d | rev | cut -f2 -d'/' | rev`;
	if [ ! -d guppy$ver/$s.$i/$n ]; then
	    echo `date` $d $n;
	    ~/src/ont-guppy_$ver/bin/guppy_basecaller --device cuda:0 \
	      -c rna_r9.4.1_70bps_$mode.cfg -m $m.json --compress_fastq --disable_pings \
	      -ri $d -s guppy$ver/$s.$i/$n;
	    minimap2 -ax map-ont $ref guppy$ver/$s.$i/$n/*.fastq.gz 2> guppy$ver/$s.$i/$n.bam.log | samtools sort --write-index -o guppy$ver/$s.$i/$n.bam; 
	    bam2stats.py guppy$ver/$s.$i/$n.bam;
	fi;
    done;
    bam2stats.py guppy$ver/$s.*/*.bam 


Calibrate the model
^^^^^^^^^^^^^^^^^^^



.. code-block:: bash

    cd ~/cluster/rna_mods/curlcake/taiyaki
    m=models.tombo/DNAyeast_xG,xA_NNNxNNN.treated_aag.treated_untreated.256.19.2
    ver=3.6.1; i=final; n=treated_untreated
    s=`basename $m`; echo $s
    ~/src/ont-guppy_$ver/bin/guppy_aligner -i guppy$ver/$s.$i/$n -a $ref -t4 \
      -s guppy$ver/$s.$i/$n.sam --disable_pings
    ~/src/taiyaki/misc/calibrate_qscores_byread.py \
      --alignment_summary guppy$ver/$s.$i/$n.sam/alignment_summary.txt \
      --fastq <(zcat guppy$ver/$s.$i/$n/*.fastq.gz)


Prepare guppy configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ideally, you'd start from exising ``fast``, ``hac`` or ``sup`` model
present in ``~/src/ont-guppy_$ver/data`` and just edit the relevant fields:

- remove (or comment) compatible flowcells - all lines starting with ``compatible_``
- copy your model file (``.json``) to /data and rename it
  ie ``template_rna_r9.4.1_70bps_m6A_hac.jsn``
- update ``qscore_offset`` and ``qscore_scale`` with the values obtained from
  `model calibration <#Calibrate the model>`_.


