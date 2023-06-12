Test datasets
=============

Below you'll find detailed information on running nanoRMS2 pipeline on test dataset.

Download test data
------------------
.. code-block:: bash

    mkdir -p test
    cd test
    wget https://public-docs.crg.es/enovoa/public/lpryszcz/src/nanoRMS2/test/ -q --show-progress -r -c -nc -np -nH --cut-dirs=6 --reject="index.html*"


Run nanoRMS2
------------
First, you'll need to copy the modification-unware model files to guppy folder,
if you haven't done so yet.

.. code-block:: bash
		
    rsync -av ~/src/nanoRMS2/data ~/src/ont-guppy_3.6.1


DNA
^^^
Running entire pipeline with on-the-fly basecalling:

.. code-block:: bash

    acc=oligopl
    ~/src/nanoRMS2/run -e -f ref/ECOLI.fa -i $acc/C600_{control_WGA,native}/ -o nanoRMS2/C600_native \
      -c dna_r9.4.1_450bps_pcr_hac.cfg --host ~/src/ont-guppy_3.6.1/bin/guppy_basecall_server \
      -b "NC_000913.3:1-100000"


The pipeline can be also execute with Fast5 files
that were previously basecalled with modification-unaware (PCR/IVT) model. 
In such case, you don't need to specify the ``--host`` and ``-c``.

If you don't have GPU, you can use CPU by specifying ``--device ""``,
but get_features step will be veeery slow (due to basecalling on CPU). 

Instead you can run all steps one-by-one as follow:

.. code-block:: bash

    acc=oligopl
    ~/src/nanoRMS2/src/get_features.py -f ref/ECOLI.fa -i $acc/C600_control_WGA/ $acc/C600_native/ -c dna_r9.4.1_450bps_pcr_hac.cfg --host ~/src/ont-guppy_3.6.1/bin/guppy_basecall_server
    ~/src/nanoRMS2/src/predict_mods.py -f ref/ECOLI.fa -i $acc/C600_control_WGA/ $acc/C600_native/ -o nanoRMS2/C600_native -b "NC_000913.3:1-100000"
    ~/src/nanoRMS2/src/encode_mods.py -f ref/ECOLI.fa -i $acc/C600_control_WGA/ $acc/C600_native/ -o nanoRMS2/C600_native/encode -m nanoRMS2/C600_native/models.lzma


RNA
^^^
If working with RNA, make sure to specify ``--rna``
and IVT model ``-c rna_r9.4.1_70bps_ivt_hac.cfg``.

.. code-block:: bash

    ~/src/nanoRMS2/run --rna -e -f ref/Yeast_sk1.fa \
     -o nanoRMS2/yeast_ime4/ -i yeast_ime4/{ko,wt}.1/ \
     -c rna_r9.4.1_70bps_ivt_hac.cfg --host ~/src/ont-guppy_3.6.1/bin/guppy_basecall_server

    
Reuse the ML modification model
-------------------------------
Now, since you have modification model, you could also detect modifications in new sample(s) as follows:

.. code-block:: bash

    ~/src/nanoRMS2/src/encode_mods.py -f ref/ECOLI.fa -m nanoRMS2/C600_native/models.lzma \
     -i $acc/C600_native_flongle/ -o nanoRMS2/C600_native/encode \
     -c dna_r9.4.1_450bps_pcr_hac.cfg --host ~/src/ont-guppy_3.6.1/bin/guppy_basecall_server


Above does ``get_feature`` step and classification on-the-fly
without storing BAM with features.
If you wish to store BAM with features, just execute ``get_feature`` before
(and encode_mods will use BAM with features automatically):

.. code-block:: bash

    ~/src/nanoRMS2/src/get_features.py -f ref/ECOLI.fa -i $acc/C600_native_flongle/ \
     -c dna_r9.4.1_450bps_pcr_hac.cfg --host ~/src/ont-guppy_3.6.1/bin/guppy_basecall_server
     
    ~/src/nanoRMS2/src/encode_mods.py -f ref/ECOLI.fa -m nanoRMS2/C600_native/models.lzma \
     -i $acc/C600_native_flongle/ -o nanoRMS2/C600_native/encode 


Note, each modification-model that you train with nanoRMS2 is sample/condition-specific.
If in a given modification has different sequence motifs across samples/conditions,
the model will be limited to particular sample/condition it was trained on.
For example, E. coli DNA model that we trained above will detect 5mC/6mA
only in E. coli dcm/dam sequence context (CCWGG and GATC, respectively).
It won't work for 5mC/6mA in other microorganisms that may have different motifs. 
However, RNA m6A model (GGACA or GGACT) we have trained on yeast Ime4
will likely work decently with vertebrates, since the m6A motif is quite conserved.


Visualisation
-------------
First, let's look at nanoRMS2 predictions. 
Since modificiation probabilites are encoded as base qualities,
you can :doc:`visualise them directly in the genome browser <igv>`.

.. image:: NC_000913.3:22,653-22,809.png
   :align: center


Compare results to guppy
------------------------
Now, let's looks at the modifications basecalled by guppy models.
To do so, you'll need to run
`modPhred <https://github.com/novoalab/modPhred>`_ as described below
and then :doc:`load resulting BAM and BED files into IGV <igv>`. 


Old 6mA/5mC model
^^^^^^^^^^^^^^^^^
Let's looks at the modifications predicted from
dam/dcm (GATC/CCWGG sequence contexts) model
that is shipped with guppy version < 4.5.2:

.. code-block:: bash

    pip install "pyguppyclient==0.0.6"
    ~/src/modPhred/run -f ref/ECOLI.fa -o modPhred/$acc \
     -i $acc/*/ -t4 --host ~/src/ont-guppy_3.6.1/bin/guppy_basecall_server

.. image:: NC_000913.3:22,653-22,809.modPhred_damdcm.png
   :align: center

As you can see, the results are very similar between nanoRMS2 and guppy dam/dcm model,
but dam/dcm guppy model:

- seems to falsely predict many 6mA bases in PCR sample
- works poorly with flongle


New 5mC model
^^^^^^^^^^^^^
Finally, let's looks at the modifications predicted from
5mC model for all sequence contexts
that is shipped with guppy version >= 4.5.2:

.. code-block:: bash

    pip install "pyguppyclient==0.1.0"
    ~/src/modPhred/run -t4 --timeout 6000 -f ref/ECOLI.fa \
     -o modPhred/$acc.5011 -i $acc/*/ \
     --host ~/src/ont-guppy_5.0.11/bin/guppy_basecall_server -c dna_r9.4.1_450bps_modbases_5mc_hac.cfg

    
.. image:: NC_000913.3:22,653-22,809.modPhred_5mC.png
   :align: center

As you can see, the results are very different from both, old 5mC/6mA and nanoRMS2.
New 5mC guppy model:

- wrongly predicts many 5mC in PCR sample
- wrongly predicts 5mC outside of CCWGG context in native sample
- and doesn't predict 5mC in the CCWGG context in native sample
