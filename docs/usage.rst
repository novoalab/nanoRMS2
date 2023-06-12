Running the pipeline
====================

NanoRMS2 can process DNA and RNA datasets.
The only difference between the two is enabling splice-aware alignments for RNA
and `the basecalling model used <https://github.com/novoalab/nanoRMS2/tree/main/data>`_:
PCR for DNA and IVT for RNA. 
By default nanoRMS2 assumes your dataset is DNA and PCR guppy model is used
``-c dna_r9.4.1_450bps_pcr_hac.cfg``.
If you wish to analyse RNA dataset, all you need to do is to
specify IVT guppy model by ``--rna -c rna_r9.4.1_70bps_ivt_hac.cfg``

The only required inputs are:

* reference FastA sequence
* two paths containing Fast5 files: reads without (PCR/IVT/KO) and with modifications (wt)

nanoRMS can be executed in three modes:

* local on-the-fly basecalling
* remote on-the-fly basecalling
* without basecalling (assuming the Fast5 files were basecalled before)

We recommend to use on-the-fly basecalling because
it's a few times faster than running basecalling and nanoRMS2 separately.

In order to see detailed description of program parameters,
just execute it with ``-h / --help``.


Local on-the-fly basecalling
----------------------------
Here, we assume, that you have guppy already installed in you system.
nanoRMS2 will start guppy_basecall_server in the background
and stop it when it isn't needed anymore.
All you need to do is to provide path to ``guppy_basecall_server`` using ``--host``

.. code-block:: bash

   ~/src/nanoRMS2/run -e -o nanoRMS2/projectName -f reference.fasta -i fast5_folder_KO fast5_folder_WT \
    --host ~/src/ont-guppy_3.6.1/bin/guppy_basecall_server


Alternatively, if guppy_basecall_server is already running in your machine,
you can provide just its port using --host localhost --port.

Remote on-the-fly basecalling
-----------------------------
Here, we assume the guppy_basecall_server is already running in the remote machine.
All you need to do is to provide IP address ``--host`` and port ``--port``

.. code-block:: bash

   ~/src/nanoRMS2/run -e -f reference.fasta \
    -o nanoRMS2/projectName -i fast5_folder_KO fast5_folder_WT \
    --host 172.21.11.186 --port 5555


Without basecalling
-------------------
Make sure your Fast5 files are basecalled with `PCR or IVT model
distributed with nanoRMS2 <https://github.com/novoalab/nanoRMS2/tree/main/data>`_. 
Then running the pipeline is as easy as:

.. code-block:: bash

   ~/src/nanoRMS2/run -e -f reference.fasta \
    -o nanoRMS2/projectName -i fast5_folder_KO fast5_folder_WT


For more usage examples, please have a look in :doc:`test dataset <test>`.


Processing (very) large datasets
--------------------------------
There are several ways of speeding up entire analysis for large datasets:

* get_features: process each sample separately. This step benefits from GPU basecalling,
  but the slowest part (by far) is resquiggle that is CPU-bound, so
  provide as many threads as possible using ``-t``
* predict_mods: this step can be both, CPU and memory hungry.
  If you expect the modifications to be distributed
  (more or less) randomly in all chromosomes, you could analyse just 1 chromosome
  or even the representative region like we do in :doc:`test examples <test>`
  to obtain decent modification model.  
* endode_mods: provide as many threads as possible using ``-t``. 

