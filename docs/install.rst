Installation
============

Manual installation
-------------------

nanoRMS2 is written in Python3 and should work in most UNIX systems.

Make sure you install all programs listed below, before runnning the pipeline.
It's recommended to create separate conda environment and use Python version 3.7.

* first the repository

.. code-block:: bash
		
   mkdir ~/src
   cd src
   git clone https://github.com/novoalab/nanoRMS2
   cd nanoRMS2

* from `conda <https://bioconda.github.io/user/install.html#install-conda>`_ (use miniconda3!)

.. code-block:: bash
		
   conda create -y -n "nanoRMS2" python=3.7
   conda activate nanoRMS2
   conda install -y cython hdf5
   conda install -y -c bioconda meme=5.3.0
   
* guppy_basecaller has to be obtained from `Nanopore Tech. Software page <https://community.nanoporetech.com/downloads>`_
  Alternatively, you can try `this for GPU <https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_3.6.1_linux64.tar.gz>`_
  or `this for CPU <https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_3.6.1_linux64.tar.gz>`_ version.
  For GPU basecalling to work, you'll need to install CUDA with NVIDIA drivers.
  Check `my blog for instructions for Ubuntu 18.04 <https://medium.com/@lpryszcz/containers-with-cuda-support-5467f393649f>`_
  or `NVIDIA CUDA website for other systems <https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html>`_.

* from `pip <https://pypi.org/project/pip/>`_

  * pyguppyclient (this will work with guppy v3.6.1,
    for other versions check `here <https://modphred.readthedocs.io/en/latest/install.html#which-pyguppyclient-version-should-i-install>`_)
    and other dependencies

.. code-block:: bash

   pip install "pip<20" "numpy<1.20"
   pip install -r requirements.txt

It's important to install numpy before v1.20 due to `ValueError: cannot convert float NaN to integer <https://github.com/nanoporetech/tombo/issues/319>`_.  
Note, if numpy is updated later on, tombo will raise ``ValueError: numpy.ndarray size changed, may indicate binary incompatibility. Expected 88 from C header, got 80 from PyObject``.  
You may see similar error if `your Python is newer than 3.7 <https://github.com/nanoporetech/tombo/issues/319#issuecomment-785180546>`_. 


Once you have all dependencies installed,
we recommend to try running it with :doc:`test dataset <test>`.
It'll be much easier to troubleshoot all potential issues this way. 


Docker
------
NOTE: the images are not public yet - you'll need to build your own as follows:

.. code-block:: bash

   cd ~/src/nanoRMS2
   docker build --pull -t lpryszcz/nanorms2:latest docker


We maintain `docker image <https://hub.docker.com/repository/docker/lpryszcz/nanorms2>`_
for guppy v3.6.1 (with pyguppyclient v0.0.6).

If you want to use it, make sure you have Docker, GPU drivers, CUDA
and nvidia-docker installed.
The easiest may be to follow `nvidia-docker installation tutorial
<https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html#docker>`_.


In order to execute :doc:`test example <test>`, all you need to do
is to adjust the version of guppy in the below command:

.. code-block:: bash

   cd test
   acc=oligopl; find $acc -name "*.bam"|xargs rm
   docker run --gpus all -u $UID:$GID -v `pwd`:/data lpryszcz/nanorms2 \
     /opt/app/run -f /data/ref/ECOLI.fa -o /data/nanoRMS2/$acc \
     -i /data/$acc/C600_control_WGA/ /data/$acc/C600_native/ \
     -c dna_r9.4.1_450bps_pcr_fast.cfg --host /usr/bin/guppy_basecall_server \
     -b "NC_000913.3:1-100000"
     
   
As you can see, the above command got a bit complicated. This is because:

- we need to enable GPU
- define user & group (otherwise all output files will be owned by root)
- bind local directory within container
- and define all input folders (because autocompletion doesn't work inside the container)

  

