# Basecalling models

This folder stores modification-unaware basecalling models
needed to run nanoRMS2 for:
- DNA: `*pcr*`
- RNA: `*ivt*`

Those models can be used directly with guppy_basecaller using `-c` parameter.
For example, in order to basecall using `ivt_hac` RNA model, just execute
(we're assuming the guppy binaries are in `~/src/ont-guppy_3.6.1/bin`):

```bash
# copy (or link) model files into guppy data directory
rsync -av ~/src/nanoRMS2/data ~/src/ont-guppy_3.6.1

# and execute guppy with ivt_hac model
~/src/ont-guppy_3.6.1/bin/guppy_basecaller -c rna_r9.4.1_70bps_ivt_hac.cfg -i fast5_input -s output_dir --fast5_out --compress_fastq -x cuda:0
```

You can find more information about the training procedure in the
[documentation]https://github.com/novoalab/nanoRMS2/docs]. 
