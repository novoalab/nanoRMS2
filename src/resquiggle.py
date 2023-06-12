#!/usr/bin/env python3
"""Resquiggle related functions.

TODO:
- try to rescue reads that didn't resquiggle ie increasing the bandwith
"""

import numpy as np, pysam, scipy
#import statsmodels.api as sm
from tombo import tombo_stats, resquiggle, tombo_helper
from tombo._default_parameters import OUTLIER_THRESH, SHIFT_CHANGE_THRESH, SCALE_CHANGE_THRESH, RNA_SAMP_TYPE, DNA_SAMP_TYPE, COLLAPSE_RNA_STALLS, COLLAPSE_DNA_STALLS, STALL_PARAMS#, FM_OFFSET_DEFAULT
from copy import deepcopy

VERSION = '1.0a'
DEFAULT_STALL_PARAMS = tombo_helper.stallParams(**STALL_PARAMS)
USE_START_CLIP_BASES = resquiggle.USE_START_CLIP_BASES

COLLAPSE_DNA_STALLS = True # False
n_windows = 2
mini_window_size = 40 # 50
MEAN_STALL_PARAMS_DNA = dict((
    ('window_size', n_windows * mini_window_size), ('threshold', 40),
    ('edge_buffer', 25),  ('min_consecutive_obs', 75),
    ('n_windows', n_windows), ('mini_window_size', mini_window_size)))
DEFAULT_STALL_PARAMS_DNA = tombo_helper.stallParams(**MEAN_STALL_PARAMS_DNA)

def get_tombo_vars(rna, sensitive=False):
    """Return tombo vars for DNA or RNA (if rna==True)."""
    if rna:
        seq_samp_type = tombo_helper.seqSampleType('RNA', True)
        rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)
        std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
        if sensitive: rsqgl_params = rsqgl_params._replace(**{"bandwidth": 900})
    else:
        seq_samp_type = tombo_helper.seqSampleType('DNA', False)
        rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)
        std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
    return seq_samp_type, std_ref, rsqgl_params

def adjust_map_res(map_res, seq_samp_type, rsqgl_params, TRIM_RNA_ADAPTER=False):
    if seq_samp_type.name == RNA_SAMP_TYPE:
        if TRIM_RNA_ADAPTER:
            # trim DNA adapter off of RNA signal
            adapter_end = tombo_stats.trim_rna(map_res.raw_signal, rsqgl_params)
            # trim off adapter
            map_res = map_res._replace(raw_signal=map_res.raw_signal[adapter_end:])
        # flip raw signal for re-squiggling
        map_res = map_res._replace(raw_signal=map_res.raw_signal[::-1])
    elif seq_samp_type.name == DNA_SAMP_TYPE and USE_START_CLIP_BASES:
        # flip raw signal, genome and start clip seqs for re-squiggling
        map_res = map_res._replace(
            raw_signal=map_res.raw_signal[::-1],
            genome_seq=map_res.genome_seq[::-1])
        
    if COLLAPSE_RNA_STALLS and seq_samp_type.name == RNA_SAMP_TYPE:
        map_res = map_res._replace(stall_ints=tombo_stats.identify_stalls(map_res.raw_signal, DEFAULT_STALL_PARAMS))
    elif COLLAPSE_DNA_STALLS and seq_samp_type.name == DNA_SAMP_TYPE:
        map_res = map_res._replace(stall_ints=tombo_stats.identify_stalls(map_res.raw_signal, DEFAULT_STALL_PARAMS_DNA))

    return map_res

def adjust_rsqgl_res(rsqgl_res, all_raw_signal, seq_samp_type, USE_START_CLIP_BASES=False):
    if seq_samp_type.name == DNA_SAMP_TYPE and USE_START_CLIP_BASES:
        # flip raw signal and events back for storage in genome direction
        rev_rsrtr = (all_raw_signal.shape[0] -
                     rsqgl_res.read_start_rel_to_raw -
                     rsqgl_res.segs[-1])
        rev_segs = -1 * (rsqgl_res.segs[::-1] - rsqgl_res.segs[-1])
        rsqgl_res = rsqgl_res._replace(
            read_start_rel_to_raw=rev_rsrtr, segs=rev_segs,
            genome_seq=rsqgl_res.genome_seq[::-1],
            raw_signal=rsqgl_res.raw_signal[::-1])

    return rsqgl_res

def get_exonic_blocks(a):
    """Return exonic blocks this is start-end reference-based 
    for consecutive exons covered by given read.
 
    Note, those are not necesarily exact exons, just exons infered from read alignment. 
    """
    blocks = []
    s = e = a.pos
    # iter read blocks
    for code, bases in a.cigar:
        # count blocks that alter reference positions (ignore ie insertions [1])
        if code in (0, 2, 7, 8): e += bases
        # exclude introns - those are reported as reference-padded alignment part
        elif code == 3:
            blocks.append([s, e])
            s = e + bases
            e = s
    # store exon after last intron (or entire transcript if no introns)
    blocks.append([s, e])
    return blocks

def map_read(a, faidx, seq_samp_type, std_ref, ref2len):
    """Get resquiggle result with read alignement info"""
    seq_data = tombo_helper.sequenceData(seq=a.seq, id=a.qname, mean_q_score=np.mean(a.query_qualities))
    # get chrom, start and end
    chrm, ref_start, ref_end = a.reference_name, a.reference_start, a.reference_end
    # store strand & number of clipped bases relative to read sequence
    if a.is_reverse:
        strand = "-"
        num_start_clipped_bases = len(seq_data.seq) - a.qend
        num_end_clipped_bases = a.qstart
    else:
        strand = "+"
        num_start_clipped_bases = a.qstart
        num_end_clipped_bases = len(seq_data.seq) - a.qend
    
    # 'ID', 'Subgroup', 'ClipStart', 'ClipEnd', 'Insertions', 'Deletions', 'Matches', 'Mismatches'
    align_info = tombo_helper.alignInfo(seq_data.id, "", num_start_clipped_bases, num_end_clipped_bases,
                                        0, 0, a.alen, 0) # this isn't used anywhere, so just don't bother computing it!
    # extract genome sequence from mappy aligner
    # expand sequence to get model levels for all sites (need to handle new
    # sequence coordinates downstream)
    start_skip = 0
    # get exonic blocks
    blocks = get_exonic_blocks(a)
    align_info.blocks = deepcopy(blocks)
    # here there is likely issue as DNA mer is 6 (2-nd central), while RNA 5 (1-st central)
    # but RNA is 3>5', thus in fact is should be written as RNA 5 (4th central), right?
    # on top of that squiggle seems to be shifted by -1 base relative to tr (and true mod)
    dnstrm_bases = std_ref.kmer_width - std_ref.central_pos - 1 
    #if ((seq_samp_type.name == RNA_SAMP_TYPE and strand == '-') or # this should be correct
    if ((seq_samp_type.name == RNA_SAMP_TYPE and strand == '+') or
        (seq_samp_type.name == DNA_SAMP_TYPE and strand == '-' and USE_START_CLIP_BASES) or
        (seq_samp_type.name == DNA_SAMP_TYPE and strand == '+' and not USE_START_CLIP_BASES)):
        if ref_start < std_ref.central_pos: # this bitch can cause problems for reads aligning at the ends of reference
            start_skip = std_ref.central_pos-ref_start
            ref_start = std_ref.central_pos
        ref_seq_start = ref_start - std_ref.central_pos
        ref_seq_end = ref_end + dnstrm_bases
    else:
        if ref_start < dnstrm_bases: # this bitch can cause problems for reads aligning at the ends of reference
            start_skip = dnstrm_bases-ref_start
            ref_start = dnstrm_bases
        ref_seq_start = ref_start - dnstrm_bases
        ref_seq_end = ref_end + std_ref.central_pos
    # update blocks start & end with kmer specific shifts - this sequence won't be saved! 
    blocks[0][0] = ref_seq_start
    blocks[-1][1] = ref_seq_end
    # get exonic sequence
    genome_seq = "".join([faidx.fetch(chrm, s, e) for s, e in blocks]) #faidx.fetch(chrm, ref_seq_start, ref_seq_end) #aligner.seq(chrm, ref_seq_start, ref_seq_end)
    # get missing bases in the end
    end_skip = 0 if blocks[-1][1]<=ref2len[chrm] else blocks[-1][1]-ref2len[chrm] #ref_seq_end-ref_seq_start-len(genome_seq)
    # enlarge genome seq by missing bits from ends with (random!) bases - As for now
    if start_skip or end_skip:
        genome_seq = "A"*start_skip + genome_seq + "A"*end_skip
    if strand == '-':
        genome_seq = tombo_helper.rev_comp(genome_seq)
    # store enlarged genome for P-value calculation, so no trimming needed later :)
    genome_seq = genome_seq.upper() #.upper() is important to correctly process soft-masked sequences
    #if "N" in genome_seq: genome_seq = genome_seq.replace("N", "A") # # dirty N>A replacement for tRNA
    align_info.refseq = genome_seq # res.genome_seq is altered during find_adaptive_assignment
    genome_loc = tombo_helper.genomeLocation(ref_start, strand, chrm)
    return tombo_helper.resquiggleResults(align_info, genome_loc, genome_seq, seq_data.mean_q_score)

def resquiggle_reads(aligner, ref, rna=False, add_si=True,  
                     outlier_thresh=OUTLIER_THRESH, max_scaling_iters=3, max_per_ref=0, 
                     valid_bases=set(list('ACGT'))):
    """Resquiggle aliged reads from aligner (basecall.basecall_and_align) 
    and yield tombo.resquiggle object and error information for every alignment. 

    Error is empty string if resquiggle was performed correctly. 
    If the error was encountered, the error information is returned 
    and None is returned instead of tombo.resquiggle object.

    In addition, tombo.resquiggle object is extended with: 
    - a: read alignment as pysam.AlignmentSegment object
    - sig: raw reads signal (without triming and scaling)
    - move: move table from the basecaller
    - trace: trace table from the basecaller
    - si: per-base mean signal intensity (if add_si==True)
    - mp: per-base modification probability (if add_si==True)
    """
    ref2c = {}
    # load tombo model & its parameters
    seq_samp_type, std_ref, rsqgl_params = get_tombo_vars(rna)    
    # process reads from multi fast5
    faidx = pysam.FastaFile(ref)
    ref2len = {r: l for r, l in zip(faidx.references, faidx.lengths)}
    for a, sig, move, trace, md in aligner:
        # process only given number of reads per reference
        if max_per_ref:
            contig = a.reference_name #map_results.genome_loc.Chrom
            if contig in ref2c:
                if ref2c[contig]>=max_per_ref: continue
            else: ref2c[contig] = 0
        # skip reads without alignment or secondary/qcfails
        if a.is_unmapped or a.is_secondary or a.is_qcfail:
            yield None, "No alignment" if a.is_unmapped else "Secondary alignment"
            continue
        # get alignment data
        map_results = map_read(a, faidx, seq_samp_type, std_ref, ref2len)
        # make sure only ACGT chars in reference!
        if set(map_results.genome_seq).difference(valid_bases):
            yield None, "Non-ACGT sequence" # instead maybe just replace by random char?
            continue
        # map results
        map_results = map_results._replace(raw_signal=sig)
        try:
            # this causes sometimes TomboError: Read event to sequence alignment extends beyond bandwidth
            map_results = adjust_map_res(map_results, seq_samp_type, rsqgl_params)
            rsqgl_res = resquiggle.resquiggle_read(map_results, std_ref, rsqgl_params, outlier_thresh)
            n_iters = 1
            while n_iters < max_scaling_iters and rsqgl_res.norm_params_changed:
                rsqgl_res = resquiggle.resquiggle_read(map_results._replace(scale_values=rsqgl_res.scale_values),
                                                       std_ref, rsqgl_params, outlier_thresh)
                n_iters += 1
        except Exception as inst:
            yield None, str(inst)
            continue
        res = adjust_rsqgl_res(rsqgl_res, sig, seq_samp_type)
        if add_si:
            # get signal intensity means
            res.si = get_norm_mean(res.raw_signal, res.segs)
            # and modification probability by comparing observed si to expected for that sequence
            res.mp = get_mod_prob(res.align_info.refseq, res.si, std_ref, res.genome_loc.Strand)
            #res.si, res.mp = si, mp
        # add alignment and read as those are needed later
        res.a, res.sig, res.move, res.trace, res.md = a, sig, move, trace, md
        # update ref counter
        if ref2c: ref2c[contig] += 1
        yield res, ""

def get_norm_mean(raw, segs): 
    """Return raw signal means for given segments."""
    # ~2x faster than .mean()
    return np.array([raw[s:e].sum() for s, e in zip(segs[:-1], segs[1:])]) / (segs[1:]-segs[:-1]) # 25.6 ms ± 179 µs

def get_mod_prob(r_seq, r_means, std_ref, strand, fm_offset=0, SMALLEST_PVAL=1e-50):
    """Return modification probability calculated by comparing observed signal
    to the expected for given reference sequence. 

    If fm_offset > 0 additionally fisher method of P-value normalisation is used.

    This is equivalent to de_novo in tombo. 
    Note, tombo by default sets fm_offset=1, while here we stick to 0.
    Also r_seq needs to be already enlarged by k-mer (this is done earlier). 
    """
    r_ref_means, r_ref_sds = std_ref.get_exp_levels_from_seq(r_seq, strand=='-')  #get_exp_levels_from_seq_with_gaps()  
    # reverse means to match genomic order
    if strand == '-':
        r_means = r_means[::-1]
    # normalise signal using regression on expected signal
    #r_means = norm_signal_on_ref(r_means, r_ref_means)
    # get probability
    z_scores = np.abs(r_means - r_ref_means) / r_ref_sds
    r_pvals = scipy.stats.norm.cdf(-z_scores) * 2.0
    if fm_offset > 0:
        r_pvals = tombo_stats.calc_window_fishers_method(r_pvals, fm_offset)
    # ignore errors in max over NAN values if fisher's method was used
    #with np.errstate(invalid='ignore'):
    #    r_pvals = np.maximum(r_pvals, SMALLEST_PVAL)
    # reverse means back to read order
    #if strand == "-":
    #    r_means = r_means[::-1]
    return r_pvals #, r_means

def norm_signal_on_ref(r_means, ref_means):
    """Return r_means normalised to ref_means using ordinary least squares"""
    X = sm.add_constant(ref_means)
    est = sm.OLS(r_means, X).fit()
    a, b = est.params
    r_pred = a+b*r_means
    return r_pred
    

