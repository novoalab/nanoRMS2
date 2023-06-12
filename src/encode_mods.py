#!/usr/bin/env python3
desc="""Encode modifications probabilites predicted using de novo models into BAM files.

TO DO:
- 1+ modification per mer
- clean-up redundant functions
- make sure the functions between get_features and encode_mods are consistent
- no modifications on negative strand. why?
- replace multiprocessing with pebble
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Mizerów, 19/11/2020
"""

import joblib, os, resource, sys, numpy as np, mappy, pysam
from array import array
from datetime import datetime
from get_features import VERSION, logger, warning, base2complement, start_guppy_server, \
    get_basecall_client, init_args, get_sam_from_resquiggle
from predict_mods import moving_average
from multiprocessing import Pool
from scipy.ndimage.filters import maximum_filter1d

# set max 3 modifications per each base
MaxModsPerBase = 3
MaxProb = 256
QUALS = "".join(map(chr, range(33, 127)))
# features should be stored withing each model!
dt_shift = 10
features = ["si", "mp", "dt0", "dt%s"%dt_shift] + ["t%s"%b for b in "ACGT"]
f2idx = {f: i for i, f in enumerate(features)}

def get_all_mers(motif):
    mers = ['']
    # TODO add motifs with degenerated sequence ie CCWTT or DRACH
    for pos in motif:
        nmers = []
        for b in pos:
            for m in mers:
                nmers.append(m+b)
        mers = nmers
    return set(mers)

def get_modseq_from_clf(a, seq, mer2clf, can2mod, rna, mer2shift={}, k=7, nn=1, minp=0.5):
    """Return modified sequence - this is used by taiyaki"""
    # get data from tags and prepare storage
    tags = {k: v for k, v in a.tags}
    data = get_read_features(a, features, tags, f2idx, dt_shift, rna)
    # reverse rseq
    rseq = seq.upper()
    if a.is_reverse: 
        rseq = "".join(base2complement[b] for b in rseq[::-1])
        data = np.flip(data, axis=1)
    # and add 3 As at both ends seq
    n = (k//2)
    nseq = np.array(list(seq))
    rseq = n*"A" + seq + n*"A"
    #'''
    pvals = np.zeros(len(nseq))
    for i in range(nn, data.shape[1]-nn):
        mer = rseq[i:i+k]
        if mer in mer2clf:
            pvals[i] = mer2clf[mer].predict_proba(data[:, i-nn:i+nn+1].reshape(1, -1))[0, 1]
    # now make sure that the mer centre corresponds to modified base
    for i in np.argwhere(pvals>=minp)[:, 0]:
        # ideally narrow down to given pattern
        mer = rseq[i:i+k]#; print(mer)
        # update i with shift info for given mer
        if mer in mer2shift:
            i += mer2shift[mer]
        # proceed only if the position hasn't been marked as modified by their neighbours already
        if nseq[i] in can2mod:
            # store modified base taking into account possible shifts
            nseq[i] = can2mod[nseq[i]]
    return "".join(nseq)

def get_MaxPhredProb(MaxModsPerBase, QUALS=QUALS):
    MaxPhredProb = int(len(QUALS)/MaxModsPerBase)
    return MaxPhredProb

def load_models(models):
    """Return modification models"""
    # load models
    mer2clf = {}
    for f in models:
        # here we could store more predictor per mer in list ie 5mC & 5hmC for AACCAGG :)
        mer2clf.update(joblib.load(f))
    # get base2mod and mod2idx
    base2mod = {b: set() for b in base2complement.keys()}
    for m, clf in mer2clf.items():
        try: b, mod = clf.mod_info[:2]
        except:
            b = m[2]; mod = "mod%s"%b # dummy modification models - this is temporary
            clf.mod_info = (b, mod, ) # mers_info
        base2mod[b].add(mod)
    mod2idx = {m: i for base_mods in base2mod.values() for i, m in enumerate(base_mods)}
    return mer2clf, base2mod, mod2idx

def get_read_features(a, features, tags, f2idx, dt_shift, rna, dtype="float16"):
    """Return features from given read"""
    data = np.empty((len(features), len(tags[features[0]])), dtype=dtype)
    data[:] = -1 # store -1 instead of 0 (trace with -1 will be skipped) # solve missing trace at deletions in the read
    # here we need to add special treatment for dt!
    if "dt0" in f2idx or dt_shift:
        # normalise dwell time by moving average and store as log2
        dt = np.array(tags["dt"])
        dt = np.log2(dt / moving_average(dt)) #dt.mean())
    # store
    for j, k in enumerate(features, 0): 
        # dwell-time for position 0
        if k=="dt0": data[j] = dt
        # shifted dwell time
        ## TBD: make it strand aware!
        elif k.startswith("dt"):
            if rna and not a.is_reverse or not rna and a.is_reverse: 
                if data.shape[1]>dt_shift: # make sure enough alignment here # and len(dt[s+dt_shift:e]):
                    data[j, :-dt_shift] = dt[dt_shift:]
            elif data.shape[1]>dt_shift: # and here as well len(dt[s:e-dt_shift]): 
                data[j, dt_shift:] = dt[:-dt_shift]
        # normalise trace
        elif k.startswith("t"): data[j] = np.array(tags[k], dtype=dtype)/255
        # and remaining features
        else: data[j] = tags[k]
    return data

def get_sequence_and_cigar(a, tags, faidx):
    """Return reference sequence and cigar string"""
    # get reference sequence (spliced for RNA)
    if "bs" in tags and len(tags["bs"])>2:
        # process alignment blocks and store cigar and seq
        seq, cigar, pe = [], [], 0
        blocks = np.array(tags["bs"]).reshape(-1, 2)
        for s, e in blocks: # NEED TESTING
            if pe: cigar.append("{}N".format(s-pe, )) # store intron cigar
            cigar.append("{}M".format(e-s,)) # store exon cigar
            seq.append(faidx.fetch(a.reference_name, s, e)) # store exon seq
            pe = e # update previous exon end
        # join exonic sequence and cigar entries
        seq, cigar = "".join(seq), "".join(cigar)
    else:
        seq = faidx.fetch(a.reference_name, a.reference_start, a.reference_end)
        cigar = "{}M".format(len(seq))
    return seq, cigar

def encode_mod_prob_as_quals(rseq, data, mer2clf, mod2idx, MaxPhredProb, nn, k):
    """Return modification probabilities encoded as PHRED quality scores"""
    tested = nmods = 0
    quals = np.zeros(data.shape[1], dtype="int")
    for i in range(nn, data.shape[1]-nn):
        # rseq is already extended with 3 bases on both sides
        mer = rseq[i:i+k] # this is ok, but why shouldn't this be shifted by -1?
        if mer in mer2clf: # here a 7-mer can have more than 1 modification ie AACCAGG can by 5mC and 5hmC
            tested += 1
            # get clf
            clf = mer2clf[mer]
            # get mod_probability
            # RF: 4.17 ms ± 10.6 µs; GradientBoostingClassifier is ~10x faster
            p = clf.predict_proba(data[:, i-nn:i+nn+1].reshape(1, -1))[0, 1] # 10% faster than .flatten()[np.newaxis, :]
            if p>0.5: nmods += 1 #(p>=0.5).sum()
            # and store Phred-scaled probability 
            quals[i] = int(round(p*MaxPhredProb))+mod2idx[clf.mod_info[1]]*MaxPhredProb
    return quals, nmods, tested

def get_alg_with_mod_prob(a, faidx, mer2clf, mod2idx, features, f2idx, dt_shift, rna,
                          MaxPhredProb, nn, k):
    """Return alignment object with modification probability encoded as base qualities
    and updated sequence and cigar. 
    All tags are stripped by default.
    """
    # get data from tags and prepare storage
    tags = {k: v for k, v in a.tags}
    data = get_read_features(a, features, tags, f2idx, dt_shift, rna)
    seq, cigar = get_sequence_and_cigar(a, tags, faidx)    
    # reverse rseq
    rseq = seq.upper()
    if a.is_reverse: 
        rseq = "".join(base2complement[b] for b in rseq[::-1])
        data = np.flip(data, axis=1)
    # and add 3 As at both ends
    rseq = (k//2)*"A" + rseq + (k//2)*"A"
    # encode modifications as base qualities
    quals, nmods, tested = encode_mod_prob_as_quals(rseq, data, mer2clf, mod2idx, MaxPhredProb, nn, k)
    # reverse qualities if rev-complement
    if a.is_reverse: quals = quals[::-1]
    # store seq, qualities and fixed cigar
    a.query_sequence, a.query_qualities, a.cigarstring = seq, array("B", quals), cigar
    # strip tags other than polyA tail length estimation
    a.tags = [(k, v) for (k, v) in a.tags if k.startswith("a")]
    return a, nmods, tested

def bam2encode(outdir, indir, fasta, rna, models, nn, k, dt_shift, mapq, chrom="",
               debug=False, conf=("", "", "")):
    """Predict modifications using de novo model and features from BAM file
    and store PHRED encoded modifications in new BAM file
    """
    bam = indir.rstrip(os.path.sep) + ".bam"
    # open out BAM file
    outfn = os.path.join(outdir, os.path.basename(bam))
    if os.path.isfile(outfn):
        logger(" file exists: %s"%outfn)
        return outfn
    # load models
    mer2clf, base2mod, mod2idx = load_models(models)    
    t0 = datetime.now() # features should be stored withing each model!
    features = ["si", "mp", "dt0", "dt%s"%dt_shift] + ["t%s"%b for b in "ACGT"]
    f2idx = {f: i for i, f in enumerate(features)}
    # get maxPhredProb
    MaxPhredProb = get_MaxPhredProb(MaxModsPerBase) #31
    # open FastA and alignment file
    faidx = pysam.FastaFile(fasta)
    if os.path.isfile(bam):
        sam = pysam.AlignmentFile(bam)
        header = sam.header.as_dict()
        sort = False
        mapped = sum(s.mapped for s in sam.get_index_statistics() if not chrom or chrom and s.contig==chrom)
    else:
        sam = get_sam_from_resquiggle(indir, fasta, rna, conf)
        header = pysam.AlignmentHeader.from_references(faidx.references, faidx.lengths).to_dict()
        sort = True
        mapped = "unknown"
        chrom = "" # disable chromosome fetching - it can't be done here
    # strip PG fiels
    if "PG" in header: header.pop("PG")
    # and add mod info
    header["CO"] = ["%s:%s"%(b, " ".join(m)) for b, m in base2mod.items() if m]
    out = pysam.AlignmentFile(outfn, "wb", header=header)
    # process alignment
    ri = tot_bases = tot_mods = tot_tested = 0
    for ri, a in enumerate(sam.fetch(chrom) if chrom else sam, 1):
        if debug and ri>1000: break
        dt = datetime.now()-t0
        if a.mapq<mapq: continue
        # get alg object with mod_prob encoded as baseQ, fixed seq and cigar
        a, nmods, tested = get_alg_with_mod_prob(a, faidx, mer2clf, mod2idx, features, f2idx, dt_shift, rna, MaxPhredProb, nn, k)
        tot_bases += len(a.seq)
        tot_mods += nmods
        tot_tested += tested
        if not ri%1000 and dt.seconds:
            sys.stderr.write(" {} / {} {:.1f} kbps; {:.3f}% tested mers modified  \r".format(ri, mapped, tot_bases/1000/dt.seconds, 100*tot_mods/tot_tested))
        # and store
        out.write(a)
    out.close()
    # sort if needed (live basecalling)
    if sort:
        os.rename(outfn, outfn+".unsorted")
        pysam.sort("-o", outfn, outfn+".unsorted")
        os.unlink(outfn+".unsorted")
    # index outfn
    _ = pysam.index(outfn)
    logger(" {} with {:,} bases ({:.3f}% tested mers modified) in {:,} alignments".format(bam, tot_bases, 100*tot_mods/tot_tested, ri))
    return outfn

def mod_encode(outdir, indirs, fasta, models, threads=1, rna=True, mapq=15, chrom="",
               config="", host="", port=5555, device="cuda:0", 
               dt_shift=10, nn=1, k=7, debug=False):
    """Process multiple directories from Fast5 files"""
    # start guppy server if needed
    conf = (config, host, port)
    if host:
        guppy_proc, host, port = start_guppy_server(host, config, port, device)
        # try connecting to guppy server first
        conf = (config, host, port)
        client, _get_read_data = get_basecall_client(*conf)
    else:
        guppy_proc = None
        logger("We'll use basecall information from Fast5 files...")
    # load reference for mappy
    logger("Loading reference index from %s..."%fasta)
    aligner = mappy.Aligner(fasta, preset="spliced" if rna else "map-ont", k=13)
        
    if not os.path.isdir(outdir): os.makedirs(outdir)    
    logger("Encoding modification probabilities in %s BAM file(s)..."%len(indirs))    
    # process input BAM files #mer2clf, base2mod, mod2idx
    args = [(outdir, fn, fasta, rna, models, nn, k, dt_shift, mapq, chrom, debug, conf)
            for fn in indirs]
    #if debug: bams = [bam2encode(*arg) for arg, bam in zip(args, indirs)]
    # else:
    # start Pool of processes
    p = Pool(threads if threads < len(indirs) else len(indirs), maxtasksperchild=1,
             initializer=init_args, initargs=(aligner, ))        
    bams = [bam for bam in p.starmap(bam2encode, args)]
    # close pool
    p.close(); p.join()
    # and guppy basecall server if it was started by this process
    if guppy_proc: guppy_proc.terminate()
    return bams

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-o", "--outdir", required=1, help="output directory")
    parser.add_argument("-i", "--indirs", nargs="+", help="input directory with Fast5 files")
    #parser.add_argument("-i", "--input", nargs="+", help="input BAM file(s)")
    parser.add_argument("--rna", action='store_true', help="project is RNA sequencing [DNA]")
    parser.add_argument("-f", "--fasta", required=1, help="reference FASTA file")
    parser.add_argument("-m", "--models", nargs="+", help="de novo model file(s)")
    parser.add_argument("-t", "--threads", default=6, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("--mapq", default=15, type=int, help="min mapping quality [%(default)s]")
    parser.add_argument("--debug", action="store_true", help="debugging mode")
    parser.add_argument("--chr", help="process only one chromosome [all]")
    guppy = parser.add_argument_group("Basecalling options") #mutually_exclusive
    guppy.add_argument("-c", "--config", default="rna_r9.4.1_70bps_ivt_fast.cfg", help="guppy model [%(default)s]")
    guppy.add_argument("--host", "--guppy_basecall_server", default="",
                        help="guppy server hostname or path to guppy_basecall_server binary [use basecall information from Fast5]")
    guppy.add_argument("--port", default=5555, type=int,
                        help="guppy server port (this is ignored if binary is provided) [%(default)s]")
    guppy.add_argument("--device", default="cuda:0", help="CUDA device to use [%(default)s]")
    #guppy.add_argument("--timeout", default=60*30, help="timeout in seconds to process each Fast5 file [%(default)s]")

    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    bamfiles = mod_encode(o.outdir, o.indirs, o.fasta, o.models, o.threads, o.rna, o.mapq, o.chr,
                          o.config, o.host, o.port, o.device, debug=o.debug)
        
if __name__=='__main__': 
    t0 = datetime.now()
    os.setpgrp() # create new process group, become its leader    
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except Exception as err:
        import signal, traceback
        sys.stderr.write(traceback.format_exc()+"\n")
        os.killpg(0, signal.SIGTERM) # terminate all processes in my group
    finally:
        dt = datetime.now()-t0
        sys.stderr.write("#Time elapsed: %s    \n"%dt)


############################################################
# BAM processing functions
def fasta2bases(fastafn, ref, start, end, strands="+-", n=3):
    """Generator of individual bases from FastA file.

    The output consists of: 
    - position in reference (1-based)
    - strand integer (0 for plus, 1 for minus)
    - strand as +/-
    - base (complement for -)
    - 7-mer (+/- n bases surrounding given position)
    """
    fasta = pysam.FastaFile(fastafn)
    ref2len = {r: l for r, l in zip(fasta.references, fasta.lengths)}
    if ref not in ref2len: #fasta.references:
        raise StopIteration
    for pos, refbase in enumerate(fasta.fetch(ref, start, end), start+1):
        refbase = refbase.upper()
        # combine before start NNN (if needed) sequence from ref and after start NNN (if needed)
        mer = "N"*(n-pos+1) + "".join(fasta.fetch(ref, pos-n-1 if pos>n+1 else 0, pos+n)) + "N"*(pos-ref2len[ref]+n)
        mer = mer.upper() # need to be upper case
        for si, strand in enumerate(strands):
            if si:
                refbase = base2complement[refbase]
                mer = "".join(base2complement[b] for b in mer[::-1])
            yield pos, si, strand, refbase, mer

def bam2data(bam, ref, start, end, maxDepth=100000, nn=1, dt_shift=10, rna=True, n=6):
    """
    Generator of data for consecutive positions from ref:start-end region
    More memory efficient.
    Having same dtype is crucial as it speeds up unloading int16 ~20-80% (int32 - float/double).
    """
    sam = pysam.AlignmentFile(bam)#; print(ref, start, end)
    # define dimensions
    dim = (maxDepth, end-start) 
    # sis, tr, baseQ, mod_prob_tombo, dts dts_shifted
    data = np.zeros((n, maxDepth, end-start), dtype="float32")
    strands = np.zeros(dim, dtype="int8") # 0: +/FOR; 1: -/REV; -1: no alignment/read
    strands[:] = -1
    readi = 0
    s_keys, t_keys, m_keys = [], [], []
    for a in sam.fetch(ref, start, end):
        # make sure first position of read always corresponds to first pos in data
        while a.pos>start:
            for strand in range(2):
                flt = strands[:readi, nn] == strand
                # report
                yield data[:, :readi][:, flt, :2*nn+1] #[d[:readi][flt, :2*nn+1] for d in data]
            # strip position 0 from arrays
            data = data[:, :, 1:] #[d[:, 1:] for d in data]
            strands = strands[:, 1:]
            start += 1
        # define read start & end and strand (1>-, 0>+)
        s, e = start-a.pos, a.aend-a.pos if a.aend<=end else end-a.pos
        # and region end
        re = e-s 
        #print("%s-%s"%(start, end), a.pos, a.aend, s, e)
        strands[readi][:re] = 1 if a.is_reverse else 0
        # get tags
        tags = {k: v for k, v in a.tags}
        data[0][readi, :re] = np.array(tags["si"][s:e]) / 1000
        data[1][readi, :re] = np.array(tags["tr"][s:e]) / 255
        data[2][readi, :re] = a.query_qualities[s:e]
        data[3][readi, :re] = tags["mp"][s:e]
        # normalise dwell time by moving average and store as log2
        dt = np.array(tags["dt"])
        dt = np.log2(dt / moving_average(dt)) # get moving average
        data[-2][readi, :re] =  dt[s:e] #np.log2(dt[s:e] / dt[s:e].mean())
        # shift dt by 10 bases
        if e-s+dt_shift>0:
            if rna and not a.is_reverse or not rna and a.is_reverse: 
                if len(dt[s+dt_shift:e]):
                    data[-1][readi, :re-dt_shift] = dt[s+dt_shift:e] #np.log2(dt[s+dt_shift:e] / dt[s+dt_shift:e].mean())
            elif len(dt[s:e-dt_shift]):
                data[-1][readi, :re-dt_shift] = dt[s:e-dt_shift] #np.log2(dt[s:e-dt_shift] / dt[s:e-dt_shift].mean())
        readi += 1
        if readi>=maxDepth:
            logger("[WARNING] maxDepth reached for %s:%s-%s @ %s"%(ref, start, end, bam))
            break
    # report last bit from region
    for pos in range(strands.shape[1]):
        for strand in range(2):
            flt = strands[:readi, pos+nn] == strand
            # report
            yield data[:, :readi][:, flt, pos:pos+2*nn+1] #[d[:readi][flt, pos:pos+2*nn+1] for d in data]
