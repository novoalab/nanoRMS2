#!/usr/bin/env python3
desc="""Compare two samples (KO/IVT/PCR and wt) using features encoded in BAM files
and report likely modified positions. 

"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Mizerów/Barcelona/Cap Salou, 23/11/2020
"""
import matplotlib; matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
import csv, gzip, h5py, joblib, os, pysam, scipy, sys, warnings
import matplotlib.pyplot as plt, numpy as np, pandas as pd, seaborn as sns
import xml.etree.ElementTree as ET
from datetime import datetime
from scipy.stats import stats
from sklearn.neighbors import KNeighborsClassifier
from sklearn.semi_supervised import LabelPropagation
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from multiprocessing import Pool
from pebble import ProcessPool, ProcessExpired
from get_features import get_features, VERSION, logger, warning, base2complement, get_results_from_pool
# skip convergence warnings
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings(action='ignore', category=ConvergenceWarning)

alphabet = "NACGT" # first base should be always N!!!
base2index = {b: i for i, b in enumerate(alphabet)}
for i, b in enumerate(alphabet.lower()):
    base2index[b] = i

# init random with seed, so results are deterministic
random_state = 12345
rng = np.random.default_rng(random_state)
old_settings = np.seterr(under='ignore') # only needed for scaling

# disable sklearn paralelism: https://scikit-learn.org/stable/computing/parallelism.html#parallelism
#env_vars = ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "BLIS_NUM_THREADS"]
#for v in env_vars: os.environ[v] = "1"
# doesn't work https://joblib.readthedocs.io/en/latest/parallel.html#avoiding-over-subscription-of-cpu-resources
# disabling through joblib also seems not working

def get_KS_for_one_mer(mer_data):
    """Calculate KS statistics for data from given mer"""
    # store KS and P-value
    ks = np.zeros((mer_data[0].shape[1], 2))
    # ignore all errors
    with np.errstate(all='ignore'): 
        for fi in range(mer_data[0].shape[1]):
            ks[fi] = [*stats.ks_2samp(*[a[:, fi] for a in mer_data])]
    # return flattened array: KS, P for each feature
    return tuple(ks.flatten())

def get_KS(mer2data, threads=6, verbose=True, p=None):
    """Return KS statistics for all mers for all features"""
    # single-thread version also uses lots of memory!
    #return [(mer, *get_KS_for_one_mer(mer2data[mer][:2])) for mer in mer2data]
    t0 = datetime.now()
    ks = []
    # get sorted mers
    all_mers = list(sorted(mer2data.keys()))
    args = (mer2data[mer][:2] for mer in all_mers)
    # process all mers
    '''
    p = ProcessPool(max_workers=threads)#, max_tasks=1000)
    future = p.map(get_KS_for_one_mer, args, chunksize=10)#, timeout=timeout)
    for i, r in enumerate(get_results_from_pool(future, all_mers)):
    '''
    if not p: p = Pool(threads)
    for i, r in enumerate(p.imap(get_KS_for_one_mer, args)):    
        if verbose and not i%100: sys.stderr.write(" {:,} / {:,} {}\r".format(i, len(all_mers), datetime.now()-t0))
        ks.append((all_mers[i], *r)) # and store mer and KS stats for each feature
    # close pool of workers; skip join as it takes lots of memory
    #p.close()
    return ks

def init_args(*args):
    """Share globals with pool of workers"""
    global sams, faidx
    bams, fasta = args
    sams = [pysam.AlignmentFile(bam) for bam in bams]
    faidx = pysam.FastaFile(fasta)

def load_data_worker(args): #ref, start, end, fasta, bams, features, nn=1, maxReads=10000, rna=True):
    ref, start, end, fasta, bams, features, nn, maxReads, rna, minReads, mers, mapq = args
    mers = set(mers)
    mer2data = {}
    if start<nn: start=nn
    # this is global value initialized by init_args - this saves ~80ms per region on SSD!
    #sams = [pysam.AlignmentFile(bam) for bam in bams] 
    parsers = [bam2data(sam, ref, start-nn, end+nn+1, rna, nn, features, mapq=mapq, 
                        maxDepth=maxReads*3, verbose=False) for sam in sams]
    # faidx is global var as well
    #faidx = pysam.FastaFile(fasta)
    refparser = fasta2bases(faidx, ref, start, end) # report 7 or 9-mers
    for ((pos, _, strand, refbase, mer), *calls) in zip(refparser, *parsers):
        if "N" in mer or mers and mer not in mers: continue
        # here we'd need to flip -1, 0, +1 for strand -
        if strand=="-": calls = [np.flip(c, axis=2) for c in calls]
        sample2data = [np.hstack(c) for c in calls]
        cov = tuple(map(len, sample2data))
        max_reads = min(cov)#; print(ref, pos, strand, mer, max_reads)
        if max_reads<minReads: continue # this is faster since only one sample has choice done
        elif max_reads>maxReads: max_reads = maxReads
        #'''
        # this takes ~25% of time - let's skip random selection and just take first N reads?
        #in 6-thread run loading 12m vs 11m
        merdata = [sd if len(sd)==max_reads
                   else sd[rng.choice(np.arange(len(sd)), max_reads, replace=False)]
                   for sd in sample2data]
        '''
        merdata = [sd[:max_reads] for sd in sample2data]
        #'''
        pos_info = "{}\t{}\t{}\t{}".format(ref, pos, strand, max_reads)
        if mer not in mer2data: # create new entry for mer
            mer2data[mer] = [*[[merdata[i]] for i in range(len(cov))], [pos_info]]
        else: # add new position for existing mer
            for i in range(len(cov)): mer2data[mer][i].append(merdata[i])
            mer2data[mer][-1].append(pos_info)
    return mer2data

def load_data(fasta, bams, features, nn=1, rna=True, regions=[], minReads=10, threads=6, mers=[], 
              maxReads=100, mapq=15, maxtasksperchild=100, p=None):
    """Return features for positions of interest"""
    mer2data = {}
    t0 = datetime.now()
    args = ((ref, start, end, fasta, bams, features, nn, maxReads, rna,
             minReads, mers, mapq) for ref, start, end in regions)
    '''
    p = ProcessPool(max_workers=threads, #max_tasks=maxtasksperchild, # this is actually worse
                    initializer=init_args, initargs=(bams, fasta))
    future = p.map(load_data_worker, args)#, timeout=timeout)
    for ri, md in enumerate(get_results_from_pool(future, regions), 1):
    '''
    if not p: p = Pool(threads, initializer=init_args, initargs=(bams, fasta))
    for ri, md in enumerate(p.imap(load_data_worker, args), 1): 
        if not ri%10: sys.stderr.write(" %s / %s  %s\r"%(ri, len(regions), datetime.now()-t0))
        for mer, merdata in md.items(): #(ko, wt, pos)
            # create new entry for mer or add new position for existing mer
            if mer not in mer2data:
                mer2data[mer] = merdata
            else:
                for i in range(len(merdata)):
                    mer2data[mer][i] += merdata[i]
            # memory overhead will be ~0.3MB for 1000 elements for each mer (6GB in total)
            # it may cause memory exploding after vstack, but can't understand fully why
            if len(mer2data[mer][0])>1000:
                md = mer2data[mer]
                mer2data[mer] = [[np.vstack(md[i]), ] for i in range(len(bams))] + [[";".join(md[-1]), ]]
            #'''
    # final combine lists of elements
    for mer in mer2data:
        md = mer2data[mer]
        mer2data[mer] = [np.vstack(md[i]) for i in range(len(bams))] + [";".join(md[-1])]
    # close pool of workers; skip join as it takes lots of memory
    #p.close()
    return mer2data

def load_data_single(fasta, bams, features, nn=1, rna=True, regions=[], minReads=10,
                     mers=[], maxReads=100, mapq=15):
    """Return features for positions of interest"""
    mers = set(mers)
    t0 = datetime.now()
    mer2data = {}
    faidx = pysam.FastaFile(fasta)
    sams = [pysam.AlignmentFile(bam) for bam in bams]
    for ri, (ref, start, end) in enumerate(regions, 1):
        if not ri%10: sys.stderr.write(" %s / %s  %s\r"%(ri, len(regions), datetime.now()-t0))
        if start<nn: start=nn
        parsers = [bam2data(sam, ref, start-nn, end+nn+1, rna, nn, features, mapq=mapq, 
                            maxDepth=maxReads*3, verbose=False) for sam in sams]
        refparser = fasta2bases(faidx, ref, start, end) # report 7-mers
        for ((pos, _, strand, refbase, mer), *calls) in zip(refparser, *parsers):
            if "N" in mer or mers and mer not in mers: continue
            # here we'd need to flip -1, 0, +1 for strand -
            if strand=="-": calls = [np.flip(c, axis=2) for c in calls]
            sample2data = [np.hstack(c) for c in calls]
            cov = tuple(map(len, sample2data))
            max_reads = min(cov)
            if max_reads<minReads: continue # this is faster since only one sample has choice done
            elif max_reads>maxReads: max_reads = maxReads
            #'''
            # this takes ~25% of time - let's skip random selection and just take first N reads?
            merdata = [sd if len(sd)==max_reads
                       else sd[rng.choice(np.arange(len(sd)), max_reads, replace=False)]
                       for sd in sample2data]
            '''
            merdata = [sd[:max_reads] for sd in sample2data]
            #'''
            pos_info = "{}\t{}\t{}\t{}".format(ref, pos, strand, max_reads)
            # create new entry for mer or add new position for existing mer
            if mer not in mer2data:
                mer2data[mer] = [*[[merdata[i]] for i in range(len(cov))], [pos_info]]
            else:
                for i in range(len(cov)): mer2data[mer][i].append(merdata[i])
                mer2data[mer][-1].append(pos_info)
    # combine lists of elements
    for mer in mer2data:
        md = mer2data[mer]
        mer2data[mer] = [np.vstack(md[i]) for i in range(len(bams))] + [";".join(md[-1])]
    return mer2data

def load_bed(fname):
    """Return regions from BED file. If not a file, try to unload region(s) from a string"""
    regions = []
    if os.path.isfile(fname) or os.path.islink(fname):
        for l in open(fname):
            if l.startswith('#') or not l[:-1]:
                continue
            ldata = l[:-1].replace(',','').split('\t')#; print ldata
            if len(ldata) >= 3:
                ref, start, end = ldata[:3]
            else:
                ref, se = ldata[0].split(':')
                start, end = se.split('-')
            start, end = map(int, (start, end))
            regions.append((ref, start, end)) #yield ref, start, end
    else:
        for region in fname.split():
            if not region: continue
            ref, se = region.replace(',','').split(':')
            start, end = se.split('-')
            start, end = map(int, (start, end))
            regions.append((ref, start, end))
    # & split regions into 1kb windows!
    return regions

def get_revcomp(bases):
    """Return reverse comlement"""
    return "".join(base2complement[b] for b in bases[::-1])

def fasta2bases(faidx, ref, start, end, n=3):
    """Generator of individual bases from FastA file.

    The output consists of: 
    - position in reference (1-based)
    - strand integer (0 for plus, 1 for minus)
    - strand as +/-
    - base (complement for -)
    - 7-mer (+/- n bases surrounding given position)
    """
    if ref not in faidx: raise StopIteration
    # get uppercase sequence extended with N if needed at the ends
    k = 2*n+1
    seq = faidx.fetch(ref, start-n if start>=n else 0, end+n).upper()
    seq = "N"*(n-start)+seq+"N"*n
    cseq = "".join(base2complement[b] for b in seq) # doing this speeds-up everything 2x
    for pos, i in enumerate(range(end-start), start+1):
        mer = seq[i:i+k]
        rc_mer = cseq[i:i+k][::-1]
        # yield positive and negative strand
        yield pos, 0, "+", mer[n], mer
        yield pos, 1, "-", rc_mer[n], rc_mer

#@numba.jit # this is 4x faster
def moving_average(a, n=5):
    """Calculate moving average including first n-1 objects"""
    ret = np.cumsum(a, dtype=float) # for numba remove dtype and turn into float64 before
    ret[n:] = ret[n:] - ret[:-n]
    ret[n-1:] *= 1 / n
    ret[:n-1] *= 1 / np.arange(1, n)
    return ret

def bam2data(sam, ref, start, end, rna=True, nn=1, features=["si", "tr"],
             maxDepth=1000, mapq=20, dtype="float16", verbose=1):
    """Generator of data for consecutive positions from ref:start-end region"""
    #sam = pysam.AlignmentFile(sam)#; print(ref, start, end)
    # we can skip the region if few reads
    #sam.count(ref, start, end)
    # get dt_shift
    f2idx = {f: i for i, f in enumerate(features)}
    dt_shift_keys = list(filter(lambda k: k.startswith("dt") and k!="dt0", f2idx.keys()))
    dt_shift = 0 if not len(dt_shift_keys) else int(dt_shift_keys[0][2:]) # dt10 > 10
    # update end position with shift
    end += dt_shift # here for DNA it's a bit complicated as we'd need to do start-=dt_shift
    # this is needed later
    requested_tags = list(filter(lambda f: not f.startswith("dt"), features))
    if dt_shift or "dt0" in features: requested_tags.append("dt")
    # here only SI & MP # here dt_shift should be read from the feature
    data = np.empty((len(features), maxDepth, end-start), dtype=dtype)
    # solve missing trace at deletions in the read
    data[:] = -1 # store -1 instead of 0 (trace with -1 will be skipped)
    strands = np.zeros((maxDepth, end-start), dtype="int8") # 1: +/FOR; -1: -/REV; 0: no alignment/read
    readi = 0
    for a in sam.fetch(ref, start, end):
        # filter by mapq
        if is_qcfail(a, mapq): continue
        # make sure first position of read always corresponds to first pos in data
        while a.pos>start: # consider skipping first/last 5-15 bases
            # report data for reads from + & - strand separately
            for strand in (1, -1):
                flt = strands[:readi, nn] == strand
                yield data[:, :readi][:, flt, :2*nn+1] 
            # strip position 0 from arrays
            data = data[:, :, 1:]
            strands = strands[:, 1:]
            start += 1
        # get data from tags
        tags = {k: v for k, v in a.tags}
        # skip read if frac of positions that have P-value reported by tombo < 0.01
        if (np.array(tags['mp']) < 0.01).mean() > 0.3: continue
        # define read start & end for current data view
        s, e = start-a.pos, a.aend-a.pos if a.aend<=end else end-a.pos
        # and region end
        re = e-s
        # turn exonic blocks back into genomic coordinates
        if "bs" in tags and len(tags["bs"])>2:
            # get blocks as 2D array (start-end) with exonic intervals of the read
            blocks = np.array(tags["bs"]).reshape(-1, 2)-tags["bs"][0]
            # take care only about requested features
            _tags = {}
            for f in requested_tags:
                # storing 1s is importand as dt is log2 obs/exp, thus it can't be 0s
                _tags[f] = np.ones(a.reference_length, dtype=dtype)
            # mark exonic block in strands
            read_strands = np.zeros(a.reference_length, dtype="int8")
            # store block info
            pe = 0
            for bs, be in blocks:
                # mark exonic block in read_strands
                read_strands[bs:be] = -1 if a.is_reverse else 1
                # store exon block into new tags
                blen = be-bs
                for f in requested_tags:
                    #print(f, bs, be, be-bs, pe, be-pe)
                    _tags[f][bs:be] = tags[f][pe:pe+blen]
                pe += blen
            # replace tags & udpate exonic strands
            tags = _tags
            strands[readi, :re] = read_strands[s:e]
        else:
            # mark entire read as stand
            strands[readi, :re] = -1 if a.is_reverse else 1
        # here we need to add special treatment for dt!
        if "dt0" in f2idx or dt_shift:
            # normalise dwell time by moving average and store as log2
            dt = np.array(tags["dt"])
            dt = np.log2(dt / moving_average(dt)) #dt.mean())
        # store
        for j, k in enumerate(features, 0): #for k, j in f2idx.items(): #
            # dwell-time for position 0
            if k=="dt0": data[j, readi, :re] = dt[s:e]
            # shifted dwell time
            elif k.startswith("dt"):
                if rna and not a.is_reverse or not rna and a.is_reverse: 
                    if e>s+dt_shift: # make sure enough alignment here # and len(dt[s+dt_shift:e]):
                        data[j, readi, :re-dt_shift] = dt[s+dt_shift:e]
                elif e-dt_shift>s: # and here as well len(dt[s:e-dt_shift]): 
                    data[j, readi, :re-dt_shift] = dt[s:e-dt_shift]
            # normalise trace - this isn't needed cause we'll do min-max anyway
            elif k.startswith("t"): data[j, readi, :re] = np.array(tags[k][s:e], dtype=dtype)/255
            # and remaining features
            else:
                data[j, readi, :re] = tags[k][s:e]
        readi += 1
        # clean-up only if maxDepth reached
        if readi>=maxDepth:
            if verbose: logger("[INFO] maxDepth reached for %s:%s-%s @ %s"%(ref, start, end, bam))
            # get algs that still carry data
            ## here read has strand over from 0 to end (not excluding introns)
            ne = strands[:, 0]!=0 # np.all(strands!=0, axis=0)#?
            readi = ne.sum() # update readi
            if readi>=maxDepth: # if all reads are still aligned, remove random 25% of algs
                ne[rng.integers(0, len(ne), int(0.25*maxDepth))] = False
                readi = ne.sum() # update readi
            # update strands & data
            _strands, _data = np.zeros_like(strands), np.zeros_like(data)
            _strands[:readi] = strands[ne]
            _data[:, :readi] = data[:, ne]
            strands, data = _strands, _data
    # report last bit from region
    for pos in range(strands.shape[1]-nn):
        # report data for reads from + & - strand separately
        for strand in (1, -1):
            flt = strands[:readi, pos+nn] == strand
            yield data[:, :readi][:, flt, pos:pos+2*nn+1]

def min_max_norm(X):
    """Return (X-min(X))/(max(X)-min(X))"""
    Xmax, Xmin = X.max(axis=0), X.min(axis=0)
    sel = Xmin!=Xmax
    if sel.sum():
        X[:, sel] = (X[:, sel] - Xmin[sel])/(Xmax[sel] - Xmin[sel]) # here if min==max div by 
    return X

def get_kmer(faidx, chrom, pos, strand, r2l, k=7):
    """Return base and k-mer centering at the positions of interest."""
    # positions are 0-based, right?
    sseq = eseq = ""
    if pos<k: 
        s = 0 
        sseq = "N"*(k-pos)
    else:
        s = pos-k
    if pos>r2l[chrom]-k:
        e = r2l[chrom]
        eseq = "N"*(r2l[chrom]-k+1)
    else:
        e = pos+k+1
    seq = sseq+faidx.fetch(chrom, s, e)+eseq
    seq = seq.upper() # get rid of lowercase
    if strand=="-": seq = get_revcomp(seq)
    return seq[k], seq

def get_freq(y_pred, cov):
    freq = []
    ps = 0
    for c in cov:
        freq.append(y_pred[ps:ps+c].mean())
        ps+=c
    return freq

def get_classifier(alldata, sel, clf0, clf, labelProp, fnames, probath=0.8):
    """Retrun semi-supervised classifier trained on sel examples."""
    # here we could actually limit the training set size if it's too large!
    ## as it'll take ages to get the classifier for larger datasets
    # get X and y
    data = [d[sel] for d in alldata]
    cov = list(map(len, data))
    X = np.vstack(data)
    y = np.zeros(len(X), dtype="int8") # KO
    y[cov[0]:] = 1 # WT - here many may be unmodified
    # train dirty classifier
    clf0.fit(X, y) # here we train and predict on the same dataset (dirty classifier) # 2ms
    proba = clf0.predict_proba(X) # get probability of being 1 (modified) # 9.3s for 24k
    # establish certain labels
    y[len(data[0]):][proba[len(data[0]):, 1]<probath] = -1 # store all WT with prob1 < threshold as unknown (-1)
    if (y==0).sum()==0 or (y==1).sum()==0:
        return None # (mer, None), []
    # infer unknown (uncertain) labels
    with np.errstate(all='ignore'): 
        y2 = labelProp.fit(X, y).transduction_ # infer unknown labels # 9.3s for 24k
    # here actually we should fix y2 from KO to 0 again!
    #if KO: y2[:len(y2)//2] = 0
    # fit final classifier using inferred (semi-supervised) labels
    clf.fit(X, y2) # GBC 7s vs HistGBC 538 ms, but HistGBC is very slow on single sample (like RF)
    return clf

def get_calls(args, clf0=KNeighborsClassifier(n_jobs=1), clf=GradientBoostingClassifier(random_state=random_state), # RandomForestClassifier(),
              labelProp=LabelPropagation(kernel='knn', max_iter=1000, n_jobs=1),
              probath=0.8, train_frac=0.5, max_train_samples=10000):
    """Return modification frequency predicted by GBC trained on semi-supervised labels. 
    
    By default (KO==True) only modified reads from WT are propagated, while 
    all reads from KO are forced to be unmodified. That's safer. 
    If KO==False, LabelPropagation estimates both, unmodified reads in KO and modified reads in WT.
    In addition, we'll try to guess which group is less modified if there is trace (TR) information. 
    Otherwise, first sample is always assumed to have less modified reads. 
    """
    # unload variables
    (fasta, mer, mer_data, fnames, minModDiff) = args
    alldata, pos_data = mer_data[:2], mer_data[2]
    # need at least several samples to work with...
    if len(alldata[0])<10:
        logger("[WARNING] Very few samples for %s: %s"%(mer, len(alldata[0])))
        return (mer, None), []
    # init faidx
    faidx = pysam.FastaFile(fasta)
    r2l = {r: l for r, l in zip(faidx.references, faidx.lengths)}
    # get semi-supervised classifier
    # select randomly 50% of the reads - here we could just take odd/even reads, right?
    sel = np.zeros(len(alldata[0]), dtype="bool")
    step = 2*len(sel)//max_train_samples
    if step < 2: step = 2
    sel[::step] = True
    # here we could use random 5k reads for large datasets or half of the samples if N<10k
    clf = get_classifier(alldata, sel, clf0, clf, labelProp, fnames, probath)
    if clf==None:
        return (mer, None), []
    # and evaluate on test set
    data_test = [d[~sel] for d in alldata]
    cov = list(map(len, data_test))
    X_test = np.vstack(data_test)
    freq3 = get_freq(clf.predict(X_test), cov)
    clf.freqs = freq3
    mod_diff = np.abs(np.diff(freq3))[0]
    # skip if difference in mod_freq between KO & WT for given mer is < than threshold
    ## here we could develop something to at least recover a bunch of important positions
    ## ie high coverage and high mod_freq difference
    ## but how to do it?
    if mod_diff<minModDiff: sig_mer = False
    else: sig_mer = True
    # estimate modification frequency on the reads that were not used for training
    pe = 0
    outdata = []
    # iterate over positions
    if not pos_data:
        return (mer, clf), outdata
    for pos_info in pos_data.split(";"):
        chrom, pos, strand, cov = pos_info.split()
        pos, cov = int(pos), int(cov)
        # for very low coverage sample it can happen that there is no samples after ~sel
        # this can be done better! maybe predict across all positions before and here just calculate mod_freq?
        # yes, we can predict for all and then change sel to NAN and calculate .mean() on for each position
        if len(alldata[0][pe:pe+cov][~sel[pe:pe+cov]])<1:
            pe += cov
            continue
        # subset all data 
        mod_freq = [clf.predict(alldata[si][pe:pe+cov][~sel[pe:pe+cov]]).mean() for si in range(2)]
        pe += cov
        # if mer is significant, take all positions for that mer having enough mod_freq difference
        # if mer isn't significant, take only positions having enogh mod_freq difference and high coverage
        if np.abs(np.diff(mod_freq))[0]<minModDiff or not sig_mer and cov<50: continue
        # get 21-mer sequence
        base, mer_seq = get_kmer(faidx, chrom, pos-1, strand, r2l, k=10)
        outdata.append((chrom, pos, base, strand, mer_seq, cov, *mod_freq))
    return (mer, clf), outdata

def is_qcfail(a, mapq=20):
    """Return True if read fails QC"""
    if a.mapq<mapq: return True

def get_coverage_for_ref(sam, ref, mapq, reflen): 
    """Return coverage from sam (pysam.AlignmentFile)"""
    coverage = np.zeros(reflen, dtype='uint16')
    for a in sam.fetch(ref):
        if is_qcfail(a, mapq): continue
        if a.blocks:
            for s, e in a.blocks: coverage[s:e] += 1
        else: coverage[a.pos:a.aend] += 1
    return coverage

def get_regions_with_mers(faidx, ref, s, e, mers, k=7):
    """Return positions that contain certain mers"""
    if not mers: return [(ref, s, e), ]
    regions = []
    seq = faidx.fetch(ref, s, e)
    for mi in range(len(seq)-k):
        mer = seq[mi:mi+k]
        if mer in mers:
            # store mer center positions
            regions.append((ref, s+mi+k//2, s+mi+k//2+1))
    return regions
    
def get_covered_regions(bams, fasta, step=1000, minCov=10, mapq=20,
                        maxdist=100, threads=6, mers=set()):
    """Return regions covered by at least minCov reads.
    
    Here, we only analyse the file with lower number of reads.
    """
    # iterate chromosomes
    sams = [pysam.AlignmentFile(bam, threads=threads) for bam in bams]
    # file with fewer alignemnts first
    sams = list(sorted(sams, key=lambda x: x.mapped))
    sam = sams[0]
    faidx = pysam.FastaFile(fasta)
    regions = []
    logger("Retrieving regions covered by %s+ reads..."%minCov)
    # first run idxstats and skip references with <10 reads altogether - usefull for transcript alignements
    ref2algs = {s.contig: s.mapped for s in sam.get_index_statistics()}
    ref2len = {r: l for r, l in zip(sam.references, sam.lengths)}
    bases = 0
    for ri, ref in enumerate(sam.references, 1):
        if ref not in faidx.references or ref not in ref2algs or ref2algs[ref]<minCov: continue
        sys.stderr.write(" %s / %s %s ...\r"%(ri, len(sam.references), ref))
        # get min coverage from bam file that has fewer mapped reads
        coverage = get_coverage_for_ref(sam, ref, mapq, ref2len[ref])
        # get regions with minCov coverage
        covered = np.where(coverage>=minCov)[0]
        for positions in get_consecutive(covered, maxdist):
            if len(positions)<1: continue
            s, e = positions[0]+1, positions[-1]+1
            bases += e-s
            # further split regions for max windows
            while s < e-step:
                regions += get_regions_with_mers(faidx, ref, s, s+step, mers)
                s += step
            regions += get_regions_with_mers(faidx, ref, s, e, mers)
    logger(" {:,} bases in {:,} regions to process.\n".format(bases, len(regions)))
    return regions

def get_consecutive(data, stepsize=1):
    """Return consecutive windows allowing given max. step size"""
    return np.split(data, np.where(np.diff(data) > stepsize)[0]+1)

def get_regions(fasta, step=10000):
    """Return windows of size step with sliding windows of step from FastA file"""
    faidx = pysam.FastaFile(fasta)
    regions = []
    for chrom, size in zip(faidx.references, faidx.lengths):
        for s in range(0, size, step):
            e = s+step if size>s+step else size
            regions.append((chrom, s, e))
    return regions

def get_significant_mers(df, Pval, KSfn, ignoreKSscores=False, logger=sys.stderr.write):
    """Return significant mers from KS test and save plots."""
    mers = set()
    bins = np.arange(0, 1, 1/50)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    cols = [c for c in df.columns if c.startswith("KS ")]
    for c, color in zip(cols, sns.color_palette(n_colors=len(cols))):
        # plot hist
        x = df[["mer",c]].melt("mer")["value"].to_numpy()
        n, _, _ = ax.hist(x, bins, label=c, color=color, alpha=0.5)
        # make sure that KS P-value is below P-value cut-off!
        sel = df[c.replace("KS", "P")]<Pval
        if not ignoreKSscores:
            # and that KS statistics that is different (P-value cut-off of 0.001)
            # from normal distribution of KS across all mers for given feature
            pval = 2.0 * scipy.stats.norm.cdf(-np.abs(df[c]-df[c].mean())/df[c].std())
            sel = sel&(pval<Pval)
            # add P-value cut-off to the histogram - two-sided, maybe should be one-sided?
            pvalx = scipy.stats.norm.ppf(1-Pval/2) * df[c].std() + df[c].mean()
            ax.vlines(pvalx, 0, max(n), label="%s P=%s"%(c, Pval), color=color, ls="--") 
        new_mers = df.loc[sel, "mer"]
        logger("  %s significant mers for %s"%(len(new_mers), c))
        mers.update(new_mers)
        # plot hist of significant mers
        ax2.hist(x[sel], bins, label=c, color=color, alpha=0.5)

    logger(" overall %s unique mers are different at P-value cut-off of %s"%(len(mers), Pval))
    if len(mers)==0:
        logger("There is no mers significantly different between provided samples.")
        logger("Please relax significance threashold (-p / --Pval) or try other samples...")
        sys.exit(0)
    # add plot vars
    name = KSfn.split(os.path.sep)[-2]
    ax.legend(); ax.set_xlabel("KS statistics"); ax.set_ylabel("Count"); ax.set_title("%s histogram of KS for %s features"%(name, len(cols)))
    ax2.legend(); ax2.set_xlabel("KS statistics"); ax2.set_ylabel("Count"); ax2.set_title("%s histogram: %s significant mers"%(name, len(mers)))
    if os.path.isdir(os.path.dirname(KSfn)):
        for ext in ("pdf", "png"): fig.savefig("%s.hist.%s"%(KSfn, ext))
    return mers

def check_bam_index(bams, logger=sys.stderr.write):
    """Warn if index file is older than BAM. Exit if no index exists."""
    for bam in bams:
        # always use csi if exists
        idx = bam+".csi" if os.path.isfile(bam+".csi") else bam+".bai"
        if not os.path.isfile(idx):
            logger("[WARNING] No index for %s  \n"%bam)
            sys.exit(1)
        # warn if BAM index is older than BAM
        diff = os.path.getmtime(bam)-os.path.getmtime(idx)
        if diff>1:
            logger("[WARNING] Index is %s seconds older than BAM file for %s   \n"%(int(diff), bam))

def mod_report(outdir, bam, fasta, threads, regionsfn, features, mapq=15,
               minReads=10, maxReads=100, minModFreq=0.1, nn=1, rna=True, 
               Pval=0.001, ignoreKSscores=False, logger=sys.stderr.write,
               window=1000, timeout=60*60):
    """Report likely modified positions from BAM files"""
    if len(bam)!=2:
        logger("[ERROR] Provide 2 BAM files!\n")
        sys.exit(1)
    # make sure index files exist
    check_bam_index(bam, logger)
    # skip if outfile exists and isn't empty (gz is 141 even if empty)
    outfn = os.path.join(outdir, "denovo.tsv.gz")
    if os.path.isfile(outfn) and os.stat(outfn).st_size>141:
        logger(" %s exists"%outfn)
        return outfn
    # start pool of processes
    p = Pool(threads, initializer=init_args, initargs=(bam, fasta))
    # get pool of workers 
    pp = ProcessPool(max_workers=threads) #, max_tasks=10)    
    # prepare output directory
    if not os.path.isdir(outdir): os.makedirs(outdir)
    # get regions to process
    if regionsfn: regions = load_bed(regionsfn)
    elif rna: regions = get_covered_regions(bam, fasta, window, minReads, mapq, threads=threads)
    else: regions = get_regions(fasta, window)
    if not regions:
        logger("No regions to process. Consider changing -d / --minDepth.")
        sys.exit(1)
    # get KS statistics for 7-mers
    KSfn = os.path.join(outdir, 'KS.tsv.gz')
    if not os.path.isfile(KSfn) or os.stat(KSfn).st_size<=141:
        # first get significant mers using only SI0 and TR0
        features0 = ["si", "tr"] if "tA" in features or "tr" in features else features 
        logger("Loading per-mer data for %s ..."%" & ".join(features0))
        feature_names = ["%s_%s"%(f.upper(), i) for f in features0 for i in range(1)]
        if threads>1: mer2data = load_data(fasta, bam, features0, 0, rna, regions, minReads, threads, maxReads=maxReads, mapq=mapq, p=p)
        else: mer2data = load_data_single(fasta, bam, features0, 0, rna, regions, minReads, maxReads=maxReads, mapq=mapq)

        logger("Reporting per-mer KS statistics to %s ..."%KSfn)
        ks = get_KS(mer2data, threads, p=p)
        # note above with pop all elements from mer2data!
        df = pd.DataFrame(ks, columns=["mer",]+["%s %s"%(s, f) for f in feature_names for s in ("KS", "P")])
        df.to_csv(KSfn, sep="\t", index=False)
    else:
        df = pd.read_csv(KSfn, sep="\t")
        
    # get significant mers - here relying on P-value isn't enough
    mers = get_significant_mers(df, Pval, KSfn, ignoreKSscores, logger)
    
    # later get all features for subset of mers
    feature_names = ["%s_%s"%(f.upper(), i) for f in features for i in range(-nn, nn+1)]
    logger("Loading %s features for %s mers..."%(len(feature_names), len(mers)))
    if threads>1: mer2data = load_data(fasta, bam, features, nn, rna, regions, minReads, threads, mers, maxReads=maxReads, mapq=mapq, p=p)
    else: mer2data = load_data_single(fasta, bam, features, nn, rna, regions, minReads, mers, maxReads=maxReads, mapq=mapq)

    #def report(): move below to separate function
    #"""Report modiifed baess"""
    logger("Reporting per-position modification frequency to %s ..."%outfn)
    # write header
    out = gzip.open(outfn, "wt")
    tsv = csv.writer(out, delimiter='\t')
    header = ["chrom", "pos", "base", "strand", "mer", "min coverage"] + ["%s mod_freq"%b for b in bam]
    tsv.writerow(header)
    
    t0 = datetime.now()
    freqs, skipped, mer2clf = [], [], {}
    posN = 0
    # get list of mers to process
    moi = [mer for mer in sorted(mers) if mer in mer2data]
    args = ((fasta, mer, mer2data[mer], feature_names, minModFreq) for mer in moi)
    # process mers with timeout, because rarely it gets stucked
    future = pp.map(get_calls, args, timeout=timeout)
    for i, ((mer, clf), data) in enumerate(get_results_from_pool(future, moi, timeout), 1):
        sys.stderr.write(" {:,} / {:,} {} \r".format(i, len(mers), datetime.now()-t0))
        if data:
            tsv.writerows(data)
            mer2clf[mer] = clf
        else:
            skipped.append(mer)
            continue
        posN += len(data)
        freqs.append((mer, *clf.freqs, *clf.feature_importances_)) #.steps[-1][-1]
    # close pool of workers; skip join as it takes lots of memory
    #pp.close()
    out.close()
    logger(" {:,} positions processed".format(posN))
    # and store clf
    modelsfn = os.path.join(outdir, "models.lzma")
    logger("Saving ML models to %s ..."%modelsfn)
    joblib.dump(mer2clf, modelsfn)
    # save features
    cols = ["mer", "KO freq", "WT freq"]+["%s importance"%f for f in feature_names]
    df1 = pd.DataFrame(freqs, columns=cols)
    fifn = os.path.join(outdir, 'feature_importances.tsv.gz')
    df1.to_csv(fifn, sep="\t", index=False)
    logger("Per-mer estimated modification frequency and feature importances saved to %s"%fifn)
    return outfn

def get_highest_from_consecutive(df2, minModFreq):
    """Return the positions with the highest modification frequency difference 
    from consecutive positions. 
    """
    # get difference between KO and WT
    df2["diff"] = df2[df2.columns[-2:]].diff(axis=1)[df2.columns[-1]]#.abs()
    # filter out those with small difference
    df2.drop(index=df2[df2["diff"]<minModFreq].index, inplace=True)
    # this needs to be strand-specific!!!
    # get consecutive positions
    df2["group"] = (df2["pos"].diff(1)!=1).astype('int').cumsum()
    # and keep only positions with max difference across each group of consecutive positions
    idx = df2.groupby(["group"])["diff"].transform(max) == df2["diff"]
    return df2[idx].reset_index()

def get_motifs(outfn, fasta, threads, minModFreq, nmotifs=5, Evalue=1e-05, logger=sys.stderr.write):
    """Run meme on significant motifs"""
    outdir = os.path.dirname(outfn)
    memedir = os.path.join(outdir, "meme")
    motiffn = os.path.join(outdir, 'denovo.motifs.tsv.gz')
    if not os.path.isdir(memedir): 
        # build genome model if not exist
        if not os.path.isfile("%s.ooc"%fasta):
            logger("Generating markov model...")
            os.system("fasta-get-markov %s %s.ooc"%(fasta, fasta))
        logger("Loading candidate positions...")
        # read reported modified positions
        df2 = pd.read_csv(outfn, sep="\t").sort_values(["chrom", "strand", "pos"])
        logger(" {:,} positions loaded".format(len(df2)))
        # store only positions with highest modification difference between KO & WT
        # from each group of consecutive positions
        df2 = get_highest_from_consecutive(df2, minModFreq)
        logger("  {:,} significant positions with max modification difference per group".format(len(df2)))
        if len(df2)<2: # meme needs at least 2 sequences to compute motifs
            logger("Nothing else to do.")
            return memedir
        logger("Running MEME...")
        motiffa = os.path.join(outdir, "motifs.fa")
        with open(motiffa, "wt") as out:
            # this sometimes causes out of memory error!
            #out.write(df2.to_string(columns=["chrom", "pos", "strand", "mer"], index=False, header=False, 
            #                        formatters=[">{}".format, ":{}".format, "{}".format, "\n{}".format]).replace(" ", ""))
            for idx, r in df2[["chrom", "pos", "strand", "mer"]].iterrows(): out.write(">%s:%s%s\n%s\n"%tuple(r))
        cmd = "meme -nostatus -brief 1000000000 -p %s -dna -mod zoops -nmotifs %s -minw 2 -maxw 10 -bfile %s.ooc -o %s %s"%(threads, nmotifs, fasta, memedir, motiffa)
        os.system(cmd)
        # create symlink
        #os.symlink("%s/meme.html %s"%(memedir, outdir))
        # parse meme xml
        tree = ET.parse(os.path.join(memedir, "meme.xml"))
        root = tree.getroot() #['training_set', 'model', 'motifs', 'scanned_sites_summary']
        # get significant motifs
        sig_motifs = [m for m in root[2] if float(m.attrib["e_value"])<Evalue]
        if len(sig_motifs)<1:
            logger("No significant motifs were found!")
            return memedir
        logger("Annotating positions with %s significant motif(s)..."%len(sig_motifs))
        for m in sig_motifs:
            logger(" {} {} with {} sites and E-value {}".format(*[m.attrib[k] for k in ["id", "name", 'sites', "e_value"]]))
            motif_name = " ".join([m.attrib[k] for k in ["id", "name"]]) + " P-value"
            # get idx of sequences with given motif - sequence_id starts from sequence_0
            idx = [int(e.attrib["sequence_id"].split("_")[-1]) for e in m[-1]]
            pvals = np.array([float(e.attrib["pvalue"].split("_")[-1]) for e in m[-1]])
            # store motifs p_values
            df2[motif_name] = 1
            df2.loc[idx, motif_name] = pvals
        # finally filter-out positions without motif and store file with annotated motifs
        motif_cols = [c for c in df2.columns if c.startswith("motif_")]
        sel = np.any([df2[c]!=1 for c in motif_cols], axis=0)
        df2 = df2[sel]
        for c in motif_cols:
            m = c.split()[1]
            df2.loc[df2[c]!=1, "motif"] = m
        # store before renaming columns
        df2.to_csv(motiffn, sep="\t", index=False)
        logger("{:,} positions with motifs reported to {}".format(len(df2), motiffn))
    else:
        logger(" %s exists"%memedir)
        if not os.path.isfile(motiffn): return memedir
        # just load final list of positions
        df2 = pd.read_csv(motiffn, sep="\t", low_memory=False)
    # rename columns
    df2.rename(inplace=True, columns={c: n for n, c in zip(("ko", "wt"), [c for c in df2.columns if c.endswith("mod_freq")])})
    # melt and plot catplot
    df2melt = df2[["motif", "strand", "ko", "wt"]].melt(["motif", "strand",], var_name="sample", value_name="mod_freq"); 
    g = sns.catplot(x="sample", y="mod_freq", hue="strand", col="motif", kind="box", data=df2melt) #kind="violin"
    fig = g.fig; fig.suptitle("Per-motif modification frequency")
    fig.savefig(os.path.join(outdir, "per_motif.mod_freq.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(outdir, "per_motif.mod_freq.png"), bbox_inches="tight")
    return memedir

def store_bed(outdir):
    """Save positions as bed file"""
    def get_motif(r): return r.motif
    def get_base(r): return "mod%s"%r.base

    nfiles = 0
    for fn in ("denovo.tsv.gz", "denovo.motifs.tsv.gz"):
        fn = os.path.join(outdir, fn)
        outfn = fn[:-7]+".bed" # denovo.tsv.gz > denovo.bed
        if os.path.isfile(outfn):
            logger(" %s exists"%outfn)
            nfiles += 1
            continue
        if not os.path.isfile(fn): return nfiles
        nfiles += 1
        df = pd.read_csv(fn, sep="\t")
        df.sort_values(list(df.columns[:2]), inplace=True)
        if "diff" not in df:
            df["diff"] = (df[df.columns[-1]]-df[df.columns[-2]]).abs()    
        bi = 0
        bedMline = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
        with open(outfn, "w") as out:
            if "motif" in df: get_name = get_motif
            else: get_name = get_base
            for idx, r in df.iterrows():
                s, e, cov, d = r.pos-1, r.pos, r["min coverage"], r["diff"]
                out.write(bedMline%(r.chrom, s, e, get_name(r), d, r.strand, s, e, 
                                    "%s,0,0"%int(round(d*255)), cov, int(round(d*100))))
                bi += 1
        logger(" %s mods saved in %s"%(bi, outfn))
    return nfiles

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-i", "--indirs", nargs=2, help="input directories for PCR/IVT/KO and native sample")
    parser.add_argument("-r", "--recursive", action='store_true', help="recursive processing of input directories [%(default)s]")
    parser.add_argument("--rna", action='store_true', help="project is RNA sequencing [DNA]")
    parser.add_argument("--sensitive", action='store_true', help="use sensitive mapping parameters ie tRNA")
    parser.add_argument("-o", "--outdir", required=True, help="output directory")
    parser.add_argument("-f", "--fasta", required=1, help="reference FASTA file")
    parser.add_argument("-t", "--threads", default=6, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("-e", "--encode", action='store_true', help="encode modifications in BAM")
    parser.add_argument("--mapq", default=15, type=int, help="min mapping quality [%(default)s]")
    parser.add_argument("-d", "--minDepth", default=10, type=int, help="min depth of coverage (for each strand separately) [%(default)s]")
    parser.add_argument("--maxDepth", default=100, type=int, help="max depth of coverage (for each strand separately) [%(default)s]")
    parser.add_argument("--minModFreq", default=0.01, type=float, help="min modification frequency difference [%(default)s]")
    parser.add_argument("--Pval", default=0.01, type=float, help="P-value cut-off for KS [%(default)s]") #& motifs 
    parser.add_argument("--ignoreKSscores", action="store_true", help="ignore KS scores [use KS scores]")
    parser.add_argument("-b", "--bed", help="BED file with regions to analyse [optionally]")
    parser.add_argument("-n", "--nn", default=1, type=int, help="neighbours to consider [%(default)s]")
    parser.add_argument("--features", default=["si", "mp", "dt0", "dt10", "tA", "tC", "tG", "tT"], nargs="+",
                        help="features to use [%(default)s]")
    
    guppy = parser.add_argument_group("Basecalling options") #mutually_exclusive
    guppy.add_argument("-c", "--config", default="dna_r9.4.1_450bps_pcr_hac.cfg",
                       help="guppy modification-unaware model [%(default)s]")
    guppy.add_argument("--host", "--guppy_basecall_server", default="",
                        help="guppy server hostname or path to guppy_basecall_server binary [use basecall information from Fast5]")
    guppy.add_argument("--port", default=5555, type=int,
                        help="guppy server port (this is ignored if binary is provided) [%(default)s]")
    guppy.add_argument("--device", default="cuda:0", help="CUDA device to use [%(default)s]")
    guppy.add_argument("--timeout", default=60*60, type=int, help="timeout in seconds to process each Fast5 file [%(default)s]")
    
    meme = parser.add_argument_group("Motif-enrichment options")
    meme.add_argument("--evalue", default=1e-05, type=float, help="E-value cut-off [%(default)s]")
    meme.add_argument("--nmotifs", default=5, type=int, help="number of motifs to report [%(default)s]")

    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))
    
    if o.indirs[0]==o.indirs[1]:
        logger("Two identical files provided")
        sys.exit(1)
        
    bams = get_features(o.indirs, o.fasta, o.threads, o.rna, o.sensitive, 
                        o.config, o.host, o.port, o.recursive, o.device, o.timeout)
    
    # process everything only if outfile does not exists
    logger("Retrieving likely modified motifs & positions...")
    outfn = mod_report(o.outdir, bams, o.fasta, o.threads, o.bed, o.features, o.mapq,
                       o.minDepth, o.maxDepth, o.minModFreq, o.nn, o.rna, o.Pval, o.ignoreKSscores,
                       logger=logger, timeout=o.timeout)

    # calculate MEME motif enrichment - meme v5+ (conda install -U meme)
    logger("Calculating enriched motifs...")
    memedir = get_motifs(outfn, o.fasta, o.threads, o.minModFreq, o.nmotifs, o.evalue, logger)

    logger("Reporting BED files...")
    if store_bed(o.outdir):

        # encode modifications in bam
        if o.encode:
            edir = os.path.join(o.outdir, "encode")
            if not os.path.isdir(edir):
                import encode_mods as em
                models = [os.path.join(o.outdir, "models.lzma")]
                em.mod_encode(edir, o.indirs, o.fasta, models, o.threads, o.rna, o.mapq)#, "", o.config, o.host, o.port, o.device)
            else:
                logger(" %s exists"%edir)
    logger("Done!")
    
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
