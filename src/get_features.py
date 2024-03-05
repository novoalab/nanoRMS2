#!/usr/bin/env python3
desc="""Basecall Fast5 (or read basecalled Fast5), align & 
requiggle reads and store features in BAM file. 

For all reference bases we store (as BAM comments):
- normalised signal intensity mean [tag si:B:f]
- base probability [tag tr:B:C] retrieved from guppy (trace scaled 0-255)
- dwell time [tag dt:B:C] in signal step capped at 255
- probability of modification calculated by tombo [tag mp:B:f]

Starting from v0.11a all features are matched versus padded reference sequnce blocks 
ie excluding introns and large (padded) deletions from reference. 
Those blocks (2-D array of start & ends) are stored as flattened 1-D array [tag bs:B:i]
ie. exons [(8114, 8244), (8645, 8797)] will be stored as array('I', [8114, 8244, 8645, 8797]). 

Also, --rna will: 
- automatically enable spliced alignments 
- and estimates of polyA length [tag a1, a2, a3, a4] using mean, median, median_A and find_peaks

TBD:
- clean tmp .bam files (only if all finished properly)
- combine logs
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Cologne/Barcelona/Mizerów/Cap Salou, 18/02/2021
"""

import itertools, json, os, sys, pysam, traceback
import mappy
import numpy as np
from datetime import datetime
from pathlib import Path
from pebble import ProcessPool, ProcessExpired
from concurrent.futures import TimeoutError
from array import array
from basecall import init_args, start_guppy_server, get_basecall_client, basecall_and_align
from resquiggle import resquiggle_reads
from common import VERSION, logger, warning

# only DNA bases as in SAM U is always referred as T
bases = "ACGT"
base2idx = {b: i for i, b in enumerate(bases)}
#base2idx["N"] = base2idx["A"] # dirty N>A replacement for tRNA
base2complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"} 
# add lower-case for get_aligned_pairs as it reports substitutions as lower-case
for b, i in list(base2idx.items()): base2idx[b.lower()] = i
for b, c in list(base2complement.items()): base2complement[b.lower()] = c
            
def get_trace_for_reference_bases(a, move, trace, rna):
    """Return reference-aligned trace for tr (ref base), tA, tC, tG, tT"""
    def get_bidx_fwd(b): return base2idx[b]
    def get_bidx_rev(b): return base2idx[base2complement[b]]
    # trace for tA/C/G/T, insertion and tr (reference-base)
    trace = trace.astype("float64") # working on float64 is a bit faster
    atrace = np.zeros((5, a.reference_length))#, dtype="uint8")
    # trace and move data from read
    move_pos = np.append(np.argwhere(move==1).flatten(), len(trace)) # add end of trace
    # here we need to remember that DNA 5'>3', but RNA 3'>5'
    # plus the strand matters
    if a.is_reverse: # for REV alg
        get_bidx = get_bidx_rev # take complement base
        if not rna: move_pos = move_pos[::-1] # reverse move_pos for DNA
    else: # for FWD alg
        get_bidx = get_bidx_fwd # take base
        if rna: move_pos = move_pos[::-1]     # reverse move_pos for RNA
    # fucking cols in atrace has to respond to rev alignments as well!
    ## but we really do want to get complement here as it'll ease feature extraction ;)
    # process aligned bases - that's quite elegant, right? :P
    ## with_seq require MD tags: in minimap2 use --MD and -Y (soft-clip supplementary)
    ref_base_ai = np.zeros(atrace.shape[1], dtype="int")
    sums = np.ones(a.reference_length)
    for qi, ri, b in a.get_aligned_pairs(matches_only=True, with_seq=True):
        if b=="N": continue
        s, e = move_pos[qi], move_pos[qi+1]
        if s>e: e, s = s, e
        ai = ri-a.reference_start   # get ref position
        # tA/C/G/T: probabilities for A, C, G & T/U
        atrace[:4, ai] = trace[s:e].sum(axis=0)
        sums[ai] = e-s
        # and store baseidx - this isn't even needed as it can be computed by bam2data
        ref_base_ai[ai] = get_bidx(b)
    # tr: probability for the reference base
    atrace[-1] = atrace[ref_base_ai, np.arange(atrace.shape[1])]
    # this way we save ~25% time compared to .sum()/(e-s) and even more compare to mean
    atrace /= sums
    return atrace.astype("uint8")

def get_sam_from_resquiggle(fast5, ref, rna=False, conf=('', '', ''), return_res=False,
                            add_si=True, max_per_ref=0, only_forward=False, verbose=True):
    """Yield alignment (pysam.AlignmentSegment) with features encoded as tags
    for reads from Fast5 file. 

    If fast5 is a directory, process all Fast5 within that directory. 

    If return_res=True, it'll also yield resquiggle object from tombo.
    """
    # get all fast5 file from within directory
    fnames = sorted(map(str, Path(fast5).rglob('*.fast5'))) if os.path.isdir(fast5) else [fast5, ]
    i = 0
    errors = {}
    mssg = "  %s - %s reads skipped: %s \r"
    for fn in fnames:
        sam = basecall_and_align(fn, ref, rna, conf, only_forward=only_forward)
        for i, (res, err) in enumerate(resquiggle_reads(sam, ref, rna, add_si=add_si,
                                                        max_per_ref=max_per_ref), i+1):
            if verbose and not i%100:
                sys.stderr.write(mssg%(i, sum(errors.values()), str(errors)))
            if err:
                if err not in errors: errors[err] = 1
                else: errors[err] += 1
                continue
            # get pysam alignment object & exonic blocks
            a, blocks = res.a, res.align_info.blocks
            # catch problems - here exonic seq will have different length
            #if len(si)!=sum([e-s for s, e in blocks]): continue
            # get dwell times capped at 255 to fit uint8 (1 byte per base)
            si, mp = res.si, res.mp
            dt = res.segs[1:]-res.segs[:-1]
            dt[dt>255] = 255
            # get reference-aligned base probabilities: tA, tC, tG, tT & tr (ref base)
            tr = get_trace_for_reference_bases(a, res.move, res.trace, rna) # this takes 189µs (>50%) of time!
            if a.is_reverse: si, dt = si[::-1], dt[::-1]
            # and finally set tags matching refseq
            ## but if alignment reaches seq end the end signal/probs will be wrong!
            ## same at exon-intron boundaries
            a.set_tag("bs", array("i", np.array(blocks).flatten()))
            a.set_tag("si", array("f", si))
            a.set_tag("mp", array("f", mp))
            a.set_tag("dt", array("B", dt))
            # tA/C/G/T are strand oriented, but tr correspond to reference base
            # get exonic tr
            exonic_pos = np.concatenate([np.arange(s, e) for s, e in blocks])
            tr = tr[:, exonic_pos-a.pos]
            for ci, tag in enumerate(("tA", "tC", "tG", "tT", "tr")): #"ti",
                a.set_tag(tag, array("B", tr[ci]))
            yield (a, res) if return_res else a

def store_features(args):
    """Process individual Fast5 files"""
    ref, fn, conf, rna = args
    outfn = "%s.bam"%fn
    if os.path.isfile(outfn): return outfn
    # use some tmpname
    tmpoutfn = outfn+"_"
    faidx = pysam.FastaFile(ref)
    # disable pysam verbosity https://github.com/pysam-developers/pysam/issues/939#issuecomment-669016051
    pysam.set_verbosity(pysam.set_verbosity(0))
    # get output bam
    bam = pysam.AlignmentFile(tmpoutfn, "wb", reference_names=faidx.references, reference_lengths=faidx.lengths)    
    # store read alignment with additional info
    ## instead we could just store them in list and save sorted to reduce IO at the cost of memory
    ## the code for that is in modPhred branch
    for a in get_sam_from_resquiggle(fn, ref, rna, conf):
        bam.write(a)
    # close, sort & clean-up
    bam.close() # this is not sorted
    pysam.sort("-o", outfn, tmpoutfn)
    os.unlink(tmpoutfn)
    # write error report
    #with open('%s.json'%outfn, 'w') as f:
    #    errors["Alignements"] = i # store number of alignements
    #    f.write(json.dumps(errors))
    return outfn

def get_results_from_pool(future, args, timeout=0):
    """Generator of output from pool of workers catching errors and time-outs."""
    i = 0
    iterator = future.result()
    while True:
        try:
            yield next(iterator)
            sys.stderr.write(" %s / %s \r"%(i+1, len(args)))
        except StopIteration:
            break
        except TimeoutError as error:
            #print(error.args, args) # error.args[1] worked before instead timeout, but now it causes an error
            logger("[WARNING] worker took longer than %d seconds for %s "%(timeout, args[i]))
        except ProcessExpired as error:
            logger("%s. Exit code: %d for %s "%(error, error.exitcode, args[i]))
        except Exception as error:
            logger("worker raised %s for %s "%(error, args[i]))
            # this causes AttributeError: 'TypeError' object has no attribute 'traceback' in Python v3.7.9
            logger(traceback.format_exc())#error.traceback)  # Python's traceback of remote process
        i += 1

def get_features(indirs, fasta, threads=1, rna=True, sensitive=False, 
                 config="dna_r9.4.1_450bps_pcr_fast", host="localhost", port=5555,
                 recursive=False, device="cuda:0", timeout=60*30):
    """Process multiple directories from Fast5 files"""
    logger("Retrieving features from Fast5 files in %s directories...\n"%len(indirs))    
    # skip feature storing if BAM files exists
    bams = []
    for indir in indirs:
        bam = indir.rstrip(os.path.sep) + ".bam"
        if os.path.isfile(bam): bams.append(bam)
    if len(bams)==len(indirs):
        logger(" BAM files already exists: %s !"%", ".join(bams))
        return bams
        
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
    kwargs = {"preset": "spliced" if rna else "map-ont", "k": 13}
    if sensitive:
        kwargs.update({"k": 6, "w": 3, "scoring": [1, 1, 1, 1, 32, 0], 
                       "min_cnt": 1, "min_chain_score": 13, "min_dp_score": 20})
    aligner = mappy.Aligner(fasta, **kwargs)
    # start pool of workers
    # it's important to initialise the pool with aligner object as it can't be pickled
    p = ProcessPool(max_workers=threads, initializer=init_args, initargs=(aligner, )) #, max_tasks=10)
    # process directories
    bams = []
    for indir in indirs:
        # skip if out bam exists
        bam = indir.rstrip(os.path.sep) + ".bam"
        bams.append(bam)
        if os.path.isfile(bam):
            logger(" %s exists"%bam)
            # pysam has issues with .csi index https://github.com/pysam-developers/pysam/issues/1033
            if os.path.isfile(bam+".csi"):
                os.rename(bam+".csi", bam+".csi0")
                pysam.index("-@ %s"%threads, bam)
            continue
        # make sure the indir exists
        if not os.path.isdir(indir):
            sys.stderr.write("[ERROR] No such directory: %s\n"%indir)
            sys.exit(1)
        # get Fast5 files
        if recursive:
            fnames = sorted(map(str, Path(indir).rglob('*.fast5')))
        else:
            fnames = sorted(map(str, Path(indir).glob('*.fast5')))
        logger(" %s with %s Fast5 file(s)...\n"%(indir, len(fnames)))
        # basecall individual fast5 on-the-fly
        args = [(fasta, fn, conf, rna) for fn in fnames]
        future = p.map(store_features, args, timeout=timeout)
        _bams = [b for b in get_results_from_pool(future, fnames, timeout)]
        # merge & index BAM files
        if _bams:
            logger("  merging %s BAM files > %s ...\n"%(len(_bams), bam))
            pysam.merge("-@ %s"%threads, "--write-index", bam, *_bams)
            # clean-up tmp fast5.bam files
        else:
            logger("[WARNING] No BAM files to merge!")
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
    parser.add_argument("-i", "--indirs", nargs="+", help="input directory with Fast5 files")
    parser.add_argument("-r", "--recursive", action='store_true', help="recursive processing of input directories [%(default)s]")
    parser.add_argument("--rna", action='store_true', help="project is RNA sequencing [DNA]")
    parser.add_argument("--sensitive", action='store_true', help="use sensitive mapping parameters ie tRNA")
    parser.add_argument("-f", "--fasta", required=1, help="reference FASTA file")
    parser.add_argument("-t", "--threads", default=6, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("-c", "--config", default="dna_r9.4.1_450bps_pcr_fast.cfg", help="guppy model [%(default)s]")
    parser.add_argument("--host", "--guppy_basecall_server", default="",
                        help="guppy server hostname or path to guppy_basecall_server binary [use basecall information from Fast5]")
    parser.add_argument("-p", "--port", default=5555, type=int,
                        help="guppy server port (this is ignored if binary is provided) [%(default)s]")
    parser.add_argument("-d", "--device", default="cuda:0", help="CUDA device to use [%(default)s]")
    parser.add_argument("--timeout", default=60*60, type=int, help="timeout in seconds to process each Fast5 file [%(default)s]")
    # 30 minute time-out was causing problems with dna damage
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    get_features(o.indirs, o.fasta, o.threads, o.rna, o.sensitive, 
                 o.config, o.host, o.port, o.recursive, o.device, o.timeout)
        
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

