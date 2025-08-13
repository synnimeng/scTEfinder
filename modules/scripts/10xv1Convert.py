#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# # # # # # # # # # # # 
"""
╭───────────────────────────────────────╮ 
│ 10xv1Convert.py   2025/03/02/-19:45 
╰───────────────────────────────────────╯ 
│ Description:
    Transfer 10x-v1 pattern like [R1: [Barcode]{repeat}, R2: [Reads][UMI]{repeat}] for STAR inputs
""" # [By: Synni Meng & yiwen Hu]

# region |- Import -|
from time import time
import concurrent.futures as futures
import argparse
import gzip
import os
# endregion


# region |- Util Function -|
def getRecord(fileHandler):
    try:
        seqID, seqInfo = next(fileHandler).rstrip().split(' ')
        seq = next(fileHandler).rstrip()
        symbol = next(fileHandler).rstrip()
        qual = next(fileHandler).rstrip()
        return seqID, seqInfo, seq, symbol, qual
    except StopIteration:
        return None

def concatRecord(seq1, seq2):
    # seqID, seqInfo, seq, symbol, qual
    return seq1[0], seq1[1], seq1[2]+seq2[2], seq1[3], seq1[4]+seq2[4] 

def saveRecord(seqID, seqInfo, seq, symbol, qual, fileHandler):
    return fileHandler.write(f'{seqID} {seqInfo}\n{seq}\n{symbol}\n{qual}\n')
# endregion


# region |- Handle Function -|
def runChunk(r1, r2, o1, o2, start, n_record, f=4, jobname='Unnamed Job', pp_gap=25000):
    # region |- input & output Reads Handler -|
    bc_Handler = gzip.open(r1, 'rt')
    read_Handler = gzip.open(r2, 'rt')
    R1_Handler = gzip.open(o1, 'wt')
    R2_Handler = gzip.open(o2, 'wt')
    # endregion
    # Skip to the start line
    for _ in range(start*f):
        next(bc_Handler)
        next(read_Handler)
        next(read_Handler)
    # Handle
    stampT = time()
    # Run main
    i = 0
    bc = getRecord(bc_Handler)
    while i < n_record and bc != None:
        if i % pp_gap == 0:
            print(f'\r{jobname}: {i}/{n_record} seqs processed, {i/(time()-stampT):.2f}/s', end="")
        # seqID:0, seqInfo:1, seq:2, symbol:3, qual:4
        read = getRecord(read_Handler)
        umi = getRecord(read_Handler)
        # check seq ID the same
        alignID = (bc[0] == read[0] == umi[0])
        if not alignID:
            print(f'{jobname}: Error In seq [{start+i}] ID aligned: \nBC:\t{bc.id}\nR:\t{read.id}\nUMI:\t{umi.id}')
            continue
        # concat bc + umi > out R1.fq.gz
        bc_umi = concatRecord(bc, umi)
        _=saveRecord(*bc_umi, R1_Handler)
        _=saveRecord(*read, R2_Handler)
        # next iter
        i += 1
        bc = getRecord(bc_Handler)
    # result
    stampT = time() - stampT
    print(f'{jobname}: Done, Elapsed Time: {stampT/60:.2f} min')
    # region |- Clean & Close gzip file -|
    bc_Handler.close()
    read_Handler.close()
    R1_Handler.close()
    R2_Handler.close()
    # endregion
    return stampT

# endregion



if __name__ == '__main__':
    exam = [
        '-r1', 'ACGCGGAA_L001_R1.fastq.gz',
        '-r2', 'ACGCGGAA_L001_R2.fastq.gz',
        '-o1', 'Demo_R1.fastq.gz',
        '-o2', 'Demo_R2.fastq.gz',
    ]
    # argparse
    parser = argparse.ArgumentParser(description='Convert 10x v1 FASTQ files to STAR inputs.')
    parser.add_argument('-r1', '--Read1', required=True, type=str, help='Input R1 file (barcode only)')
    parser.add_argument('-r2', '--Read2', required=True, type=str, help='Input R2 file ([Reads][UMI])')
    parser.add_argument('-o1', '--OutRead1', required=True, type=str, help='Output R1 file (barcode+umi)')
    parser.add_argument('-o2', '--OutRead2', required=True, type=str, help='Output R2 file (reads)')
    parser.add_argument('-n', '--nThread', default=8, required=False, type=int, help='n Threads cpu')
    args = parser.parse_args(  )
    # get job name
    jobR1 = args.OutRead1[::-1].split('.', maxsplit=1)[-1][::-1]
    jobR2 = args.OutRead2[::-1].split('.', maxsplit=1)[-1][::-1]
    # counts in 30s
    nRecord = sum(True for _ in gzip.open(args.Read1, 'rt'))//4
    gapRecord = 100000   # Interval roughly corresponding to processing capacity
    n_gap = round(args.nThread/2*(args.nThread+1))
    avgRecord = (nRecord-gapRecord*n_gap)//args.nThread
    job_pp = 0
    jobRegions = []
    for i in range(args.nThread):
        add_njob = avgRecord+(args.nThread-i)*gapRecord
        jobRegions.append(
            (job_pp, add_njob)
        )
        job_pp += add_njob
    #
    jobRegions[-1] = (jobRegions[-1][0], nRecord-jobRegions[-1][0]+gapRecord)
    print(f"Total {nRecord}:\n{jobRegions}")
    # timer
    stampT = time()
    # submit job 
    with futures.ThreadPoolExecutor(max_workers=args.nThread) as executor:
        jobs = [executor.submit(runChunk, 
                                r1=args.Read1, r2=args.Read2, 
                                o1=f"{jobR1}.tmp{i}.gz", o2=f"{jobR2}.tmp{i}.gz", 
                                start=P[0], n_record=P[1], 
                                jobname=f"Chunk {i}") 
                                for i, P in enumerate(jobRegions)]
        # Wait for all tasks to complete and collect results
        results = []
        for job in futures.as_completed(jobs):
            try:
                result = job.result()  # obtain task result
                results.append(result)
                print(f"Result: {result}")
            except Exception as e:
                print(f"Task failed with exception: {e}")
    print(f"All tasks completed. Mean Elapsed time: {sum(results)/len(results):.2f}")
    # merge all!
    print("Merge ...")
    tmpR1 = [f"{jobR1}.tmp{i}.gz" for i in range(args.nThread)]
    os.system(f"cat {' '.join(tmpR1)} > {args.OutRead1}")
    tmpR2 = [f"{jobR2}.tmp{i}.gz" for i in range(args.nThread)]
    os.system(f"cat {' '.join(tmpR2)} > {args.OutRead2}")
    # remove tmps
    os.system(f"rm {' '.join(tmpR1)}")
    os.system(f"rm {' '.join(tmpR2)}")
    # result
    stampT = time() - stampT
    print(f"Done All, Elapsed Time: {stampT/60:.2f} min")
