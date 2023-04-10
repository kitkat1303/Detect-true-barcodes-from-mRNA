from Bio.Seq import Seq
import os
import math
import gzip
import numpy as np
import pandas as pd
import os, errno
import multiprocessing as mp


fileName = "10_5_candidates.fastq.gz"
#data = fastqToDF(fileName)
#findAnchorSeq(data)


#convert file into data frame
def fastqToDF(filename, size=100000):
    """Convert fastq to dataframe.
        size: limit to the first reads of total size
        Returns: dataframe with reads
    """

    ext = os.path.splitext(filename)[1]
    if ext=='.fastq' or ext=='.gz':
        fastq_parser = SeqIO.parse(gzopen(filename, "rt"), "fastq")
    else:
        fastq_parser = SeqIO.parse(open(filename, "r"), "fastq")
    i=0
    res=[]
    for fastq_rec in fastq_parser:
        #print (fastq_rec.seq)
        i+=1
        if i>size:
            break
        res.append([fastq_rec.id, str(fastq_rec.seq)])
    df = pd.DataFrame(res, columns=['id','seq'])
    df['length'] = df.seq.str.len()
    return df

#find anchor sequence and return a new data frame
#def findAnchorSeq(df):
  
  #find anchor seq and if found at to previous 30 characters to data frame
  # return df