#!/usr/bin/env python
import os
import time
import random
import numpy as np
import multiprocessing as mp
from datetime import datetime
from collections import OrderedDict
from argparse import ArgumentParser
import localcider
from localcider.sequenceParameters import SequenceParameters
import math
import scipy.stats as stats

def get_kappa(seq,type1,type2):
    blobsz=[5]
    kappab=[]
    for b in blobsz:
        # Get full sequence asymmetry
        count1=0
        for res in type1: 
            count1=count1+seq.count(res)
        
        count2=0    
        for res in type2: 
            count2=count2+seq.count(res)
            
        sigAll=((count1/len(seq))-(count2/len(seq)))**2/((count1/len(seq))+(count2/len(seq)))
        
        # Get asymmetry for each blob
        sigX=[]
        for x in range(0,len(seq)-b+1):
            subseq=seq[x:x+b]
            
            count1=0
            for res in type1:
                count1=count1+subseq.count(res)
                
            count2=0
            for res in type2:
                count2=count2+subseq.count(res)
            
            if count1+count2==0:
                sigX.append(0)
            else:
                sigX.append(((count1/b)-(count2/b))**2/((count1/b)+(count2/b)))
        
        asym=[]
        for x in range(0,len(sigX)):
            asym.append((sigX[x]-sigAll)**2)
        
        kappab.append(np.mean(asym))
   
    kappa=np.mean(kappab)
    return kappa 

    
def get_omega(seq,type1):
    blobsz=[5]
    omegab=[]
    for b in blobsz:
        # Get full sequence asymmetry
        count=0
        for res in type1: 
            count=count+seq.count(res)
            
        sigAll=((count/len(seq))-(1-(count/len(seq))))**2
        
        # Get asymmetry for each blob
        sigX=[]

        for x in range(0,len(seq)-b+1):
            count=0
            subseq=seq[x:x+b]
            for res in type1:
                count=count+subseq.count(res)
                
            sigX.append(((count/b)-(1-(count/b)))**2)
        
        asym=[]
        for x in range(0,len(sigX)):
            asym.append((sigX[x]-sigAll)**2)
        
        omegab.append(np.mean(asym))
   
    omega=np.mean(omegab)
    return omega
        
def get_org_seq_vals(myseq,typeall,fracsall):
    org_seq_arr = np.zeros((len(typeall),len(typeall)))
    
    for count1 in range(0,len(typeall)):
        type1 = typeall[count1]

        for count2 in range(count1,len(typeall)):
            type2 = typeall[count2]

            if type1 == type2 and fracsall[count1]>0.10:
                org_seq_arr[count1, count2]=get_omega(myseq,type1)
                
            if type1 != type2 and fracsall[count1]>0.10 and fracsall[count2]>0.10:
                org_seq_arr[count1, count2]=get_kappa(myseq,type1,type2)
    
    org_seq_1d=org_seq_arr.reshape([1, len(typeall)**2])
    
    return org_seq_1d

def get_scramble_seqs_vals(myseq,num_seqs,typeall,fracsall):
    
    currseq=[]
    allseqs=[]
    scr_vals=np.zeros((num_seqs,len(typeall)**2))

    random.seed(42)
    np.random.seed(42)
    
    for x in range(0,num_seqs):
        currseq=''.join(random.sample(myseq,len(myseq)))

        scr_seq_arr = np.zeros((len(typeall),len(typeall)))
    
        for count1 in range(0,len(typeall)):
            type1 = typeall[count1]

            for count2 in range(count1,len(typeall)):
                type2 = typeall[count2]

                if type1 == type2 and fracsall[count1]>0.10:
                    scr_seq_arr[count1, count2]=get_omega(currseq,type1)

                if type1 != type2 and fracsall[count1]>0.10 and fracsall[count2]>0.10:
                    scr_seq_arr[count1, count2]=get_kappa(currseq,type1,type2)
        
        scr_vals[x,0:len(typeall)**2] = scr_seq_arr.reshape([1, len(typeall)**2])
        allseqs.append(currseq)
    
    #fit to a gamma distribution and obtain mean and variance
    alpha=[]
    beta=[]
    amean=[]
    avar=[]
    
    scr_vals_t=scr_vals.transpose()
    scr_vals_row = scr_vals_t.shape[0]
    for i in range(0,scr_vals_row):   
        fit_alpha, fit_loc, fit_beta = stats.gamma.fit(scr_vals_t[i,:])

        cmean = stats.gamma.mean(fit_alpha,fit_loc,fit_beta)
        cvar = stats.gamma.var(fit_alpha,fit_loc,fit_beta)
        alpha.append(fit_alpha)
        beta.append(fit_beta)
        amean.append(cmean)
        avar.append(cvar)
 
    return [alpha,amean,avar,scr_vals,allseqs]

def nardini(job_num, seq):
    #print("seq # = " + str(seq))
    #nard_seq = SequenceParameters(seq)
    #nard_val = nard_seq.calculate_zscore(num_scrambles=100000, random_seed=None)
    #oseq, sseq, index, zm, sm = nard_val['seq1']

    myseq = seq
    num_seqs=100000
    pol=['S','T','N','Q','C','H']
    hyd=['I','L','M','V']
    pos=['R','K']
    neg=['E','D']
    aro=['F','W','Y']
    ala=['A']
    pro=['P']
    gly=['G']
    typeall=[pol,hyd,pos,neg,aro,ala,pro,gly]
    type_checklist = []

    zvecdb=np.zeros((1,len(typeall)**2))

    fracsall=[]
    for type1 in typeall:
        mycount=0
        for res in type1:
            mycount=mycount+myseq.count(res)
        fracsall.append(mycount/len(seq))
    
    for i in range(len(typeall)):
        if fracsall[i] > 0.1:
            type_checklist.append(True)
        else:
            type_checklist.append(False)

    myarr=get_org_seq_vals(seq,typeall,fracsall)

# Returns mean of scrambles, std of scrambles, all values in a number of scramble x 64 list, all scramble sequences, and kappa values of scrambles
    [alpha,amean,avar,allscrvals,allscrseqs]=get_scramble_seqs_vals(seq,num_seqs,typeall,fracsall)
    
    for x in range(0,myarr.shape[1]):
        if myarr[0,x]==0:
            zvecdb[0,x]=0
        else:
            zvecdb[0,x]=(myarr[0,x]-amean[x])/math.sqrt(avar[x])

    zm = zvecdb ###?????###
    zm_list = np.reshape(zm,(1,64))

    # it may be better to save the results rather than have them returned
    np.save('zscores/job-{:04d}.npy'.format(job_num), zm_list)
    return zm_list

def main():
    parser = ArgumentParser()
    parser.add_argument('-o', '--seqs', help='The filepath containing the sequences.', type=str)
    parser.add_argument('-p', '--num-processes', help='The number of parallel processes to use in the multiprocessing pool.', type=int, default=1)
    args = parser.parse_args()


    # check if the arrays directory exists - create if nonexistent
    if not os.path.exists('zscores'):
        os.makedirs('zscores')

    seqs = np.loadtxt(args.seqs, dtype=str)

    # Create a worker pool which will process N jobs in M amounts. This means that at
    # any given time, a maximum of M processes will be performing work.
    pool = mp.Pool(processes=args.num_processes)

    # As we're consuming N sequences, that is a better iterator than the number of jobs (removed)
    results = OrderedDict()
    for job_num, npy_seq in enumerate(seqs, start=1):
        seq = str(npy_seq)
        # Store the job for asynchronous execution.
        # Here we store the `pool.apply_async` result (which has a return type of `AsyncResult`)
        # as we can later consult the jobs and query them as they compute and the results "arrive" (if desired).
        async_result = pool.apply_async(nardini, args=(job_num, seq,))  # the extra `,` in `args` is not a typo!
        start_time = datetime.now()
        results[seq] = (job_num, start_time, async_result)

    pool.close()  # Tell the pool that we're no longer accepting jobs
    pool.join()   # Now, perform the asynchronous calculations. The results will arrive out of order (as expected).

    # To get the results as they complete, we need to use `AsynResult.get()`. This must be done
    # after the loop has been closed, otherwise each call will block and there would be no
    # performance gains - it would be as though you had only 1 processor
    actual_results = OrderedDict()
    for seq in results:
        job_num, start_time, async_result = results[seq]
        result = async_result.get()
        current_time = datetime.now()
        actual_results[seq] = (job_num, start_time, current_time, result)

    np.save('zscores/all-analysis.npy', actual_results)


if __name__ == '__main__':
    main()
