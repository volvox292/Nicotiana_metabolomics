#uses 22 cores and the in-silico library split into 22 parts use sample.mgf as input argument
# will run dbsearch.py as subprocess

import os                                                                       
from multiprocessing import Pool   
import pandas as pd
import time
import sys
start_time = time.time()
processes=[]
#sample="GAQ-06-01-28012022-py.mgf"
sample=str(sys.argv[1])
for i,n in enumerate(range(22)):
    processes.append("dbsearch.py "+str(n)+" {}".format(sample))
    
processes=tuple(processes)                                   
                                                  
                                                                                
def run_process(process):                                                             
    os.system('python {}'.format(process))                                       
                                                                                
if __name__ == '__main__':
    with Pool(processes=22) as pool:
        pool.map(run_process, processes)
        try:
            while pool.is_alive()==True:
                time.sleep(60)
        except:
            print("Workers are done!")
            print("Now merging results...")
            results_names=["C:\\Users\\Public\\results\\p{}-results.tsv".format(n)for i,n in enumerate(range(22))]
            df = pd.concat([pd.read_csv(f, low_memory=False, sep="\t") for f in results_names],ignore_index=True)
            df=df.rename(columns={"cousine-score Hit 1": "cosine"})
            agg_func_selection = {"cosine":['first'], "m/z":["first"],"Library Hit 1":["first"],"m/z 1":["first"],"numoffpeaks 1":["first"],"Smiles 1":["first"]}
            result=df.sort_values(["cosine"],ascending=False).groupby("Feature_id").agg(agg_func_selection)
            result.to_csv("C:\\Users\\Public\\results\\results.tsv", sep="\t")
            print("Done!")
            print("--- %s minutes ---" % ((time.time() - start_time)/60))