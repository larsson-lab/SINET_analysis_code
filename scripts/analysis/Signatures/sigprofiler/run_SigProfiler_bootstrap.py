from SigProfilerAssignment import Analyzer as Analyze
import os
import shutil
import multiprocessing
def runAssignment(i):

    samples = '/home/alice/Timing_publication_remake/code_repository/results/matrixgenerator/bootstrap/run' + str(i) + '/SINET.SBS96.all'
    output = '/home/alice/Timing_publication_remake/code_repository/results/sigprofiler/bootstrap/run' + str(i) + '/'
    Analyze.cosmic_fit(samples, output,
                    cosmic_version=3.2, 
                    genome_build="GRCh38", 
                    make_plots=False)
    
    shutil.rmtree(output + '/Assignment_Solution/Signatures')
    shutil.rmtree(output + '/Assignment_Solution/Solution_Stats')
    os.remove(output + '/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.txt')
    print(str(i) + ' version 3.2 is finished')

pool_obj = multiprocessing.Pool(processes=4)
ans = pool_obj.map(runAssignment,range(1, 501))
pool_obj.close()
