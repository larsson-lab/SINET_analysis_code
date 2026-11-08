#!/bin/bash
set -euox pipefail

for i in {1..500}
do
        echo "Run" $i
	path='/home/alice/Timing_publication_remake/code_repository/results/matrixgenerator/bootstrap/run'$i'/'
	SigProfilerMatrixGenerator matrix_generator 'SINET' 'GRCh38' $path
	mv /home/alice/Timing_publication_remake/code_repository/results/matrixgenerator/bootstrap/run$i/output/SBS/SINET.SBS96.all /home/alice/Timing_publication_remake/code_repository/results/matrixgenerator/bootstrap/run$i/ 
	rm -r /home/alice/Timing_publication_remake/code_repository/results/matrixgenerator/bootstrap/run$i/logs/
	rm -r /home/alice/Timing_publication_remake/code_repository/results/matrixgenerator/bootstrap/run$i/output/
	rm -r /home/alice/Timing_publication_remake/code_repository/results/matrixgenerator/bootstrap/run$i/input
done