from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from functions import read_FASTA
import time

oriC = 'gatctatttatttagagatctgttctattgtgatctcttattaggatcgcactgccctgtggataacaaggatccggcttttaagatcaacaacctggaaaggatcattaactgtgaatgatcggtgatcctggaccgtataagctgggatcagaatgaggggttatacacaactcaaaaactgaacaacagttgttctttggataactaccggttgatccaagcttcctga'.upper()
print(oriC)

seq = read_FASTA('test_fastas\Escherichia_coli_K_12.fna')[1][3920600:3927040]
print(seq)

# aligner = PairwiseAligner()
# # aligner.mode = "local"
# # aligner.match_score = 1.0
# alignments = aligner.align(oriC, seq)

# for alignment in alignments:
#     print(alignment)