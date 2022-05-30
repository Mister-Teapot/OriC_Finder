from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from oriC_Finder import read_FASTA
import time

DnaA = Seq("MENILDLWNQALAQIEKKLSKPSFETWMKSTKAHSLQGDTLTITAPNEFARDWLESRYLHLIADTIYELTGEELSIKFVIPQNQDVEDFMPKPQVKKAVKEDTSDFPQNMLNPKYTFDTFVIGSGNRFAHAASLAVAEAPAKAYNPLFIYGGVGLGKTHLMHAIGHYVIDHNPSAKVVYLSSEKFTNEFINSIRDNKAVDFRNRYRNVDVLLIDDIQFLAGKEQTQEEFFHTFNTLHEESKQIVISSDRPPKEIPTLEDRLRSRFEWGLITDITPPDLETRIAILRKKAKAEGLDIPNEVMLYIANQIDSNIRELEGALIRVVAYSSLINKDINADLAAEALKDIIPSSKPKVITIKEIQRVVGQQFNIKLEDFKAKKRTKSVAFPRQIAMYLSREMTDSSLPKIGEEFGGRDHTTVIHAHEKISKLLADDEQLQQHVKEIKEQLK")

start = time.time()
_, frame_1 = read_FASTA('test_fastas\Bacillus_subtilis_168.fna')
frame_2 = frame_1[1:] + frame_1[0]
frame_3 = frame_2[1:] + frame_2[0]
read = time.time() - start
print('Done reading:', read)

start = time.time()
test_1 = Seq(frame_1)
test_2 = Seq(frame_2)
test_3 = Seq(frame_3)

aa_1 = test_1.translate()
aa_2 = test_2.translate()
aa_3 = test_3.translate()

convert = time.time() - start
print('Done converting to aa', convert)

aligner = PairwiseAligner()
aligner.mode = "local"
# aligner.match_score = 1.0

start = time.time()
alignments = aligner.align(aa_1, DnaA)
align = time.time() - start
print('Done aligning', align)

for alignment in alignments:
    print(alignment)