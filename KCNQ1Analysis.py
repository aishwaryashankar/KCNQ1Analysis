from Bio import AlignIO
from Bio import SeqIO
import subprocess
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator

# NOTE: This Python module, KCNQ1Analysis.py, was written and documented by Aishwarya Shankar
#       Jennifer Benbow researched and collected the needed KCNQ1 sequence data (kcnq1sequences.txt) to use as the input file for this program


# The muscle variable refers to the MUSCLE executable file which will be used to run MSA on chosen sequences in a child process
muscle = "muscle"

# The text file containing the desired sequences should be in the same directory as this program and the MUSCLE executable
input = "./sequences.txt"

# MSA result will be output to the alignment.fasta file
output = "alignment.fasta"

# Conducting MSA thru the MUSCLE tool in a child process (subprocess)
muscle = subprocess.check_output([muscle,"-in", input,"-out",output])

# Printing the MSA alignment (summary) on command prompt
print("\n\n")
print("MSA Alignment\n")
alignments = AlignIO.parse(open('alignment.fasta'), 'fasta')
for each in alignments:
    print(each)

# Writing MSA results from fasta format to a phylip file (align.phylip) to be used to construct phylogenetic trees
file = SeqIO.parse('alignment.fasta','fasta')
count = SeqIO.write(file, 'align.phylip','phylip')

# Reading from phylip file to construct the distance matrix for UPGMA and NJ tree-building
alignment = AlignIO.read('align.phylip','phylip')
distanceCalculator = DistanceCalculator('identity')
distance_matrix = distanceCalculator.get_distance(alignment)
print("\n\nTHIS IS THE DISTANCE MATRIX TO CONSTRUCT THE PHYLOGENETIC TREES\n")
print(distance_matrix)

# Constructing and drawing the UPGMA and Neighbor-Joining (distance-based) trees using the previously derived distance matrix
treeConstructor = DistanceTreeConstructor()
upgma_phyl_tree = treeConstructor.upgma(distance_matrix)
print("\n\nUPGMA PHYLOGENETIC TREE DERIVED FROM MUSCLE MSA\n")
Phylo.draw_ascii(upgma_phyl_tree)
nj_phyl_tree = DistanceTreeConstructor().nj(distance_matrix)
print("\n\nNEIGHBOR-JOINING (NJ) PHYLOGENETIC TREE DERIVED FROM MUSCLE MSA\n")
Phylo.draw_ascii(nj_phyl_tree)