import sys
import numpy as np

if len(sys.argv) != 4:
    raise SystemExit(
        f"Usage: {sys.argv[0]} <genome-fasta-file> <input-peak-seq-fasta> <output-control-seq-fasta>"
    )


genome_fasta = sys.argv[1]
input_peak_fasta = sys.argv[2]
output_control_fasta = sys.argv[3]


A_count = 0
T_count = 0
G_count = 0
C_count = 0

with open(genome_fasta, "r") as genome_file:
    for line in genome_file:
        if line[0] != ">":
            for char in line:
                if char == "A":
                    A_count += 1
                elif char == "T":
                    T_count += 1
                elif char == "G":
                    G_count += 1
                elif char == "C":
                    C_count += 1

print(A_count)
print(T_count)
print(G_count)
print(C_count)

total_nts = A_count + T_count + G_count + C_count
A_frac = A_count/total_nts
T_frac = T_count/total_nts
G_frac = G_count/total_nts
C_frac = C_count/total_nts

nts = ["A","T","G","C"]

def create_background_seq_file(control_file, peak_seq_file):
    counter = 0
    for line in peak_seq_file:
        if line[0] != ">":
            counter += 1
            control_file.write(">%s\n" % counter)
            control_file.write(str("%s\n" % "".join(np.random.choice(nts,len(line.strip()),p=[A_frac,T_frac,G_frac,C_frac]))))

with open(output_control_fasta, "w") as out_control_file:
    with open(input_peak_fasta, "r") as in_peak_file:
        create_background_seq_file(out_control_file, in_peak_file)
