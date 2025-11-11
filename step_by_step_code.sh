# Step 53

cd /your/working/directory
mkdir fastq
cp /source/directory/*.fasta.gz fastq
gunzip fastq/*


# Step 54a

mkdir STARv2.7.11b_genome_index_mm10
module load star/2.7.11b

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir STARv2.7.11b_genome_index_mm10 --genomeFastaFiles GRCm38.p6.genome.fa --sjdbGTFfile gencode.vM25.basic.annotation.gtf


# Step 54b

mkdir STAR_output

STAR --genomeDir STARv2.7.11b_genome_index_mm10 --runThreadN 12 --readFilesIn fastq/${sample_name}.fastq --outFilterMultimapNmax 1 --outFileNamePrefix STAR_output /${sample_name} --outSAMtype SAM --outReadsUnmapped Fastx


# Step 55

mkdir stranded_sam
module load samtools/1.21

samtools view -f 0x10 STAR_output/${sample_name}Aligned.out.sam  > stranded_sam/${sample_name}.sam_pos

samtools view -F 0x10 STAR_output/${sample_name}Aligned.out.sam  > stranded_sam/${sample_name}.sam_neg


# Step 56

mkdir wiggles

# Make wiggle file for negative strand. Specify desired display color (e.g., red), paired-end data flag (n/y), and log10 normalization flag (n/y)

perl bigsam_to_wig_mm10_wcigar2.pl stranded_sam/${sample_name}.sam_neg wiggles/${sample_name}_neg red n n 50

gzip ${sample_name}_neg.wig

# Make wiggle file for positive strand. Specify desired display color (e.g., blue), paired-end data flag (n/y), and log10 normalization flag (n/y)

perl bigsam_to_wig_mm10_wcigar2.pl stranded_sam/${sample_name}.sam_pos wiggles/${sample_name}_pos blue n n 50

gzip ${sample_name}_pos.wig


# Optional:

perl scale_wiggle_rpm.pl wiggles/${sample_name}_neg.wig <total-reads> wiggles/${sample_name}_neg_RPM.wig
gzip wiggles/${sample_name}_neg_RPM.wig

perl scale_wiggle_rpm.pl wiggles/${sample_name}_pos.wig <total-reads> wiggles/${sample_name}_pos_RPM.wig
gzip wiggles/${sample_name}_pos_RPM.wig


# Step 60

module load star/2.7.11b
module load samtools/1.21

# Use UNIX to obtain a comma separated list of the names of all IgG fastq files, “rip_filenames”, with all white spaces and newline characters removed
rip_filenames=$(ls -m *.fastq | sed -r 's/\s+//g' | tr -d '\n')

# Use STAR to align to the reference genome all IgG replicates in the rip_filenames list 
STAR --genomeDir /your/STARv2.7.11b_genome_index_folder --runThreadN 12 --outFilterMultimapNmax 1 --outFileNamePrefix igg_ --readFilesIn $rip_filenames

# Split the filtered SAM file by strand
samtools view -h -F 0x10 igg_Aligned.out.sam > igg_Aligned_neg.out.sam
samtools view -h -f 0x10 igg_Aligned.out.sam > igg_Aligned_pos.out.sam


# Step 61

# Use UNIX to obtain a comma separated list of the names of all RIP FASTQ files, “rip_filenames”, with all white spaces and newline characters removed 
rip_filenames=$(ls -m *.fastq | sed -r 's/\s+//g' | tr -d '\n')

# Use STAR to perform a single alignment of all RIP replicates merged together, and assign the prefix of the output file as “yourRIPname_”
STAR --genomeDir /your/STARv2.7.11b_genome_index_folder --runThreadN 12 --outFilterMultimapNmax 1 --outFileNamePrefix yourRIPname_ --readFilesIn $rip_filenames

# The next series of commands will perform STAR alignment on each individual replicate of RIP data in the directory
# Count the total number of replicates (assuming each FASTQ file represents one separate biological replicate), and save that number as “rip_file_num”
rip_file_num=$(ls *.fastq | wc -l)

if [ $rip_file_num -gt 1 ]
# Only perform these subsequent steps if the replicate number is greater than 1, and not when there is only 1 replicate

then
# Convert rip_filenames into an array for the loop in the subsequent step
    IFS=',' read -r -a rip_arr <<< "$rip_filenames" 

# Loop on indices of array
for i in ${!rip_arr[@]}; do 
        # For each file name
        rip_name="${rip_arr[$i]}" 
        # Save the FASTQ name and its corresponding replicate number
        echo The rip${i} file is $rip_name. >> yourRIPname_rip_info.txt 
        # Index+1 so that it's 1-based for the file name
        j=$(($i + 1)) 

        # Align each RIP file to the reference genome, and store each alignment using a filename structure “yourRIPname_rip${j}_”
        STAR --genomeDir /your/STARv2.7.11b_genome_index_folder --runThreadN 12 --outFilterMultimapNmax 1 --outFileNamePrefix yourRIPname_rip${j}_ --readFilesIn $rip_name

    done
fi


# Step 62

# Split the filtered SAM by strand
samtools view -h -F 0x10 yourRIPname_Aligned.out.sam > yourRIPname_Aligned_neg.out.sam
samtools view -h -f 0x10 yourRIPname_Aligned.out.sam > yourRIPname_Aligned_pos.out.sam

# Load randomization Perl script macs_strand_rand_sam.pl
cp /path/to/macs_strand_rand_sam.pl .

# Randomize strand value (+ or -) for each aligned read within each stranded file. This is performed because MACS is designed to read non-strand-specific ChIP-seq data and uses the average distance between positive- and negative-stranded alignments to center its peak calls.

perl macs_strand_rand_sam.pl yourRIPname_Aligned_pos.out.sam yourRIPname_Aligned_pos.out.rand.sam
perl macs_strand_rand_sam.pl yourRIPname_Aligned_neg.out.sam yourRIPname_Aligned_neg.out.rand.sam

# Convert sam to bam to avoid errors in MACS peak calling
samtools view -S -b yourRIPname_Aligned_pos.out.rand.sam > yourRIPname_Aligned_pos.out.rand.bam
samtools view -S -b yourRIPname_Aligned_neg.out.rand.sam > yourRIPname_Aligned_neg.out.rand.bam

# Use MACS to identify potential peak regions from the bam files

module load macs/2.2.7.1

macs2 callpeak -t yourRIPname_Aligned_pos.out.rand.bam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir yourRIPname_pos_peaks -n yourRIPname_pos
macs2 callpeak -t yourRIPname_Aligned_neg.out.rand.bam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir yourRIPname_neg_peaks -n yourRIPname_neg 


# Step 63

module load subread/2.0.6

# Use the output .broadPeak from MACS peak calling (Step 62)
# Convert the broadPeak file to a SAF file in preparation for input into featureCounts
cat yourRIPname_neg_peaks.broadPeak | awk -v OFS="\t" '{print $1,$2,$3,"-"}' > yourRIPname_saf.txt
cat yourRIPname_pos_peaks.broadPeak | awk -v OFS="\t" '{print $1,$2,$3,"+"}' >> yourRIPname_saf.txt

# Add column names to the first line of the file and save it with the suffix .saf 
cat yourRIPname_saf.txt | awk -v OFS="\t" '{print NR,$0}' | sed '1 i\GeneID\tChr\tStart\tEnd\tStrand' > yourRIPname_peaks.saf

# Count reads under each potential MACS peak using the total merged replicate file
featureCounts -s 2 -F SAF -a yourRIPname_peaks.saf -o yourRIPname_fc yourRIPname_Aligned.out.sam

# Count reads under each potential MACS peak using the merged IgG control file
# Using the path in Step 60: /your/path/to/igg_data_folder
featureCounts -s 2 -F SAF -a yourRIPname_peaks.saf -o yourRIPname_igg_fc /your/path/to/igg_data_folder/igg_Aligned.out.sam

# Extract the read counts under each potential MACS peak from the IgG file. These counts will be appended to yourRIPname_fc later, along with other information 
cut -f7 yourRIPname_igg_fc > yourRIPname_igg_counts.txt

if [ $rip_file_num -gt 1 ]
# Proceed only if there is more than 1 RIP replicate 
then
    # Loop on indices of array rip_arr
for i in ${!rip_arr[@]}; do 
        # Index+1 so that it's 1-based for the file name
        j=$(($i + 1)) 
        # Count reads under each MACS peak for individual replicates
        featureCounts -s 2 -F SAF -a yourRIPname_peaks.saf -o yourRIPname_rip${j}_fc yourRIPname_rip${j}_Aligned.out.sam
        
        # Extract the read counts under each potential MACS peak
        # These counts will be appended to yourRIPname_fc later
        cut -f7 yourRIPname_rip${j}_fc > yourRIPname_rip${j}_counts.txt
        
        # Prepare for concatenation of count files
        # Build a string called “rip_count_filename” that contains all of the replicates _counts.txt files, separated by spaces
        if [ $j == 1 ]
        then
        	# Initialize the string properly for the first replicate
        	rip_count_filename="yourRIPname_rip1_counts.txt"
        else
        	# Append the rest of the replicate file names
        	rip_count_filename="${rip_count_filename} yourRIPname_rip${j}_counts.txt"
        fi
    done

    # Make summary file by concatenating in the following order: counts from all RIP replicates merged together, counts from the merged IgG control, counts from the individual RIP replicates
    paste yourRIPname_fc yourRIPname_igg_counts.txt $rip_count_filename | sed 1d > yourRIPname_fc.txt

else
# For situations in which there is only 1 replicate
# Make a summary file by concatenating: counts from the single RIP replicate, counts from the merged IgG control 
    paste yourRIPname_fc yourRIPname_igg_counts.txt | sed 1d > yourRIPname_fc.txt

fi


# Step 64

# RPM conversion for IgG control:
# Copy over all IgG FASTQ files to count total read number
cp /your/path/to/igg_data_folder/*igg*.fastq .
# Count the total line number across all *igg*.fastq files
fastq_lin_num_igg=$(wc -l *igg*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
# Total read count is the total line number divided by 4
(( reads_igg = fastq_lin_num_igg / 4 ))

# Count the total line number across all *yourRIPname*.fastq files
fastq_lin_num=$(wc -l *yourRIPname*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
# Total read count is the total line number divided by 4
(( reads = fastq_lin_num / 4 ))

# RPM conversion for all replicates, $7 is the 7th column that contains the raw count
cat yourRIPname_fc.txt | sed 1d | awk -v OFS="\t" -v reads=$reads '{print $0, (($7*1000000)/reads)}' > yourRIPname_fc_allRpm.txt

# RPM conversion for individual replicates:

if [ $rip_file_num -gt 1 ]
# Proceed if there is more than 1 replicate for the RIP-of-interest

then
    # Loop on indices of array rip_arr
for i in ${!rip_arr[@]}; do
        # Index+1 so that it's 1-based for the file name 
        j=$(($i + 1)) 
        # Count total line number of each individual replicate FASTQ 
        fastq_lin_num=$(wc -l ${rip_arr[$i]} | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
        # Total read count is the total line divided by 4
        (( reads = fastq_lin_num / 4 ))
        # RPM conversion for each replicate, $1 is the 1st column that contains the raw count
        cat yourRIPname_rip${j}_counts.txt | sed 1,2d | awk -v reads=$reads '{print (($1*1000000)/reads)}' > yourRIPname_rip${j}_rpm.txt

        # Prepare to concatenate rpm files and file header
        if [ $j == 1 ]
        then
        	# Initiate the table properly for the first replicate
        	rip_count_header="yourRIPname_rip1_counts"
        	rip_rpm_header="yourRIPname_rip1_rpm"
        	rip_rpm_filename="yourRIPname_rip1_rpm.txt"
        else
		rip_count_header="${rip_count_header}\tyourRIPname_rip${j}_counts"
        	rip_rpm_header="${rip_rpm_header}\tyourRIPname_rip${j}_rpm"
        	rip_rpm_filename="${rip_rpm_filename} yourRIPname_rip${j}_rpm.txt"
        fi
    done

    # Make summary file with the following information: geneID, chr, start, end, strand, length, read count across all RIP replicates, total read count for IgG control, read counts for each RIP replicate, rpm across all RIP replicates, IgG rpm, and rpm for each RIP replicate
    paste yourRIPname_fc_allRpm.txt yourRIPname_igg_rpm.txt $rip_rpm_filename | sed "1 i\GeneID\tChr\tStart\tEnd\tStrand\tLength\tyourRIPname_counts\tigg_counts\t${rip_count_header}\tyourRIPname_rpm\tigg_rpm\t${rip_rpm_header}" > yourRIPname_fc_rpm.txt

else
    # If there is only 1 replicate of the RIP-of-interest
    # Make summary file with the following information: chr, start, end, strand, length, RIP read counts, IgG read counts, RIP rpm, IgG rpm
    paste yourRIPname_fc_allRpm.txt yourRIPname_igg_rpm.txt | sed "1 i\GeneID\tChr\tStart\tEnd\tStrand\tLength\tyourRIPname_counts\tigg_counts\tyourRIPname_rpm\tigg_rpm" > yourRIPname_fc_rpm.txt

fi


# Step 65

if [ ${rip_file_num} -gt 1 ]
# Proceed if there is more than 1 replicate

then
# Append a new column (last column) and store 0 for all lines
# This step creates a new column in the file that
# will be used to denote for each peak, 
# the number of replicates where RIP-over-IgG RPM values are >2 
    cat yourRIPname_fc_rpm.txt | sed 1d | awk -v OFS="\t" '{print $0, 0}' > yourRIPname_ripi_2igg.txt

((igg_col = 10 + rip_file_num)) # Column number of igg_rpm

    for i in $(seq 1 $rip_file_num) # For each replicate
    do
        ((rip_col = 10 + rip_file_num + i)) # Column number of rip_rpm
        # (rip_rpm > 2*igg_rpm) -> if TRUE, increment count in the last column
        awk -i inplace -v OFS="\t" -v c1=$rip_col -v c2=$igg_col  '{if($c1 > 2 * $c2) {$NF=$NF+1;print} else {print $0}}' yourRIPname_ripi_2igg.txt
    done
    
    # Find peaks whose RIP-over-IgG RPM values are >2 
    # in >=2 replicates
    cat yourRIPname_ripi_2igg.txt | awk -v OFS="\t" '{if($NF > 1) {print $0}}' > yourRIPname_fc_rpm_2igg_filtered.txt

    # Write BED files based on yourRIPname_fc_rpm_2igg_filtered.txt which identify the peaks whose RIP-over-IgG RPM values are >2 in >=2 replicates

    cat yourRIPname_fc_rpm_2igg_filtered.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > yourRIPname_peaks_2igg_2reps.bed
  
fi


# Step 66

module load bedtools/2.31.1

bedtools getfasta -s -fi <input-genome-fasta> -bed <input-peak-bed-file> -fo <output-peak-fasta>


# Step 67

python3 create_bkg_seq.py genome-fasta-file.fa input-peak-seq-fasta.fa output-control-seq-fasta.fa

# Step 68

module load meme/5.5.7

streme --oc <output-name-prefix-only> -n <control-seq-fasta> --maxw 12 --minw 4 --nmotifs 20 --rna --p <peak-seq-fasta>
