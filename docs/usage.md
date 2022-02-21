# Usage

### 1. File preparation
Input
1. fasta file of paired-end reads
2. fasta file of assembled contigs

Output: BAM file and pileup file

Map paired-end reads to assembled contigs to generate BAM file and use samtools to generate pileup file
```
bwa index $contig_file
bwa mem -a -t 8 $contig_file $read1 $read2 | samtools view -h -q 10 -m 50 -F 4 -b | samtools sort > $bam_file
samtools mpileup -C 50 -A -f $contig_file $bam_file |  awk '$3 != "N"' > $pileup_file
```

### 2. Feature extraction

Input: BAM file, Pileup file, assembled contigs

#### For metagenomics

```
# Step 1: extract features [output file: feature_matrix/window_fea_matrix.txt,feature_matrix/contig_fea_matrix.txt]

metaMIC extract_feature --bam $bam_file -c $contig_file -o $output_dir --pileup $pileup_file -m meta
```
#### For isolate genomes

```
# Step 1: extract features [output file: feature_matrix/window_fea_matrix.txt]

metaMIC extract_feature --bam $bam_file -c $contig_file -o $output_dir --pileup $pileup_file -m single
```

### 3. Predict and correct misassemblies 

#### For metagenomics

```
# Step 2: misassembly breakpoint identification and correction;
# output directory must be same as the above $output_dir
# [output file: metaMIC_corrected_contigs.fa, misassembly_breakpoint.txt, anomaly_score.txt]

metaMIC predict -c $contig_file -o $output_dir -a MEGAHIT -m meta
```

#### For isolate genomes

```
# Step 2: misassembly breakpoint identification and correction;
# output directory must be same as the above $output_dir
# [output file: metaMIC_corrected_contigs.fa, misassembly_breakpoint.txt, anomaly_score.txt]

metaMIC predict -c $contig_file -o $output_dir -m single
```

### 4. Training on new datasets

If you want to generate a new training model on a novel dataset. Contig labels and name of new training models should be provided. The Step 1 is the same as above, then the contig_fea_matrix.txt will be used as the training datasets.


```
# Step 1: extract features [output file: feature_matrix/window_fea_matrix.txt,feature_matrix/contig_fea_matrix.txt]

metaMIC extract_feature --bam $bam_file -c $contig_file -o $output_dir --pileup $pileup_file -m meta

# Step 2: Generating new training models
# output directory must be same as the above $output_dir
# format of contig_label file should be:
# contig1\t0
# contig2\t1
# contig3\t0
# ...
# contig100\t0

metaMIC train -o $output_dir -a $New_model_name --label $contig_label 
```
