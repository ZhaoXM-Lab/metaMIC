# metaMIC: Reference-free Misassembly Identification and Correction of metagenomic assemblies
metaMIC is a fully automated tool for identifying and correcting misassemblies of (meta)genomic assemblies with the following three steps. Firstly, metaMIC extracts various types of features from the alignment between paired-end sequencing reads and the assembled contigs. Secondly, the features extracted in the first step will be used as input of a random forest classifier for identifying misassembled metagenomic assemblies. Thirdly, metaMIC will localize misassembly breakpoints for each misassembled contig and then correct misassemblies by splitting into parts at the breakpoints.


## Requirements and Installation
Make sure you have the dependencies below installed and accessible in your $PATH.

### Prepare dependencies

- [python 3.6-3.9](https://www.python.org/downloads/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [pysam](https://pypi.org/project/pysam/0.8.4/)
- [biopython](https://pypi.org/project/biopython/)
- [bwa 0.7.17](https://sourceforge.net/projects/bio-bwa/files/)
- [samtools 1.9](https://sourceforge.net/projects/samtools/files/samtools/)
- [jellyfish](http://www.cbcb.umd.edu/software/jellyfish/)
- [scikit-learn](https://scikit-learn.org/stable/)
- [scipy](https://www.scipy.org/)

1. install python modules: pandas, numpy, pysam
```
conda install -c bioconda pandas numpy pysam biopython  
```
or

```
pip install pysam pandas numpy biopython
```

2. download and install samtools, bwa, jellyfish
```
conda install -c bioconda samtools bwa jellyfish  
```
or

```
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -jxvf samtools-1.9.tar.bz2
cd samtools-1.9
./configure
make
make install
export PATH=`pwd`:$PATH
```

```
wget https://sourceforge.net/projects/bio-bwa/files/latest/download/bwa-0.7.17.tar.bz2
tar -jxvf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
export PATH=`pwd`:$PATH
```

```
wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.10.tar.gz
tar zxvf jellyfish-1.1.10.tar.gz
cd jellyfish-1.1.10
./configure
make
make install
export PATH=`pwd`/bin:$PATH
```

### Installation

##### Install metaMIC via git

```
git clone https://github.com/ZhaoXM-Lab/metaMIC.git
```
- install and download [training models](https://zenodo.org/record/4781819#.YKm_omYzZTY) 

```
cd metaMIC
python setup.py install
metaMIC -h

# downloading models
metaMIC download_model
```

## Quick Start
- First map paired-end reads to assembled contigs

```
bwa index $contig_file
bwa mem -a -t 8 $contig_file $read1 $read2 | samtools view -h -q 10 -m 50 -F 4 -b | samtools sort > $bam_file
```
- generate pileup file

```
samtools mpileup -C 50 -A -f $contig_file $bam_file |  awk '$3 != "N"' > $pileup_file
```
- run metaMIC

#### For metagenomics

```
# Step 1: extract features [output file: feature_matrix/window_fea_matrix.txt,feature_matrix/contig_fea_matrix.txt]

metaMIC extract_feature --bam $bam_file -c $contig_file -o $output_dir --pileup $pileup_file -m meta

# Step 2: misassembly breakpoint identification and correction;
# output directory must be same as the above $output_dir
# [output file: metaMIC_corrected_contigs.fa, misassembly_breakpoint.txt, anomaly_score.txt]

metaMIC predict -c $contig_file -o $output_dir -a MEGAHIT -m meta
```
#### For isolate genomes

```
# Step 1: extract features [output file: feature_matrix/window_fea_matrix.txt]

metaMIC extract_feature --bam $bam_file -c $contig_file -o $output_dir --pileup $pileup_file -m single

# Step 2: misassembly breakpoint identification and correction;
# output directory must be same as the above $output_dir
# [output file: metaMIC_corrected_contigs.fa, misassembly_breakpoint.txt, anomaly_score.txt]

metaMIC predict -c $contig_file -o $output_dir -m single
```
#### Training on new datasets

If you want to generate a new training model on a novel dataset. Contig labels and name of new training models should be provided.
The Step 1 is the same as above, then the contig_fea_matrix.txt will be used as the training datasets.


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


For more details about the usage of metaMIC, [read the docs](http:)

## Example

##### example data

- R.sphaeroides_pileup.out
- R.sphaeroides.bam
- velvet_ctg.fasta

```
cd example
sh download.sh
metaMIC extract_feature --pileup R.sphaeroides_pileup.out --bam R.sphaeroides.bam -c velvet_ctg.fasta -m single -o test
metaMIC predict -c velvet_ctg.fasta -o test -m single
```



## Output
The output folder will contain
1. Misassembly score for each contig: **metaMIC_contig_score.txt** (for only metagenomics)
2. predicted misassembly breakpoints for misassembled contigs: **misassembly_breakpoint.txt**
3. Anomaly scores for each position in contigs: **anomaly_score.txt**
4. Fasta file of corrected contigs: **metaMIC_corrected_contigs.fa**
5. Some intermediate files

For more details about the output, [read the docs](http:)

## Complete option list
metaMIC:

```
usage: metaMIC [-h]  ...

Reference-free Misassembly Identification and Correction of metagenomic
assemblies

optional arguments:
  -h, --help       show this help message and exit

metaMIC subcommands:

    extract_feature
                   Extract features from inputs.
    predict        Predict.
    train          Train model.
    

usage: metaMIC extract_feature [-h] [-t THREADS] [--bam BAMFILE] [--r1 READ1]
                               [--r2 READ2] [-p READ] -c ASSEMBLIES -o OUTPUT
                               --pileup PILEUP -m MODE [-l MIN_LENGTH]
                               [--samtools SAMTOOLS] [--jellyfish JELLYFISH]

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Maximum number of threads [default: 8]
  --bam BAMFILE         index bam file for alignment
  --r1 READ1            read1
  --r2 READ2            read2
  -p READ, --r READ     smart pairing (ignoring #2 fasta/q)
  -c ASSEMBLIES, --contig ASSEMBLIES
                        fasta file of assembled contigs
  -o OUTPUT, --output OUTPUT
                        output directory for metaMIC results
  --pileup PILEUP       path to pileup file [samtools mpileup]
  -m MODE, --mode MODE  Applied to single genomic/metagenomic assemblies
                        [meta/single]
  -l MIN_LENGTH, --mlen MIN_LENGTH
                        Minimum contig length [default: 5000bp]
  --samtools SAMTOOLS   path to samtools
  --jellyfish JELLYFISH
                        path to jellyfish
                        
usage: metaMIC predict [-h] -o OUTPUT -m MODE -c ASSEMBLIES [-a ASSEMBLER]
                       [-l MIN_LENGTH] [-s SPLIT_LENGTH] [--nb BREAK_COUNT]
                       [--rb BREAK_RATIO] [--at ANOMALY_THRED]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output directory for metaMIC results
  -m MODE, --mode MODE  Applied to single genomic/metagenomic assemblies
                        [meta/single]
  -c ASSEMBLIES, --contig ASSEMBLIES
                        fasta file of assembled contigs
  -a ASSEMBLER, --assembler ASSEMBLER
                        The assembler-specific model or user-trained model
                        used for assembled fasta file [MEGAHIT/IDBA_UD/metaSPAdes/[new
                        training model specified by users]]
  -l MIN_LENGTH, --mlen MIN_LENGTH
                        Minimum contig length [default: 5000bp]
  -s SPLIT_LENGTH, --slen SPLIT_LENGTH
                        Minimum length of split fragments [default: 1000bp]
  --nb BREAK_COUNT      Threshold of read breakpoint counts for correcting
                        misassemblies 
  --rb BREAK_RATIO      Threshold of read breakpoint ratio for correcting
                        misassemblies 
  --at ANOMALY_THRED    Threshold of anomaly score for correcting
                        misassemblies 
  --st SCORE_THRED      Threshold of contig score for correcting misassemblies                        
                        
usage: metaMIC train [-h] -o OUTPUT [--label LABEL] [-a ASSEMBLER]
                     [-t THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output directory for metaMIC results
  --label LABEL         Misassembly label of contigs for training assemblies
  -a ASSEMBLER, --assembler ASSEMBLER
                        The name of the directory of the trained model.
  -t THREADS, --threads THREADS
                        Maximum number of CPUs [default: 8]                        
```
