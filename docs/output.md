# Outputs of metaMIC

The output folder will contain

- **metaMIC_contig_score.txt**: metaMIC contig score for each contig (represent the likelihood of the misassemblies exist in the contig, *only for metagenome mode*).
- **misassembly_breakpoint.txt**: predicted misassembly breakpoints for misassembled contigs
- **anomaly_score.txt**: anomaly scores for each position in contigs.
- **metaMIC_corrected_contigs.fa**: fasta file of updated contigs
- **feature_matrix**: directory of generated feature matrix of inputs contigs.
- **temp**: some intermediate files



### metaMIC_contig_score.txt

contig | metaMIC_contig_score | length
---|--- | ---
k141_100014 | 0.25335 | 5066
k141_100173 | 0.71715 | 19223

### misassembly_breakpoint.txt

contig |misassembly_breakpoint | read_breakpoint_count | read_breakpoint_ratio | anomaly_score | anomaly_thred | metaMIC_contig_score | contig_length
---|--- |--- |--- |--- |--- |--- |--- 
k141_2323 | 6234 | 438 | 0.88 | 0.99 | 0.89 | 0.92585 | 8327
k141_87949 | 12968 | 360 | 0.80 | 1 | 0.89 | 0.850 | 9817

### anomaly_score.txt
contig | start_pos | anomaly_score | anomaly_thred | read_breakpoint_ratio |
---|--- | --- | --- | --- 
k141_108860 | 300 | 0.563 | 0.89 | 0.02
k141_108860 | 400 | 0.616 | 0.89 |0.009
