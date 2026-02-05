# Biopython
Computational Detection and Homology-Based Functional Annotation of Open Reading Frames in the Nipah Virus Genome

Project Overview

This project focuses on identifying and analyzing potential protein-coding regions in the complete genome of Nipah virus using Python and Biopython. The main objective is to detect Open Reading Frames (ORFs), select the most biologically significant ones, and predict their possible function through sequence similarity analysis.
Instead of manually copying gene sequences from databases, the entire process is performed computationally. The genome sequence is downloaded from NCBI in FASTA format and analyzed locally. This approach reflects how bioinformatics pipelines are used in real research environments to study viral genomes.

Study Framework

The script begins by reading the complete Nipah virus genome from a FASTA file using the SeqIO module from Biopython. SeqIO is responsible for reading biological sequence files and converting them into Python objects that can be manipulated programmatically. Once the file is read, the nucleotide sequence is extracted and stored for further analysis.
The total length of the genome is calculated using Python’s built-in len() function. After that, the GC content is determined using gc_fraction() from Bio.SeqUtils. GC content represents the percentage of guanine (G) and cytosine (C) nucleotides in the genome. This gives an idea of nucleotide composition and basic genomic characteristics.

ORF Identification (How Coding Regions Are Found)

The project then identifies Open Reading Frames across the three forward reading frames of the genome. Since DNA can be read in three different frames depending on the starting nucleotide, the script loops through frame positions 0, 1, and 2 using a for loop with range(3).
For each frame, the nucleotide sequence is translated into an amino acid sequence using the .translate() function. Translation converts every three nucleotides (a codon) into one amino acid. The argument to_stop=False ensures that the translation continues even after encountering stop codons.
Stop codons appear as "asterisk" in translated protein sequences. By splitting the translated sequence at "asterisk", the script separates potential protein fragments. Each fragment between two stop codons represents a possible ORF.
To avoid analyzing very small fragments that are unlikely to be functional proteins, the script filters ORFs by length. Only fragments with 150 or more amino acids are retained. These are stored in a list for further analysis.

Selecting Major ORFs

Once all ORFs are collected, they are sorted by length using Python’s built-in sorted() function. The sorting is done in descending order (reverse=True), so the longest ORFs appear first.
From this sorted list, the top four longest ORFs are selected. This step ensures that only the most significant coding regions are analyzed in detail. In viral genomes, longer ORFs often correspond to structural or replication-related proteins.

Protein Property Analysis

For each selected ORF, the script calculates basic protein characteristics. The ProteinAnalysis class from Bio.SeqUtils.ProtParam is used to compute molecular weight. Molecular weight is calculated by summing the average mass of all amino acids in the sequence.
Each ORF is also saved as a separate FASTA file. This allows individual analysis or reuse of the predicted protein sequences outside the script.

Functional Annotation Using BLAST

To predict the function of each ORF, the script performs a BLASTp search using NCBIWWW.qblast(). BLAST (Basic Local Alignment Search Tool) compares the predicted protein sequence against known proteins in the NCBI database.
The database used in this project is "nr" (non-redundant protein database), which contains a wide range of curated protein sequences.
The BLAST results are saved in XML format. The XML file is then parsed using NCBIXML.read(). From the BLAST output, the script extracts:
1.The top matching protein (best hit)
2.The alignment length
3.The number of identical amino acids
4.The E-value (statistical significance)

The percent identity is calculated manually using the formula:
(Number of identical residues / Alignment length) × 100

This provides a quantitative measure of similarity between the predicted ORF and known proteins.

Output Files Generated

The project generates:
1. Individual FASTA files for each selected ORF
2. XML files containing BLAST results
3. A summary text file containing:
   3.1 ORF number
   3.2 Protein length
   3.3 Molecular weight
   3.4 Top BLAST hit
   3.5 Percent identity
   3.6 E-value

This structured output makes the results reproducible and easy to review.

Significance of the Project
This project demonstrates how viral genome annotation can be performed computationally from raw sequence data. By combining ORF detection, biological filtering, physicochemical characterization, and homology-based functional prediction, the workflow reflects the logical structure of real-world bioinformatics pipelines.

Libraries and Tools Used:

1. Bio.SeqIO – Reading FASTA sequence files
2. Bio.SeqUtils.gc_fraction – Calculating GC content
3. Bio.SeqUtils.ProtParam.ProteinAnalysis – Computing molecular weight
4. Bio.Blast.NCBIWWW – Performing BLAST searches
5. Bio.Blast.NCBIXML – Parsing BLAST XML results
6. time – Managing request intervals
