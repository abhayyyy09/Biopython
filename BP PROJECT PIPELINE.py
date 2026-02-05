# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 11:34:34 2026

@author: abhay
"""

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time

# STEP 1: Read Genome

record = SeqIO.read("nipah_genome.fasta", "fasta")
genome_sequence = record.seq

print("========================================================================")
print("Length of the Genome is:", len(genome_sequence))
gc_content = gc_fraction(genome_sequence) * 100
print("GC content of the Genome is:", round(gc_content, 2), "%")

# STEP 2: Identify ORFs for ORF Analysis (Forward 3 Frames)

orf_list = []

for frame in range(3):

    translated_sequence = genome_sequence[frame:].translate(to_stop=False)
    protein_fragments = str(translated_sequence).split("*")
    # as * means there is a stop codon there 

    for fragment in protein_fragments:
        if len(fragment) >= 150:
            orf_list.append(fragment)

print("========================================================================")
print("The total no. of ORFs detected which have more than 150 amino acids):", len(orf_list))
print("========================================================================")
# STEP 3: Ranking the ORFs by Length


sorted_orfs = sorted(orf_list, key=len, reverse=True)
# this sorts the top 4 orfs with most amino acids

top_orfs = []

count = 0
for protein in sorted_orfs:
    if count < 4:
        top_orfs.append(protein)
        count += 1

print("Number of ORFs selected is :", len(top_orfs))
print("========================================================================")

# STEP 4: Analyze Each Selected ORF

summary_file = open("orf_analysis_summary.txt", "w")
summary_file.write("Nipah Virus ORF Functional Analysis\n")
summary_file.write("====================================\n")

orf_number = 1

for protein in top_orfs:
    
    print("\nAnalyzing ORF", orf_number)
    
    # Protein length
    protein_length = len(protein)

    # Molecular weight
    analysis = ProteinAnalysis(protein)
    molecular_weight = analysis.molecular_weight()

    print("Length:", protein_length, "aa")
    print("Molecular Weight:", round(molecular_weight, 2), "Da")

    # Save FASTA file
    fasta_filename = "orf_" + str(orf_number) + ".fasta"
    fasta_file = open(fasta_filename, "w")
    fasta_file.write(">ORF_" + str(orf_number) + "\n")
    fasta_file.write(protein)
    fasta_file.close()


    # STEP 5: BLAST of different orfs

    print("Running BLAST for ORF", orf_number)

    result_handle = NCBIWWW.qblast(
        program="blastp",
        database="nr",
        sequence=protein
    )

    xml_filename = "orf_" + str(orf_number) + "_blast.xml"
    xml_file = open(xml_filename, "w")
    xml_file.write(result_handle.read())
    xml_file.close()


    # STEP 6: Analysis of the BLAST Result
    
    result_handle = open(xml_filename)
    blast_record = NCBIXML.read(result_handle)

    if len(blast_record.alignments) > 0:

        alignment = blast_record.alignments[0]
        hsp = alignment.hsps[0]

        identities = hsp.identities
        alignment_length = hsp.align_length
        percent_identity = (identities / alignment_length) * 100
        e_value = hsp.expect
        
        print("Top Hit:", alignment.title)
        print("Percent Identity:", round(percent_identity, 2), "%")
        print("E-value:", e_value)
        print("========================================================================")
        
        # Write to summary file
        summary_file.write("\nORF " + str(orf_number) + "\n")
        summary_file.write("Length: " + str(protein_length) + " aa\n")
        summary_file.write("Molecular Weight: " + str(round(molecular_weight, 2)) + " Da\n")
        summary_file.write("Top Hit: " + alignment.title + "\n")
        summary_file.write("Percent Identity: " + str(round(percent_identity, 2)) + "%\n")
        summary_file.write("E-value: " + str(e_value) + "\n")

    time.sleep(5)

    orf_number += 1


summary_file.close()

print("\n Analysis completed successfully.")
