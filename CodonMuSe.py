#!/usr/bin/env python
# -*- coding: utf-8 -*- 
"""
Input files must be fastas and in the format genus_species_formatted_CDS.fasta eg. Arabidopsis_thaliana_formatted_CDS.fasta
How to input a command:
python CodonMuSe.py -f Genus_species.fasta -tscan Genus_species_tRNAscan.txt -tc 1 -ind -fix_mb -par
If you have a tRNAscan you can include it as the second input of the command (as above). With or without the tRNAscan file you can also run
the script on individual genes (include the command -ind) with or without fixing the mutation_bias of each gene to the global mutation bias (from the formatted_CDS file) by adding -fix_mb.
"""
import glob
import sys
import math
import time
from scipy import optimize
import numpy as np
import random
import string
import os
import re
from scipy import stats

command_line = sys.argv
if "-f" in sys.argv:
	fasta_index = command_line.index("-f") + 1
	CDS_file = sys.argv[fasta_index]
	global species
	species, extention = os.path.splitext(CDS_file)
else:
	print "\n\t\tCodonMusel version 0.1.0\n"
	print "This software is distributed under the Univeristy of Oxford Academic Use\nLicence. For details please see the License.md that came with this software.\n"
	print "Usage:\npython CodonMuSe.py [options] -f <sequence file> -tscan <tRNAscan file>\n"	
	print "Options:"
	print "  -f <FILE>     A FASTA file of protein coding nucleotide seqeunces"
	print "  -tscan <FILE> A tRNA copy number file produced by tRNAscan"
	print "  -tc <INT>     The NCBI genetic code identifier goo.gl/ByQOau (Default = 1)"
	print "  -ind          Analysed individual genes in adition to a genomewide analysis"
	print "  -fix_mb       Fix mutation bias to genome-wide value for individual genes"
	print "  -par          Determine cost and efficiency optimality of individual genes\n"
	print "Citation:";
	print "When publishing work that uses CodonMuSe please cite:";
	print "Seward EA and Kelly S (2016) Genome Biology 17(1):226\n";
	exit()
if "-tscan" in sys.argv:
	tscan_index = command_line.index("-tscan") + 1
	tSCAN_file = sys.argv[tscan_index]
else:
	print "You have not specified the tRNAscan file using the format -tscan Genus_species_tRNAscan.txt therefore this analysis will be run without selection acting on translational effiency as an option\n"
if "-tc" in sys.argv:
	codon_code_number_index = command_line.index("-tc") + 1
	codon_code_number = int(sys.argv[codon_code_number_index])
else:
	print "\nYou haven't specified a translation table from 1-6 so the standard code (1) will be used\n";
	codon_code_number = 1
complementary_nuc = { 'A' : 'T', 'T' : 'A', 'C' :'G', 'G' : 'C'}
if codon_code_number == 1: # The standard code
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TGA', 'TAG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC', 'ATA'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'B', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
elif codon_code_number == 2: # The vertebrate mitochondrial code
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TAG', 'AGA', 'AGG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG', 'ATA'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG' , 'TGA'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'W', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'M', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'B', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'B' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
elif codon_code_number == 3: # The Yeast Mitochondrial Code 
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TAG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG'], 'M' : ['ATG', 'ATA'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGG', 'AGA', 'AGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG', 'CTT', 'CTC', 'CTA', 'CTG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG', 'TGA', ], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'W', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'T', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'T', 'CCC' : 'P', 'CAC' : 'H', 'CTA' : 'T', 'CCA' : 'P', 'CAA' : 'Q', 'CTG' : 'T', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'M', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
elif codon_code_number == 4: # The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TAG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC', 'ATA'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG', 'TGA'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'W', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
elif codon_code_number == 5: # The Invertebrate Mitochondrial Code
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TAG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG', 'ATA'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'AGA', 'AGG'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG', 'TGA'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'W', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'M', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'S', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'S', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
elif codon_code_number == 6: # The Ciliate, Dasycladacean and Hexamita Nuclear Code 
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TGA'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC', 'ATA'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG', 'TAA', 'TAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'Q', 'TGA' : 'B', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'Q', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
alphabetic_codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
symbols = {'Ala': 'A', 'Cys' : 'C', 'Asp' : 'D', 'Glu':'E', 'Phe':'F', 'Gly':'G', 'His': 'H', 'Ile': 'I', 'Lys':'K' , 'Leu': 'L', 'Met': 'M', 'Asn':'N', 'Pro':'P', 'Gln':'Q', 'Arg': 'R', 'Ser':'S' , 'Thr':'T',  'Val':'V', 'Trp':'W', 'Tyr':'Y'}
sij_values = {'AT' : 0, 'CG' : 0, 'GC' : 0, 'TA' : 0, 'TG' : 0.698, 'GT' : 0.6294, 'AG' : 1, 'CC' : 1, 'CT' : 1, 'GG' : 1, 'GA' : 1, 'TT' : 0.7, 'TC' : 0.95} #('AA' : 0.8773, 'CA' : 0.7309, 'AC' : 0.4211,)addressed in get_tAI  2014 Tuller et. al.#anticodon #codon
GC_content = {'A' : 0, 'T' : 0, 'C' : 1, 'G' : 1}
N_content = {'A' : 5, 'T' : 2, 'C' : 3, 'G' : 5} #N_content = {'A' : 1, 'T' : 0.4, 'C' : 0.6, 'G' : 1} #
En_content = {'A' : 0.899, 'T' : 0.571, 'C' : 0.671, 'G' : 0.867} #energy requirement of each nucleotide being synthesised in glucose molecules.
bases = ['A', 'T', 'C', 'G']

def Genome_wide_analysis(CDS_file):
	f1=open(species+"_GenomeWideResults.txt", "w")
	f2=open(species+"_excluded_sequences.txt", "w")				
	global codon_count
	codon_count = {}
	protein_count = {}
	for key in codon_trans_standard:
		codon_count[key] = 0
	for key in proteins:
		protein_count[key] = 0;
	codon = ""
	sequence_count = 0
	bad_sequence_count = 0
	print "Reading sequence file"
	with open(CDS_file) as CDSfile:
		for line in CDSfile:
			line = str(line.rstrip().strip())
			if line == "":
				empty = 0
			elif line[0] == ">":
				acc = line[1:]
			else:
				line = line.upper()
				temp = line.replace("A", "")
				temp = temp.replace("T", "")
				temp = temp.replace("C", "")
				temp = temp.replace("G", "")
				if len(temp) > 0:
					f2.write(acc+" contains "+temp+" ie. letters that aren't ATCG\n")
					bad_sequence_count = bad_sequence_count + 1
				elif codon_trans_standard[line[0:3]] == "M" and codon_trans_standard[line[len(line) - 3:]] == "B" and len(line)%3 == 0 and len(line)>30:
					sequence_count = sequence_count + 1
					p = 1
					for letter in line:
						if p % 3 == 0:
							p += 1
							codon = codon + letter
							codon_count[codon] = codon_count[codon] + 1
							protein_count[codon_trans_standard[codon]] = protein_count[codon_trans_standard[codon]] + 1 #print "%d %s %d %s" %(codon_count[codon], codon, protein_count[codon_trans_standard[codon]], codon_trans_standard[codon])
							codon = ""
						else:
							p += 1
							codon = codon + letter
				else:
					if len(line)<=30:
						f2.write(acc+" is "+len(line)+"bp long ie.< 30bp\n")
					elif len(line)%3 != 0:
						f2.write(acc+" is not divisible by 3 (possible frame shift error)\n")
					elif codon_trans_standard[line[len(line) - 3:]] != "B":
						f2.write(acc+" has no stop codon\n")
					elif codon_trans_standard[line[0:3]] != "M":
						f2.write(acc+" has no start codon\n")
					else:
						f2.write(acc+" has been excluded for mysterious reasons. Check why manually.\n")
					bad_sequence_count = bad_sequence_count + 1
	global relative_codon_use
	relative_codon_use = {}
	total_codons = 0
	for codon in codon_trans_standard:
		total_codons = total_codons + codon_count[codon]
		if codon_count[codon] != 0:
			relative_codon_use[codon] = "%.4f" %(float(codon_count[codon])/protein_count[codon_trans_standard[codon]])
		elif protein_count[codon_trans_standard[codon]] !=0:
			relative_codon_use[codon] = 0
			f2.write(codon+" isn't used once in the whole fasta file, is this odd?\n")
	log_likelihood = float(0)
	for codon in relative_codon_use:
		probability = float(relative_codon_use[codon])
		if probability == 0:
			probability = 0.000000000001
		probability = math.log(probability)
		log_likelihood = log_likelihood + probability*codon_count[codon]
	final_log_likelihood = "%.2f" %(log_likelihood)
	#Mutation_bias = Mb	Nitrogen Selection = Sc		Translational Efficiency = St		Energy selection = Es
	try:
		tRNA_count, bad_tRNAs = get_tAI_values(tSCAN_file)
		Options = ['Mb', 'Sc', 'St']
	except NameError:
		Options = ['Mb', 'Sc']
	print "Analysing whole genome sequence data (<1 min)"
	#AIC = 2k - 2ln(L) # minimum AIC is best
	best_log_likelihood = 10000000000000000000
	AIC = 10000000000000000000
	best = ''
	i = 0
	while i < len(Options):
		opt = Options[i]
		model = opt
		number_para = 1
		res_model, log_likelihood_model, r2_model = get_math_function(model, 0, 0, 0, 0)
		AIC_model = (2*number_para)-(2*(log_likelihood_model))
		if AIC_model < AIC:
			best_res, best_log_likelihood, best_r2, best, AIC = res_model, log_likelihood_model, r2_model, model, AIC_model
		j = i + 1
		while j < len(Options):
			opt_2 = Options[j]
			model = opt+"_"+opt_2
			number_para = 2
			res_model, log_likelihood_model, r2_model = get_math_function(model, 0, 0, 0, 0)
			AIC_model = (2*number_para)-(2*(log_likelihood_model))
			if AIC_model < AIC:
				best_res, best_log_likelihood, best_r2, best, AIC = res_model, log_likelihood_model, r2_model, model, AIC_model
			k = j + 1
			while k < len(Options):
				opt_3 = Options[k]
				model = opt+"_"+opt_2+"_"+opt_3
				number_para = 3
				res_model, log_likelihood_model, r2_model = get_math_function(model, 0, 0, 0, 0)
				AIC_model = (2*number_para)-(2*(log_likelihood_model))
				if AIC_model < AIC:
					best_res, best_log_likelihood, best_r2, best, AIC = res_model, log_likelihood_model, r2_model, model, AIC_model
				k = k + 1
			j = j + 1
		i = i + 1
	best_r2 = round(best_r2, 4)
	AIC = round(AIC, 1)
	best_log_likelihood = round(best_log_likelihood, 1)
	if ('Mb' in best):
		best_res[0] = round(best_res[0], 4)
	else:
		best_res[0] = 0
	if ('Sc' in best):
		best_res[1] = round(best_res[1], 4)
	else:
		best_res[1] = 0
	if ('St' in best):
		best_res[2] = round(best_res[2], 4)
	else:
		best_res[2] = 0
	model_used = best.replace("_", "+")
	f1.write("Genome-wide Results Summary:\nLn_L = "+str(best_log_likelihood)+"\nAIC = "+str(AIC)+"\nR2 = "+str(best_r2)+"\nMutation bias (Mb) = "+str(best_res[0])+"\nSelection on cost (Sc) = "+str(best_res[1])+"\nSelection on translational efficiency (St) = "+str(best_res[2])+"\nBest fitting model = "+model_used+"\n\n")
	f1.write("Amino acid\tCodon\tObserved frequency\tFitted frequency\n")
	for aa in proteins:
		for codon in proteins[aa]:
			fitted = get_probability_codon(codon, best, best_res[0], best_res[1], best_res[2], 0)
			fitted = round(fitted, 2)
			relative = float(relative_codon_use[codon])
			relative = round(relative, 2)
			f1.write(aa+"\t"+codon+"\t"+str(relative)+"\t"+str(fitted)+"\n")
	f1.close()
	f2.close()
	print "\n%d Sequences containing %d codons were analysed\n%d Sequences excluded (details in *_excluded_sequences.txt)" %(sequence_count, total_codons, bad_sequence_count)
	print "\n%d tRNAs were used\n%d tRNAs excluded (details in *_tRNAscan_errors.txt)" %(tRNA_count, bad_tRNAs)
	print "Results Summary:\nLn_L = %d\nAIC =\t%d\nR2 =\t%s\nMb =\t%s\nSc =\t%s\nSt =\t%s\nBest fitting model =%s\n" %(best_log_likelihood, AIC, str(best_r2), str(best_res[0]), str(best_res[1]), str(best_res[2]), model_used)
	return best, best_res[0]
	
def per_gene_analysis(CDS_file, best_model):
	print "Analysing individual genes. This takes ~ 1 second per gene.\n"
	f1=open(species+"_IndividualGenesResults.txt", "w")
	f1.write("Accession\tLog_likelihood\tMb\tSc\tSt")
	with open(CDS_file) as CDSfile:
		try:
			get_tAI_values(tSCAN_file)
		except NameError:
			print "You really should have a tRNAscan file for better results\n"
		for line in CDSfile:
			line = str(line.rstrip().strip())
			if line == "":
				empty = 0
			elif line[0] == ">":
				acc = line[1:]
				global codon_count
				codon_count = {}
				protein_count = {}
				for key in codon_trans_standard:
					codon_count[key] = 0
				for key in proteins:
					protein_count[key] = 0;
				codon = ""
			else:
				line = line.upper()
				temp = line.replace("A", "")
				temp = temp.replace("T", "")
				temp = temp.replace("C", "")
				temp = temp.replace("G", "")
				if len(temp) == 0 and codon_trans_standard[line[0:3]] == "M" and codon_trans_standard[line[len(line) - 3:]] == "B" and len(line)%3 == 0 and len(line)>30:
					p = 1
					for letter in line:
						if p % 3 == 0:
							p += 1
							codon = codon + letter	
							codon_count[codon] = codon_count[codon] + 1
							protein_count[codon_trans_standard[codon]] = protein_count[codon_trans_standard[codon]] + 1 #print "%d %s %d %s" %(codon_count[codon], codon, protein_count[codon_trans_standard[codon]], codon_trans_standard[codon])
							codon = ""
						else:
							p += 1
							codon = codon + letter
					global relative_codon_use
					relative_codon_use = {}
					for codon in codon_trans_standard:
						if codon_count[codon] != 0:
							relative_codon_use[codon] = "%.4f" %(float(codon_count[codon])/protein_count[codon_trans_standard[codon]])
						elif protein_count[codon_trans_standard[codon]] !=0:
							relative_codon_use[codon] = 0
					likelihood = float(0)
					for codon in relative_codon_use:
						probability = float(relative_codon_use[codon])
						if probability == 0:
							probability = 0.000000000001
						probability = math.log(probability)
						likelihood = likelihood + probability*codon_count[codon]
					final_likelihood = "%.2f" %(likelihood)
					res_model, likelihood_model, r2_model = get_math_function(best_model, 0, 0, 0, 0)
					r2_model = round(r2_model, 4)
					likelihood_model = round(likelihood_model, 1)
					if ('Mb' in best_model):
						res_model[0] = round(res_model[0], 4)
						if fixed_Mb != 'moveable':
							res_model[0] = round(fixed_Mb, 4)
					else:
						res_model[0] = 0
					if ('Sc' in best_model):
						res_model[1] = round(res_model[1], 4)
					else:
						res_model[1] = 0
					if ('St' in best_model):
						res_model[2] = round(res_model[2], 4)
					else:
						res_model[2] = 0
					f1.write("\n"+acc+"\t"+str(likelihood_model)+"\t"+str(res_model[0])+"\t"+str(res_model[1])+"\t"+str(res_model[2]))
	f1.close()

def get_math_function(model_type, start_Mb, start_Ns, start_Te, start_Es):
	equati = write_equation(model_type)
	def equation_real(x):
		Mb, Sc, St, Es = x
		try:
			return eval(equati)
		except OverflowError:
			return 10000000000000000000
		except ZeroDivisionError:
			return 10000000000000000000
		except ValueError:
			return 10000000000000000000
	x0 = np.asarray((start_Mb, start_Ns, start_Te, start_Es))
	result_parameters = optimize.fmin(equation_real, x0, disp=False)
	result_log_likelihood = -(equation_real(result_parameters))
	r2_value = get_r2_value(result_parameters, model_type)
	return result_parameters, result_log_likelihood, r2_value

def get_r2_value(x, model_type):
	Mb = x[0]
	Sc= x[1]
	St = x[2] 
	Es = x[3]
	top_line = 0
	bottom_line = 0
	mean = 0
	number_codons_used = 0
	for codon in relative_codon_use:
		probability_codon = get_probability_codon(codon, model_type, Mb, Sc, St, Es)
		mean = mean + probability_codon
		number_codons_used = number_codons_used + 1
	mean = mean/number_codons_used
	for codon in relative_codon_use:
		probability_codon = get_probability_codon(codon, model_type, Mb, Sc, St, Es)
		relative_codon = float(relative_codon_use[codon])
		value = probability_codon - relative_codon
		value = value*value
		top_line = top_line + value
		second_value = probability_codon - mean #mean of y will always be 0.328125 as you have 64 codon options and 21 amino acids. 21/64 = 0.328125
		second_value = second_value*second_value
		bottom_line = bottom_line + second_value
	top_over_bottom = float(top_line)/bottom_line
	return (1-top_over_bottom)
	
def get_probability_codon(input_codon, model_type, Mb, Sc, St, Es):
	prot = codon_trans_standard[input_codon]
	largest_value = 0
	for syn in proteins[prot]:
		if  "St" in model_type:
			if tAI_value[syn] > largest_value:
				largest_value = tAI_value[syn]
	upper = {}
	for syn in proteins[prot]:
		upper[syn] = ""
	parameters = model_type.split('_')
	number_parameters = len(parameters)
	for parameter in parameters:
		number_parameters -= 1
		for syn in proteins[prot]:
			if "Mb_fixed" in model_type:
				additional = get_upper_equation(syn, parameter, Mb, Sc, St, Es, 'for_probability', fixed_Mb, 'off')
			else:
				additional = get_upper_equation(syn, parameter, Mb, Sc, St, Es, 'for_probability', 'not', 'off')
			upper[syn] = upper[syn]+additional
			if number_parameters == 0:
				upper[syn] = upper[syn][:len(upper[syn])-1]
	base = "("
	for syn in proteins[prot]:
		base = base+upper[syn]+" + "
	base = base[:len(base)-3]+")"
	if model_type == "St" and prot == 'B':
		return float(relative_codon_use[input_codon])
	else:
		to_be_logged = eval(upper[input_codon]+"/"+base)
		if to_be_logged <=0.00000001:
			return 0.00000001
		else:
			return to_be_logged

def get_upper_equation(codon, parameter, Mb, Sc, St, Es, equation_type, fixed_Mb, shuffle_type):
	if equation_type == "for_probability":
		Mb = '%.6f'%Mb
		Sc = '%.6f'%Sc
		St = '%.6f'%St
		Es = '%.6f'%Es
	if parameter == 'Mb':
		if shuffle_type == 'off':
			GC_codon = '%.4f'%(GC_unshuffled[codon])
		else:
			GC_codon = '%.4f'%(GC_shuffled[codon])
		if fixed_Mb == 'not':
			return "(math.exp("+Mb+"*"+GC_codon+"))*"
		else:
			fixed_Mb = '%.8f'%fixed_Mb
			return "(math.exp("+fixed_Mb+"*"+GC_codon+"))*"
	elif parameter == 'Sc':
		if shuffle_type == 'off':
			Nitrogen_selection = '%.4f'%(Nitrogen_selection_unshuffled[codon])
		else:
			Nitrogen_selection = '%.4f'%(Nitrogen_selection_shuffled[codon])
		return "(math.exp("+Sc+"*"+Nitrogen_selection+"))*"
	elif parameter == 'St':
		if shuffle_type == 'off':
			val = '%.4f'%(tAI_value[codon] - largest_value[codon_trans_standard[codon]])
		else:
			val = '%.4f'%(tAI_shuffled[codon] - largest_value[codon_trans_standard[codon]])
		return "(math.exp("+St+"*"+val+"))*"
	elif parameter == 'Es':
		if shuffle_type == 'off':
			Energy = '%.4f'%(Energy_unshuffled[codon])
		else:
			Energy = '%.4f'%(Energy_shuffled[codon])
		return "(math.exp("+Es+"*"+Energy+"))*"
	else:
		print "What has gone wrong here? Not a recognised parameter %s\n" %parameter

def get_p_value(shuffle_type, real_likelihood, accuracy, Mb_in, Ns_in, Te_in): #
	counter = 1
	better_shuffled = 0
	while counter <= accuracy:
		shuf_equat = write_equation(shuffle_type)
		def equation_shuffle(x):
			Mb, Sc, St = x
			try:
				return eval(shuf_equat)
			except OverflowError:
				return 10000000000000000000
			except ZeroDivisionError:
				return 10000000000000000000
			except ValueError:
				return 10000000000000000000
		x0 = np.asarray((Mb_in, Ns_in, Te_in))
		res_shuf = optimize.fmin(equation_shuffle, x0, disp=False)
		likelihood = equation_shuffle(res_shuf)
		if likelihood <= real_likelihood:
			better_shuffled = better_shuffled + 1
		counter = counter + 1
	p_value = '%.2f'%(float(better_shuffled)/accuracy) #change this to %.3f if you change accuracy to 1000
	return p_value

def write_equation(model_type):
	global GC_unshuffled
	global Nitrogen_selection_unshuffled
	global Energy_unshuffled
	Energy_unshuffled = {}
	GC_unshuffled = {}
	Nitrogen_selection_unshuffled = {}
	for codon in alphabetic_codons:
		GC_unshuffled[codon] = GC_content[codon[0]] + GC_content[codon[1]] + GC_content[codon[2]]
		Nitrogen_selection_unshuffled[codon] = N_content[codon[0]] + N_content[codon[1]] + N_content[codon[2]]
		Energy_unshuffled[codon] = En_content[codon[0]] + En_content[codon[1]] + En_content[codon[2]]	
	if fixed_Mb == 'moveable':
		fixed_status = 'not'
	else:
		fixed_status = fixed_Mb
	if 'shuffle' in model_type:
		global tAI_shuffled
		global GC_shuffled
		global Nitrogen_selection_shuffled
		global Energy_shuffled
		if len(sys.argv) > 2 and "tRNAscan" in sys.argv[2]:
			tAI_shuffled = using_shuffle(tAI_value)
		GC_shuffled = using_shuffle(GC_unshuffled)
		Nitrogen_selection_shuffled = using_shuffle(Nitrogen_selection_unshuffled)
		Energy_shuffled = using_shuffle(Energy_unshuffled)
	str = "-("
	if  "St" in model_type:
		global largest_value
		largest_value = {}
		for pro in proteins:
			largest_value[pro] = 0
		for pro in largest_value:
			for syn in proteins[pro]:
				if "shuffle_Te" in model_type:
					if tAI_shuffled[syn] > largest_value[pro]:
						largest_value[pro] = tAI_shuffled[syn]
				else:
					if tAI_value[syn] > largest_value[pro]:
						largest_value[pro] = tAI_value[syn]
	upper = {}
	for codon in alphabetic_codons:
		upper[codon] = ""
	parameters = model_type.split('_')
	number_parameters = len(parameters)
	shuffle = 'off'
	for parameter in parameters:
		number_parameters -= 1
		if parameter != 'fixed':
			if parameter == 'shuffle':
				shuffle = 'on'
			else:
				for codon in alphabetic_codons:
					if shuffle == 'off':
						additional = get_upper_equation(codon, parameter, 'Mb', 'Sc', 'St',  'Es', 'for_real', fixed_status, 'off')
					else:
						additional = get_upper_equation(codon, parameter, 'Mb', 'Sc', 'St', 'Es', 'for_shuffle', fixed_status, 'on')
					upper[codon] = upper[codon]+additional
					if number_parameters == 0:
						upper[codon] = upper[codon][:len(upper[codon])-1]
				shuffle = 'off'
	for codon in alphabetic_codons:
		prot = codon_trans_standard[codon]
		base = "("
		for syn in proteins[prot]:
			base = base+upper[syn]+" + "
		base = base[:len(base)-3]+")"
		count_codon = '%i'%codon_count[codon]
		if model_type == "St" and prot == 'B':
			str = str+" + "
		else:
			Mb = 0.1
			St = 0.1
			Sc = 0.1
			Es = 0.1
			to_be_logged = eval(upper[codon]+"/"+base)
			if to_be_logged <=0.00000001:
				str = str+count_codon+"*(math.log(0.00000001)) + "
			elif to_be_logged >= 0.9999999:
				str = str+" + "
			else:
				str = str+count_codon+"*(math.log("+upper[codon]+"/"+base+")) + "
	str = str[:len(str)-3]+")"
	return str
	
def get_tAI_values(input):
	f3=open(species+"_tRNAscan_errors.txt", "w")
	anticodons = {}
	perfect_match = {}
	count = 0
	protein_index = {}
	anti_codon_index = {}
	tRNA_count = 0
	bad_tRNAs = 0
	for prot in proteins:
		protein_index[prot] = count
		anticodons[prot] = []
		if prot == "I":
			anticodons["I"].append("CAT") # This is because bacteria have an Ile-tRNA with CAT anticodon that doesn't match to ATG but wobble base pairs with ATA
		count = count + 1
		for syn in proteins[prot]:
			perfect_match[syn] = ''
			anti = complementary_nuc[syn[2]]+complementary_nuc[syn[1]]+complementary_nuc[syn[0]]
			perfect_match[syn] = anti
			anticodons[prot].append(anti);
		anti_codon_index_count = 0
		for anti in anticodons[prot]:
			anti_codon_index[anti] = anti_codon_index_count
			anti_codon_index_count = anti_codon_index_count + 1	
	anticodon_count = [[0 for antic in anticodons[p]] for p in proteins] #This assigns a matrix full of 0 for each protein and anticodon combination
	with open(input) as tRNAscan:
		for line in tRNAscan:
			line = str(line.rstrip().strip())
			bits = line.split()
			an = bits[5]
			if "tRNA" in line or "---" in line:
				line = line
			elif "SeC" in line or "Undet" in line or an == '???' or "Pseudo" in line or bits[4] == 'Sup' or 'fMet' in line:
				bad_tRNAs = bad_tRNAs + 1
				if "SeC" in line:
					f3.write(line+" includes SeC so is excluded from the analysis.\n")
				elif "Undet" in line:
					f3.write(line+" includes Undet so is excluded from the analysis\n")
				elif an == '???':
					f3.write(line+" includes ??? so is excluded from the analysis.\n")
				elif "Pseudo" in line:
					f3.write(line+" includes Pseudo so is excluded from the analysis\n")
				elif bits[4] == 'Sup':
					f3.write(line+" includes Sup so is excluded from the analysis\n")
				elif 'fMet' in line:
					f3.write(line+" includes fMet so is excluded from the analysis\n")
				else:
					f3.write(line+" excluded for unknown reason\n")
			else:
				tRNA_count = tRNA_count + 1
				amino_a = symbols[bits[4]]
				p_index = protein_index[amino_a]
				ac_index = anti_codon_index[an]
				#print amino_a, an, p_index, ac_index, anticodons[amino_a], anticodon_count[p_index][ac_index]
				if complementary_nuc[an[2]]+complementary_nuc[an[1]]+complementary_nuc[an[0]] in codon_trans_standard and an in anticodons[amino_a] and bits[8]>45:
					anticodon_count[p_index][ac_index] = anticodon_count[p_index][ac_index] + 1
				else:
					print line
	#######
	p_index = protein_index["I"]
	ac_index = anti_codon_index["CAT"]
	if anticodon_count[p_index][ac_index] == 0:
		anticodon_count[p_index][ac_index] = 1
		p_index = protein_index["M"]
		if anticodon_count[p_index][ac_index] >= 2:
			anticodon_count[p_index][ac_index] = anticodon_count[p_index][ac_index] - 1
	##### This bit is here as the current version of tRNA scan is not able to distinguish Met-tRNA(CAT) from Ile-tRNA(CAT).
	##### It would also be good to remove the fMet-tRNAs from the Met count but can't do that at the moment.
	sij_final = {}
	global tAI_value
	tAI_value = {}
	for codon in codon_trans_standard:
		tai = 0
		translated = codon_trans_standard[codon]
		p_index = protein_index[translated]
		for anticodon in anticodons[translated]:
			ac_index = anti_codon_index[anticodon]
			ideal = perfect_match[codon]
			if anticodon == ideal:
				sij_value = 0
			else:
				if ideal[1] == anticodon[1] and ideal[2] == anticodon[2]: #only wobble at position 3
					pair = anticodon[0]+complementary_nuc[ideal[0]]
					if codon_trans_standard[codon] == "R":
						sij_values["AC"] = 0.4211
						sij_values["AA"] = 0.8773
					else:
						sij_values["AC"] = 1
						sij_values["AA"] = 1
					if codon_trans_standard[codon] == "I":
						sij_values["CA"] = 0.7309
					else:
						sij_values["CA"] = 1
					sij_value = sij_values[pair]
				else:
					sij_value = 1
			tai = tai + (1- sij_value)*anticodon_count[p_index][ac_index]
		tAI_value[codon] = tai
	f3.close()
	return tRNA_count, bad_tRNAs

def using_shuffle(x):
	keys = x.keys()
	values = x.values()
	random.shuffle(values)
	return dict(zip(keys, values))

def pareto_frontier(Xs, Ys, maxX, maxY):
	myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
	#print myList
	p_front = [myList[0]]
	for pair in myList[1:]:
		if maxY: 
			if pair[1] >= p_front[-1][1]: # Look for higher values of Y…
				if pair[0] == p_front[-1][0]:
					p_front.pop()
				p_front.append(pair) # … and add them to the Pareto frontier
		else:
			if pair[1] <= p_front[-1][1]: # Look for lower values of Y…
				if pair[0] == p_front[-1][0]:
					p_front.pop()
				p_front.append(pair) # … and add them to the Pareto frontier
# Turn resulting pairs back into a list of Xs and Ys
	p_frontX = [pair[0] for pair in p_front]
	p_frontY = [pair[1] for pair in p_front]
	return p_frontX, p_frontY

def run_pareto_optimisation():
	print "Determining the optimality of gene sequences. This takes ~5 seconds per gene but varies with length.\n"
	print "When using this part of the software please also cite:\nSeward EA, Kelly S (DATE), Title, Journal\n"
	get_tAI_values(tSCAN_file)
	f1=open(species+"_OptimisationResults.txt", "w")
	f1.write("Accession\t%Both_optimised\t%Cost_optimised\t%tAI_optimised\n")
	with open(CDS_file) as CDSfile:
		for sequence in CDSfile:
			sequence = str(sequence.rstrip().strip())
			if sequence == "":
				empty = 0
			elif sequence[0] == ">":
				acc = sequence[1:]
			else:
				sequence = sequence.upper()
				temp = sequence.replace("A", "")
				temp = temp.replace("T", "")
				temp = temp.replace("C", "")
				temp = temp.replace("G", "")
				if len(temp) == 0 and codon_trans_standard[sequence[0:3]] == "M" and codon_trans_standard[sequence[len(sequence) - 3:]] == "B" and len(sequence)%3 == 0 and len(sequence)>30:
					codon = ""
					amino_seq = []
					cost = 0
					translatability = 0
					line = {}
					number = 0
					options = []
					column = 0
					p = 1
					for letter in sequence:
						if p % 3 == 0:
							p += 1
							codon = codon + letter
							if codon_trans_standard[codon] != 'B':
								amino_seq.append(codon_trans_standard[codon])
								column += 1
								cost = cost + N_content[codon[0]] + N_content[codon[1]] + N_content[codon[2]]
								translatability = translatability + tAI_value[codon]
								line[number] = proteins[codon_trans_standard[codon]]
								options.append(len(line[number]))
								number += 1
							codon = ""
						else:
							p += 1
							codon = codon + letter
					i = 0
					cost_frontier = []
					trans_frontier = []
					while (i < len(line)):
						opts = line[i]
						cost_array = []
						trans_array = []
						for codon in opts:
							cost_array.append(N_content[codon[0]] + N_content[codon[1]] + N_content[codon[2]])
							trans_array.append(tAI_value[codon])
						temp_frontier_cost = []
						temp_frontier_trans = []
						if len(cost_frontier) == 0:
							temp_frontier_cost = cost_array
							temp_frontier_trans = trans_array
						else:
							for costy in cost_frontier:
								for extra_cost in cost_array:
									temp_frontier_cost.append(costy+extra_cost)
							for trans in trans_frontier:
								for extra_trans in trans_array:
									temp_frontier_trans.append(trans+extra_trans)
						cost_frontier, trans_frontier = pareto_frontier(temp_frontier_cost, temp_frontier_trans, maxX = False, maxY = True)
						i += 1
					i = 0
					cost_frontier = np.array(cost_frontier)
					trans_frontier = np.array(trans_frontier)
					min_cost = min(cost_frontier)
					max_cost = max(cost_frontier)
					min_trans = min(trans_frontier)
					max_trans = max(trans_frontier)
					cost_frontier_worst = []
					trans_frontier_worst = []
					while (i < len(line)):
						opts = line[i]
						cost_array = []
						trans_array = []
						for codon in opts:
							cost_array.append(N_content[codon[0]] + N_content[codon[1]] + N_content[codon[2]])
							trans_array.append(tAI_value[codon])
						temp_frontier_cost = []
						temp_frontier_trans = []
						if len(cost_frontier_worst) == 0:
							temp_frontier_cost = cost_array
							temp_frontier_trans = trans_array
						else:
							for costy in cost_frontier_worst:
								for extra_cost in cost_array:
									temp_frontier_cost.append(costy+extra_cost)
							for trans in trans_frontier_worst:
								for extra_trans in trans_array:
									temp_frontier_trans.append(trans+extra_trans)
						cost_frontier_worst, trans_frontier_worst = pareto_frontier(temp_frontier_cost, temp_frontier_trans, maxX = True, maxY = False)
						i += 1
					cost_frontier_worst = np.array(cost_frontier_worst)
					trans_frontier_worst = np.array(trans_frontier_worst)
					if min_cost > min(cost_frontier_worst):
						min_cost = min(cost_frontier_worst)
					if max_cost < max(cost_frontier_worst):
						max_cost = max(cost_frontier_worst)
					if min_trans > min(trans_frontier_worst):
						min_trans = min(trans_frontier_worst)
					if max_trans < max(trans_frontier_worst):
						max_trans = max(trans_frontier_worst)
					cost_frontier_worst = (cost_frontier_worst - min_cost)/float(max_cost - min_cost)
					trans_frontier_worst = (trans_frontier_worst - min_trans)/float(max_trans - min_trans)
					cost_frontier = (cost_frontier - min_cost)/float(max_cost - min_cost)
					trans_frontier = (trans_frontier - min_trans)/float(max_trans - min_trans)
					cost = float(cost - min_cost)/(max_cost - min_cost)
					translatability = float(translatability - min_trans)/(max_trans - min_trans)
					i = 0
					distance = []
					while (i < len(cost_frontier)):
						distance.append(math.sqrt((cost - cost_frontier[i])*(cost - cost_frontier[i]) + (translatability  - trans_frontier[i])*(translatability  - trans_frontier[i])))
						i += 1
					d1 = min(distance)
					d2 = (math.sqrt((cost - cost_frontier[0])*(cost - cost_frontier[0]) + (translatability - trans_frontier[0])*(translatability - trans_frontier[0])))
					d3 = (math.sqrt((cost - cost_frontier[-1])*(cost - cost_frontier[-1]) + (translatability - trans_frontier[-1])*(translatability - trans_frontier[-1])))
					i = 0
					i = 0
					distance = []
					while (i < len(cost_frontier_worst)):
						distance.append(math.sqrt((cost - cost_frontier_worst[i])*(cost - cost_frontier_worst[i]) + (translatability  - trans_frontier_worst[i])*(translatability  - trans_frontier_worst[i])))
						i += 1
					d4 = min(distance)
					d5 = (math.sqrt((cost - cost_frontier_worst[-1])*(cost - cost_frontier_worst[-1]) + (translatability - trans_frontier_worst[-1])*(translatability - trans_frontier_worst[-1])))
					d6 = (math.sqrt((cost - cost_frontier_worst[0])*(cost - cost_frontier_worst[0]) + (translatability - trans_frontier_worst[0])*(translatability - trans_frontier_worst[0])))
					Both_optimised = float(d4*100)/(d1 + d4)
					Cost_optimised = float(d5*100)/(d2 + d5)
					tAI_optimised = float(d6*100)/(d3 + d6)
					Both_optimised = round(Both_optimised, 2)
					Cost_optimised = round(Cost_optimised, 2)
					tAI_optimised = round(tAI_optimised, 2)
					f1.write(acc+"\t"+str(Both_optimised)+"\t"+str(Cost_optimised)+"\t"+str(tAI_optimised)+"\n")
	f1.close()

global fixed_Mb
fixed_Mb = "moveable"
model_to_use, fixed_Mb = Genome_wide_analysis(CDS_file)
if "-ind" in sys.argv:
	if "-fix_mb" in sys.argv:
		print "\nMutaiton bias is will be fixed to the genome-wide Mb value for individual gene analysis.\n"
	else:
		fixed_Mb = "moveable"
		print "\nMutaiton bias is not fixed for individual genes. I would not recommend this. add -fix_mb to end of command line to run with Mb values fixed to genome-wide Mb values\n"
	per_gene_analysis(CDS_file, model_to_use)
if "-par" in sys.argv:
	run_pareto_optimisation()
	
print "CodonMuSe (1) implements the SK model (2).\nWhen publishing work that uses CodonMuSe please cite both:\n(1) Seward EA and Kelly S (2017) bioRxiv doi.org/XX.XXXX/XXXXXX\n(2) Seward EA and Kelly S (2016) Genome Biology 17(1):226\n"