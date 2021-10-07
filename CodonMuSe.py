#!/usr/bin/env python
# -*- coding: utf-8 -*- 
"""
Input files must be fastas and in the format genus_species_formatted_CDS.fasta eg. Arabidopsis_thaliana_formatted_CDS.fasta
How to input a command:
python CodonMuSe.py -f Genus_species.fasta -tscan Genus_species_tRNAscan.txt -tc 1 -ind -GC gw -par
If you have a tRNAscan you can include it as the second input of the command (as above). With or without the tRNAscan file you can also run
the script on individual genes (include the command -ind) with or without fixing the GC content of each gene to the global GC content (from the formatted_CDS file) by adding -GC gw.
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
	print "\n\t\tCodonMuse version 0.1.0\n"
	print "This software is distributed under the Univeristy of Oxford Academic Use\nLicence. For details please see the License.md that came with this software.\n"
else:
	print "\n\t\tCodonMuse version 0.1.0\n"
	print "This software is distributed under the Univeristy of Oxford Academic Use\nLicence. For details please see the License.md that came with this software.\n"
	print "Usage:\npython CodonMuSe.py [options] -f <sequence file> -tscan <tRNAscan file>\n"	
	print "Options:"
	print "  -f <FILE>     A FASTA file of protein coding nucleotide sequences"
	print "  -tscan <FILE> A tRNA copy number file produced by tRNAscan"
	print "  -tc <INT>     The NCBI genetic code identifier goo.gl/ByQOau (Default = 1)"
	print "  -ind          Analysed individual genes in addition to a genomewide analysis"
	print "  -par          Determine cost and efficiency optimality of individual genes"
	print "  -trade        Calculates trade-off between cost and efficiency"
	print "  -GC <float>   Specify a fixed value for GC in all calculations"
	print "  -Sc <float>   Specify a fixed value for Sc in all calculations"
	print "  -St <float>   Specify a fixed value for Sc in all calculations"
	print "  -GC gw        Fix GC to genome-wide value for individual genes"
	print "  -Sc gw        Fix Sc to genome-wide value for individual genes"
	print "  -St gw        Fix St to genome-wide value for individual genes"
	print "  -m <TXT>      Specify model parameters (GC, Sc, St) eg. GC_Sc_St or GC_Sc"
	print "                (Default = determine best parameter combination automatically)\n"
	print "Citation:";
	print "CodonMuSe (1) implements the SK model (2)."
	print "When publishing work that uses CodonMuSe please cite both:"
	print "(1) Seward EA and Kelly S (2017) bioRxiv doi.org/10.1101/136861"
	print "(2) Seward EA and Kelly S (2016) Genome Biology 17(1):226\n"

	exit()
if "-tscan" in sys.argv:
	tscan_index = command_line.index("-tscan") + 1
	tSCAN_file = sys.argv[tscan_index]
else:
	print "You have not specified the tRNAscan file using the format -tscan Genus_species_tRNAscan.txt\nTherefore this analysis will run without using St (selection acting on translational effiency)\n"
if "-tc" in sys.argv:
	codon_code_number_index = command_line.index("-tc") + 1
	codon_code_number = int(sys.argv[codon_code_number_index])
else:
	print "\nNo translation table specified (1-11) so the standard code (1) will be used\n";
	codon_code_number = 1
if "-m" in sys.argv:
	model_index = command_line.index("-m") + 1
	model_fix = sys.argv[model_index]
else:
	model_fix = "find"
complementary_nuc = { 'A' : 'T', 'T' : 'A', 'C' :'G', 'G' : 'C'}
if codon_code_number == 1: # The standard code
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TGA', 'TAG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC', 'ATA'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'B', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
	start_codons = ['ATG', 'CTG', 'TTG']
elif codon_code_number == 2: # The vertebrate mitochondrial code
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TAG', 'AGA', 'AGG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG', 'ATA'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG' , 'TGA'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'W', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'M', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'B', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'B', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
	start_codons = ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
elif codon_code_number == 3: # The Yeast Mitochondrial Code 
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TAG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG'], 'M' : ['ATG', 'ATA'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGG', 'AGA', 'AGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG', 'CTT', 'CTC', 'CTA', 'CTG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG', 'TGA', ], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'W', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'T', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'T', 'CCC' : 'P', 'CAC' : 'H', 'CTA' : 'T', 'CCA' : 'P', 'CAA' : 'Q', 'CTG' : 'T', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'M', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
	start_codons = ['ATG', 'ATA']
elif codon_code_number == 4: # The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TAG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC', 'ATA'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG', 'TGA'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'W', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
	start_codons = ['TTA', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
elif codon_code_number == 5: # The Invertebrate Mitochondrial Code
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TAG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG', 'ATA'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'AGA', 'AGG'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG', 'TGA'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'W', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'M', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'S', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'S', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
	start_codons = ['TTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
elif codon_code_number == 6: # The Ciliate, Dasycladacean and Hexamita Nuclear Code 
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TGA'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC', 'ATA'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG', 'TAA', 'TAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'Q', 'TGA' : 'B', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'Q', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
	start_codons = ['ATG']	
elif codon_code_number == 9: # The Echinoderm and Flatworm Mitochondrial Code
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TAG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC', 'ATA'], 'K' : ['AAG'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG'], 'N' : ['AAT', 'AAC', 'AAA'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'AGG', 'AGA'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG', 'TGA'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'W', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'N', 'AGA' : 'S', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'S', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
	start_codons = ['ATG', 'GTG']
elif codon_code_number == 10: # The Euplotid Nuclear Code
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TAG'], 'C' : ['TGT', 'TGC', 'TGA'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC', 'ATA'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'C', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
	start_codons = ['ATG']
elif codon_code_number == 11: # The Bacterial, Archaeal and Plant Plastid Code
	proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TGA', 'TAG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC', 'ATA'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG'], 'Y' : ['TAT', 'TAC']}
	codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'B', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
	start_codons = ['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
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
	global cds_filtered
	cds = {}
	cds_filtered = {}
	print "Reading sequence file"
	with open(CDS_file) as CDSfile:
		for line in CDSfile:
			line = str(line.rstrip().strip())
			if line == "":
				empty = 0
			elif line[0] == ">":
				acc = line[1:]
				cds[acc] = []
			else:
				line = line.upper()
				cds[acc].append(line)
	for acc in cds:
		line = ''.join(cds[acc])
		temp = line.replace("A", "")
		temp = temp.replace("T", "")
		temp = temp.replace("C", "")
		temp = temp.replace("G", "")
		if len(temp) > 0:
			f2.write(acc+" contains "+temp+" ie. letters that aren't ATCG\n")
			bad_sequence_count = bad_sequence_count + 1
		elif line[0:3] in start_codons and codon_trans_standard[line[len(line) - 3:]] == "B" and len(line)%3 == 0 and len(line)>30:
			cds_filtered[acc] = line
			sequence_count = sequence_count + 1
			p = 1
			for letter in line:
				if p % 3 == 0:
					p += 1
					codon = codon + letter
					codon_count[codon] = codon_count[codon] + 1
					protein_count[codon_trans_standard[codon]] = protein_count[codon_trans_standard[codon]] + 1
					codon = ""
				else:
					p += 1
					codon = codon + letter
		else:
			if len(line)<=30:
				f2.write(acc+" is "+str(len(line))+"bp long ie.< 30bp\n")
			elif len(line)%3 != 0:
				f2.write(acc+" is not divisible by 3 (possible frame shift error)\n")
			elif codon_trans_standard[line[len(line) - 3:]] != "B":
				f2.write(acc+" has no stop codon\n")
			elif line[0:3] not in start_codons:
				f2.write(acc+"\t"+str(line[0:3])+" is not a start codon\n")
			else:
				f2.write(acc+" has been excluded for mysterious reasons. Check why manually.\n")
			bad_sequence_count = bad_sequence_count + 1
	del cds
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
		else:
			relative_codon_use[codon] = 0
	log_likelihood = float(0)
	for codon in relative_codon_use:
		probability = float(relative_codon_use[codon])
		if probability == 0:
			probability = 0.000000000001
		probability = math.log(probability)
		log_likelihood = log_likelihood + probability*codon_count[codon]
	final_log_likelihood = "%.2f" %(log_likelihood)
	#GC_content = GC	Nitrogen Selection = Sc		Translational Efficiency = St		Energy selection = Es
	try:
		tRNA_count, bad_tRNAs = get_tAI_values(tSCAN_file)
		Options = ['GC', 'Sc', 'St']
	except NameError:
		Options = ['GC', 'Sc']
	print "Analysing whole genome sequence data (<1 min)"
	if model_fix == "find":
		print "Finding best model to use for this data set\nModel\t\t-LnL\tdf\tAIC\tBIC"
		f1.write("Testing different model parameters for best fit:\nModel\t-LnL\tdf\tAIC\tBIC\n")
		#AIC = 2k - 2ln(L) # minimum AIC is best
		best_log_likelihood = 10000000000000000000
		AIC = 10000000000000000000
		log_number_data_points = math.log(total_codons)
		best = ''
		i = 0
		while i < len(Options):
			opt = Options[i]
			model = opt
			number_para = 1
			res_model, log_likelihood_model, r2_model = get_math_function(model, 0, 0, 0)
			AIC_model = (2*number_para)-(2*(log_likelihood_model))
			BIC = log_number_data_points*number_para - (2*log_likelihood_model)
			log_likelihood_model = round(log_likelihood_model, 1)
			AIC_model = round(AIC_model,1)
			BIC = round(BIC, 1)
			print "%s\t\t%s\t%s\t%s\t%s" %(model, -log_likelihood_model, number_para, AIC_model, BIC)
			f1.write(model+"\t\t\t"+str(-log_likelihood_model)+"\t"+str(number_para)+"\t"+str(AIC_model)+"\t"+str(BIC)+"\n")
			if AIC_model < AIC:
				best_res, best_log_likelihood, best_r2, best, AIC = res_model, log_likelihood_model, r2_model, model, AIC_model
			j = i + 1
			while j < len(Options):
				opt_2 = Options[j]
				model = opt+"_"+opt_2
				number_para = 2
				res_model, log_likelihood_model, r2_model = get_math_function(model, 0, 0, 0)
				AIC_model = (2*number_para)-(2*(log_likelihood_model))
				BIC = log_number_data_points*number_para - (2*log_likelihood_model)
				log_likelihood_model = round(log_likelihood_model, 1)
				AIC_model = round(AIC_model,1)
				BIC = round(BIC, 1)
				model_here = model.replace("_", "+")
				print "%s\t\t%s\t%s\t%s\t%s" %(model_here, -log_likelihood_model, number_para, AIC_model, BIC)
				f1.write(model_here+"\t\t"+str(-log_likelihood_model)+"\t"+str(number_para)+"\t"+str(AIC_model)+"\t"+str(BIC)+"\n")
				if AIC_model < AIC:
					best_res, best_log_likelihood, best_r2, best, AIC = res_model, log_likelihood_model, r2_model, model, AIC_model
				k = j + 1
				while k < len(Options):
					opt_3 = Options[k]
					model = opt+"_"+opt_2+"_"+opt_3
					number_para = 3
					res_model, log_likelihood_model, r2_model = get_math_function(model, 0, 0, 0)
					AIC_model = (2*number_para)-(2*(log_likelihood_model))
					BIC = log_number_data_points*number_para - (2*log_likelihood_model)
					log_likelihood_model = round(log_likelihood_model, 1)
					AIC_model = round(AIC_model,1)
					BIC = round(BIC, 1)
					model_here = model.replace("_", "+")
					print "%s\t%s\t%s\t%s\t%s" %(model_here, -log_likelihood_model, number_para, AIC_model, BIC)
					f1.write(model_here+"\t"+str(-log_likelihood_model)+"\t"+str(number_para)+"\t"+str(AIC_model)+"\t"+str(BIC)+"\n")
					if AIC_model < AIC:
						best_res, best_log_likelihood, best_r2, best, AIC = res_model, log_likelihood_model, r2_model, model, AIC_model
					k = k + 1
				j = j + 1
			i = i + 1
		model_used = best.replace("_", "+")
		print "%s is the best fitting model" %(model_used)
	else:
		if 'GC' in model_fix and 'Sc' in model_fix and 'St' in model_fix:
			number_para = 3
			best = 'GC_Sc_St'
		elif 'GC' in model_fix and 'Sc' in model_fix:
			number_para = 2
			best = 'GC_Sc'
		elif 'GC' in model_fix and 'St' in model_fix:
			number_para = 2
			best = 'GC_St'
		elif 'St' in model_fix and 'Sc' in model_fix:
			number_para = 2
			best = 'Sc_St'
		elif 'GC' in model_fix or 'Sc' in model_fix or 'St' in model_fix:
			number_para = 1
			best = model_fix.replace("_", "")
		else:
			print "Please input the specified model using the options GC Sc St separated by _\neg. GC_Sc_St or GC_St or GC etc."
		print "\nYou have specified model %s which has %d input(s).\n" %(best, number_para)
		best_res, best_log_likelihood, best_r2 = get_math_function(best, 0, 0, 0)
		AIC = (2*number_para)-(2*(best_log_likelihood))
	best_r2 = round(best_r2, 4)
	AIC = round(AIC, 1)
	best_log_likelihood = round(best_log_likelihood, 1)
	model_used = best.replace("_", "+")
	f1.write("\nGenome-wide Results Summary:\nLnL = "+str(best_log_likelihood)+"\nAIC = "+str(AIC)+"\nR2 = "+str(best_r2)+"\nGC content (GC) = "+str(best_res[0])+"\nSelection on cost (Sc) = "+str(best_res[1])+"\nSelection on translational efficiency (St) = "+str(best_res[2])+"\nModel used = "+model_used+"\n\n")
	f1.write("Amino acid\tCodon\tObserved frequency\tFitted frequency\n")
	for aa in proteins:
		if aa != 'B':
			for codon in proteins[aa]:
				fitted = get_probability_codon(codon, best, best_res[0], best_res[1], best_res[2])
				fitted = round(fitted, 2)
				relative = float(relative_codon_use[codon])
				relative = round(relative, 2)
				f1.write(aa+"\t"+codon+"\t"+str(relative)+"\t"+str(fitted)+"\n")
	f1.close()
	f2.close()
	print "\n%d Sequences containing %d codons were analysed\n%d Sequences excluded (details in *_excluded_sequences.txt)" %(sequence_count, total_codons, bad_sequence_count)
	print "\n%d tRNAs were used\n%d tRNAs excluded (details in *_tRNAscan_errors.txt)\n" %(tRNA_count, bad_tRNAs)
	print "Results Summary:\nLn_L = %d\nAIC =\t%d\nR2 =\t%s (excluding amino acids with only one codon)\nGC =\t%s\nSc =\t%s\nSt =\t%s\nModel = %s\n" %(best_log_likelihood, AIC, str(best_r2), str(best_res[0]), str(best_res[1]), str(best_res[2]), model_used)
	return best, best_res[0], best_res[1], best_res[2]
	
def per_gene_analysis(best_model):
	print "Analysing individual genes. This takes ~ 1 second per gene.\n"
	f1=open(species+"_IndividualGenesResults.txt", "w")
	f1.write("Accession\tLog_likelihood\tGC\tSc\tSt")
	try:
		get_tAI_values(tSCAN_file)
	except NameError:
		print "You really should have a tRNAscan file for better results\n"
	for acc in cds_filtered:
		global codon_count
		codon_count = {}
		protein_count = {}
		for key in codon_trans_standard:
			codon_count[key] = 0
		for key in proteins:
			protein_count[key] = 0;
		codon = ""
		line = cds_filtered[acc]
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
		res_model, likelihood_model, r2_model = get_math_function(best_model, 0, 0, 0)
		r2_model = round(r2_model, 4)
		likelihood_model = round(likelihood_model, 1)
		f1.write("\n"+acc+"\t"+str(likelihood_model)+"\t"+str(res_model[0])+"\t"+str(res_model[1])+"\t"+str(res_model[2]))
	f1.close()

def get_math_function(model_type, start_GC, start_Ns, start_Te):
	equati = write_equation(model_type, GC_to_use, Sc_to_use, St_to_use)
	def equation_real(x):
		GC, Sc, St = x
		try:
			return eval(equati)
		except OverflowError:
			return 10000000000000000000
		except ZeroDivisionError:
			return 10000000000000000000
		except ValueError:
			return 10000000000000000000
	x0 = np.asarray((start_GC, start_Ns, start_Te))
	result_parameters = optimize.fmin(equation_real, x0, disp=False)
	if 'GC' in model_type:
		if GC_to_use == 'GC':
			result_parameters[0] = round(result_parameters[0], 4)
		else:
			result_parameters[0] = GC_to_use
	else:
		result_parameters[0] = 0
	if 'Sc' in model_type:
		if Sc_to_use == 'Sc':
			result_parameters[1] = round(result_parameters[1], 4)
		else:
			result_parameters[1] = Sc_to_use
	else:
		result_parameters[1] = 0
	if 'St' in model_type:
		if St_to_use == 'St':
			result_parameters[2] = round(result_parameters[2], 4)
		else:
			result_parameters[2] = St_to_use
	else:
		result_parameters[2] = 0
	result_log_likelihood = -(equation_real(result_parameters))
	r2_value = get_r2_value(result_parameters, model_type)
	return result_parameters, result_log_likelihood, r2_value

def get_r2_value(x, model_type):
	GC = x[0]
	Sc= x[1]
	St = x[2]
	top_line = 0
	bottom_line = 0
	mean = 0
	number_codons_used = 0
	for codon in relative_codon_use:
		probability_codon = get_probability_codon(codon, model_type, GC, Sc, St)
		mean = mean + probability_codon
		number_codons_used = number_codons_used + 1
	mean = mean/number_codons_used
	for codon in relative_codon_use:
		probability_codon = get_probability_codon(codon, model_type, GC, Sc, St)
		relative_codon = float(relative_codon_use[codon])
		value = probability_codon - relative_codon
		value = value*value
		top_line = top_line + value
		second_value = probability_codon - mean #mean of y will always be 0.328125 as you have 64 codon options and 21 amino acids. 21/64 = 0.328125
		second_value = second_value*second_value
		bottom_line = bottom_line + second_value
	top_over_bottom = float(top_line)/bottom_line
	return (1-top_over_bottom)
	
def get_probability_codon(input_codon, model_type, GC, Sc, St):
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
			additional = get_upper_equation(syn, parameter, GC, Sc, St, 'for_probability', 'off')
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

def get_upper_equation(codon, parameter, GC, Sc, St, equation_type, shuffle_type):
	if equation_type == "for_probability":
		GC = '%.6f'%GC
		Sc = '%.6f'%Sc
		St = '%.6f'%St
	if parameter == 'GC':
		if shuffle_type == 'off':
			GC_codon = '%.4f'%(GC_unshuffled[codon])
		else:
			GC_codon = '%.4f'%(GC_shuffled[codon])
		return "(math.exp("+str(GC)+"*"+GC_codon+"))*"
	elif parameter == 'Sc':
		if shuffle_type == 'off':
			Nitrogen_selection = '%.4f'%(Nitrogen_selection_unshuffled[codon])
		else:
			Nitrogen_selection = '%.4f'%(Nitrogen_selection_shuffled[codon])
		return "(math.exp("+str(Sc)+"*"+Nitrogen_selection+"))*"
	elif parameter == 'St':
		if shuffle_type == 'off':
			val = '%.4f'%(tAI_value[codon] - largest_value[codon_trans_standard[codon]])
		else:
			val = '%.4f'%(tAI_shuffled[codon] - largest_value[codon_trans_standard[codon]])
		return "(math.exp("+str(St)+"*"+val+"))*"
	else:
		print "%s is not a recognised parameter\nOnly GC, Sc and St are recognised parameters\n\n" %parameter

def get_p_value(shuffle_type, real_likelihood, accuracy, GC_in, Ns_in, Te_in): #
	counter = 1
	better_shuffled = 0
	while counter <= accuracy:
		shuf_equat = write_equation(shuffle_type)
		def equation_shuffle(x):
			GC, Sc, St = x
			try:
				return eval(shuf_equat)
			except OverflowError:
				return 10000000000000000000
			except ZeroDivisionError:
				return 10000000000000000000
			except ValueError:
				return 10000000000000000000
		x0 = np.asarray((GC_in, Ns_in, Te_in))
		res_shuf = optimize.fmin(equation_shuffle, x0, disp=False)
		likelihood = equation_shuffle(res_shuf)
		if likelihood <= real_likelihood:
			better_shuffled = better_shuffled + 1
		counter = counter + 1
	p_value = '%.2f'%(float(better_shuffled)/accuracy) #change this to %.3f if you change accuracy to 1000
	return p_value

def write_equation(model_type, GC, Sc, St):
	global GC_unshuffled
	global Nitrogen_selection_unshuffled
	GC_unshuffled = {}
	Nitrogen_selection_unshuffled = {}
	for codon in alphabetic_codons:
		GC_unshuffled[codon] = GC_content[codon[0]] + GC_content[codon[1]] + GC_content[codon[2]]
		Nitrogen_selection_unshuffled[codon] = N_content[codon[0]] + N_content[codon[1]] + N_content[codon[2]]	
	if 'shuffle' in model_type:
		global tAI_shuffled
		global GC_shuffled
		global Nitrogen_selection_shuffled
		if len(sys.argv) > 2 and "tRNAscan" in sys.argv[2]:
			tAI_shuffled = using_shuffle(tAI_value)
		GC_shuffled = using_shuffle(GC_unshuffled)
		Nitrogen_selection_shuffled = using_shuffle(Nitrogen_selection_unshuffled)
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
		if parameter == 'shuffle':
			shuffle = 'on'
		else:
			for codon in alphabetic_codons:
				if shuffle == 'off':
					additional = get_upper_equation(codon, parameter, GC, Sc, St, 'for_real', 'off')
				else:
					additional = get_upper_equation(codon, parameter, GC, Sc, St, 'for_shuffle', 'on')
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
			GC = 0.1
			St = 0.1
			Sc = 0.1
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
				try:
					tRNA_count = tRNA_count + 1
					amino_a = symbols[bits[4]]
					p_index = protein_index[amino_a]
					ac_index = anti_codon_index[an]
					#print amino_a, an, p_index, ac_index, anticodons[amino_a], anticodon_count[p_index][ac_index]
					if complementary_nuc[an[2]]+complementary_nuc[an[1]]+complementary_nuc[an[0]] in codon_trans_standard and an in anticodons[amino_a] and bits[8]>45:
						anticodon_count[p_index][ac_index] = anticodon_count[p_index][ac_index] + 1
					else:
						print line
				except KeyError:
					f3.write(line+" includes unknown anticodon so is excluded from the analysis\n")
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
	print "Determining sequence optimality (~5 seconds per gene but varies with length).\n"
	get_tAI_values(tSCAN_file)
	f1=open(species+"_OptimisationResults.txt", "w")
	f1.write("Accession\t%Both_optimised\t%Cost_optimised\t%tAI_optimised\n")
	for acc in cds_filtered:
		sequence = cds_filtered[acc]
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
		d1_distance = []
		while (i < len(cost_frontier)):
			d1_distance.append(math.sqrt((cost - cost_frontier[i])*(cost - cost_frontier[i]) + (translatability  - trans_frontier[i])*(translatability  - trans_frontier[i])))
			i += 1
		d1 = min(d1_distance)
		i = 0
		distance = []
		while (i < len(cost_frontier_worst)):
			distance.append(math.sqrt((cost - cost_frontier_worst[i])*(cost - cost_frontier_worst[i]) + (translatability  - trans_frontier_worst[i])*(translatability  - trans_frontier_worst[i])))
			i += 1
		d4 = min(distance)
		Both_optimised = float(d4*100)/(d1 + d4)
		Cost_optimised = float((1-cost)*100) #This is cost optimisation not relative cost ie. a relative cost of 1 would be the worst optimised.
		tAI_optimised = float(translatability*100)
		Both_optimised = round(Both_optimised, 2)
		Cost_optimised = round(Cost_optimised, 2)
		tAI_optimised = round(tAI_optimised, 2)
		f1.write(acc+"\t"+str(Both_optimised)+"\t"+str(Cost_optimised)+"\t"+str(tAI_optimised)+"\n")
	f1.close()

global GC_to_use
global Sc_to_use
global St_to_use
GC_to_use = 'GC'
Sc_to_use = 'Sc'
St_to_use = 'St'
if "-GC" in sys.argv:
	GC_index = command_line.index("-GC") + 1
	GC_specified = sys.argv[GC_index]
	if GC_specified != "gw":
		GC_specified = float(GC_specified)
		GC_to_use = GC_specified
else:
	GC_specified = "No"
if "-Sc" in sys.argv:
	Sc_index = command_line.index("-Sc") + 1
	Sc_specified = sys.argv[Sc_index]
	if Sc_specified != "gw":
		Sc_specified = float(Sc_specified)
		Sc_to_use = Sc_specified
else:
	Sc_specified = "No"
if "-St" in sys.argv:
	St_index = command_line.index("-St") + 1
	St_specified = sys.argv[St_index]
	if St_specified != "gw":
		St_specified = float(St_specified)
		St_to_use = St_specified
else:
	St_specified = "No"

model_to_use, GC_gw, Sc_gw, St_gw = Genome_wide_analysis(CDS_file)

if Sc_specified == "gw":
	Sc_to_use = Sc_gw
	print "\nSc fixed to the genome-wide Sc (%d) for individual gene analysis." %(Sc_gw)
if St_specified == "gw":
	St_to_use = St_gw
	print "\nSt fixed to the genome-wide St (%d) for individual gene analysis." %(St_gw)

if "-ind" in sys.argv:
	if GC_specified == "No":
		GC_to_use = 'GC'
		print "\nGC content not fixed for individual genes.\nI would not recommend this.\nadd -GC gw to command line to run with genome-wide GC values\n"
	elif GC_specified == "gw":
		GC_to_use = GC_gw
		print "\nGC content fixed to the genome-wide GC (%f) for individual gene analysis." %(GC_gw)
	else:
		GC_to_use = GC_specified
		print "\nGC content fixed to %f for individual gene analysis.\n" % (GC_specified)
	per_gene_analysis(model_to_use)
if "-par" in sys.argv:
	run_pareto_optimisation()
	
def get_linear_relationship_tAI_nitrogen_atoms(File):
	print "Calculating the linear relationship between codon cost and efficiency\n"
	Nit = []
	Te = []
	bad_translate = 0
	f1=open(species+"_tAI_vs_nitrogen_rel.txt", "w")
	f1.write("Amino\tCodon\tRel_N_atoms\tRel_tAI\n")
	f2=open("Not_fully_translatable.txt", "w")
	#f3=open(species+"_tAI_vs_nitrogen_ab.txt", "w")
	#f3.write("Amino\tCodon\tN_atoms\ttAI\n")
	for prot in proteins:
		if prot != 'B': #and prot != 'M' and prot != 'W':
			for codon in proteins[prot]:
				#St_value = tAI_value[codon]
				#Sc_value = N_content[codon[0]] + N_content[codon[1]] + N_content[codon[2]]
				#f3.write(prot+"\t"+codon+"\t"+str(Sc_value)+"\t"+str(St_value)+"\n")
				overall_nit = []
				overall_Te = []
				for syn in proteins[prot]:
					St_value = tAI_value[syn]
					overall_Te.append(St_value)
					Sc_value = N_content[syn[0]] + N_content[syn[1]] + N_content[syn[2]]
					overall_nit.append(Sc_value)
				if tAI_value[codon] == 0:
					St_value = 0
					print codon, prot
					f2.write(species+"\t"+codon+"\t"+prot+"\n")
					bad_translate = bad_translate + 1
				else:
					St_value = float(tAI_value[codon])/max(overall_Te)
					St_value = round(St_value, 3)
				Sc_value = float(N_content[codon[0]] + N_content[codon[1]] + N_content[codon[2]])/max(overall_nit)
				Sc_value = round(Sc_value,3)
				Nit.append(Sc_value)
				Te.append(St_value)
				f1.write(prot+"\t"+codon+"\t"+str(Sc_value)+"\t"+str(St_value)+"\n")
	slope, intercept, r_value, p_value, std_err = stats.linregress(Nit,Te)
	r_squared = r_value*r_value
	slope = round(slope,3)
	intercept = round(intercept,3)
	r_squared = round(r_squared, 3)
	p_value = round(p_value,3)
	std_err = round(std_err, 3)
	f1.write("Species\tSlope\tIntercept\tR_squared\tp_value\tstd_err\n")
	f1.write(species+"\t"+str(slope)+"\t"+str(intercept)+"\t"+str(r_squared)+"\t"+str(p_value)+"\t"+str(std_err)+"\n")
	f1.close
	f2.close
	#f3.close

if "-trade" in sys.argv:
	get_linear_relationship_tAI_nitrogen_atoms(tSCAN_file)
	
print "CodonMuSe (1) implements the SK model (2).\nWhen publishing work that uses CodonMuSe please cite both:\n(1) Seward EA and Kelly S (2017) bioRxiv doi.org/10.1101/136861\n(2) Seward EA and Kelly S (2016) Genome Biology 17(1):226\n"
