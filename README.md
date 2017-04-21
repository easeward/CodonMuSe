# CodonMuSe
**_CodonMuSe_** is a method for analysing the factors that affect codon bias both at the genome-wide and individual gene-level. 

**_CodonMuSe_** estimates
1) The strength of selection acting on nucleotide biosynthesis cost.
2) The strength of selection acting on translational efficiency.
3) The proportion of codon bias attributable to mutational processes.

**_CodonMuSe_** determines the extent to which individual gene sequences are
1) Optimised for cost efficiency.
2) Optimised for translational efficiency.
3) Jointly optimised for both cost and translational efficiency.

## How _CodonMuSe_ works
CodonMuSe implements the SK model for inferring the impact on **Codon** bias of **Mu**tation and **Se**lection acting on nucleotide cost and translational efficiency. The SK model was first described in:

**Seward EA, Kelly S (2016)** Dietary nitrogen alters codon bias and genome composition in parasitic microorganisms. **_Genome Biology_** 17(1):226.

The CodonMuSe implementation of the SK model was first published in:

**Seward EA, Kelly S (in prep.)** The evolutionary economics of a gene, bacteria balance the cost and efficiency of mRNA.

## Installing _CodonMuSe_
CodonMuSe is written in python and requires **scipy**. Up-to-date instructions on how to install scipy are provided here: http://www.scipy.org/install.html. Once scipy is installed, the CodonMuSe source code can be downloaded from this repository by clicking the _Clone or download_ link at the top of the page and executed directly as described below.


## Running _CodonMuSe_

`python CodonMuSe.py [OPTIONS] -f <SEQUENCE FILE> -tscan <tRNAscan FILE>`

	-f <FILE>	A FASTA file of protein coding nucleotide sequences
	-tscan <FILE>	A tRNA copy number file produced by tRNAscan
	-tc <INT>	The NCBI genetic code identifier goo.gl/ByQOau (Default = 1)
	-ind		Analyse individual genes in addition to a genomewide analysis
	-par 		Determine the cost and efficiency optimality of individual genes
	-trade          Calculates trade-off between cost and efficiency
	-Mb <float>  	Specify a fixed value for Mb in all calculations
	-Sc <float>  	Specify a fixed value for Sc in all calculations
	-St <float>  	Specify a fixed value for Sc in all calculations
	-Mb gw       	Fix Mb (mutation bias) to genome-wide value for individual genes
	-Sc gw       	Fix Sc (seleciton on cost) to genome-wide value for individual genes
	-St gw       	Fix St (selection on translational efficiency) to genome-wide value for individual genes
	-m <TXT>	Specify model parameters (Mb, Sc, St) eg. Mb_Sc_St or Mb_Sc
	                (Default = determine best parameter combination automatically)

## Input Files

Two input files are required to run CodonMuSe

**1) A FASTA file containing nucleotide sequences for protein coding genes**

Typically this is the full set of protein coding genes from a single species. However, the method can be applied to sequences obtain from transcriptome or metagenome sequencing. By default CodonMuSe only analyses nucleotide sequences that begin with a start codon, end with stop codon, and do not contain in-frame stop codons. Sequences not meeting all of these criteria are ignored and reported in an error file.

**2) tRNA copy number file generated by tRNAscan**

This is a text file containing the output generated by running tRNAscan on the complete genome of the species in question. If you do not include this optional file then CodonMuSe cannot infer the strength of selection acting on codon translational efficiency or determine the extent to which gene sequences have been optimised for translational efficiency.

**Schattner P, Brooks AN, Lowe TM (2005)** The tRNAscan-SE, snoscan and snoGPS web servers for the detection of tRNAs and snoRNAs. **_Nucleic Acids Res_** 33:686–689


## Output files

**\<SEQUENCE FILE\>\_GenomeWideResults.txt** 

This file contains the results of analysing the genome-wide codon bias patterns. The file contains the log likelihood (LnL), AIC, and R2 values of the model as well as the fitted values for each parameter (Mutation bias [Mb], Selection acting on nucleotide cost [Sc] and Selection acting on translational efficiency [St]). Omitted parameters are given a value of zero. The file also contains a tab separated table of observed vs fitted codon frequencies. These frequencies correspond to synonymous codons ie sum to 1 for a given amino acid.

**\<SEQUENCE FILE\>\_IndividualGenesResults.txt**

This file contains the results of analysing each individual gene's codon bias patterns. The file contains a tab separated table which records the gene accession, model log likelihood value and fitted gene-specific values for each parameter (Mutation bias [Mb], Selection acting on nucleotide cost [Sc] and Selection acting on translational efficiency [St]). Omitted parameters are given a value of zero.

**\<SEQUENCE FILE\>__OptimisationResults.txt**

This file contains the results of analysing each individual gene's percentage optimisation. The file contains a tab separated table which records the gene accession, percentage optimisation for both transcript cost and transcript translational efficiency, percentage optimisation for transcript cost (individually) and percentage optimisation for translational efficiency (individually). 

**\<SEQUENCE FILE\>__tRNAscan_errors.txt**

This file contains information about any tRNAs that were not used in the CodonMuse analysis. For example tRNAscan identifies pseudogenes (_pseudo_) and tRNAs for Seleno cystein (_SeC_) that are not used by CodonMuse.

**\<SEQUENCE FILE\>__excluded_sequences.txt**

This file contains information about any sequences that were not used in the CodonMuse analysis. For example very shot genes (\<30nt), genes lacking start codons, genes lacking stop codons, and genes with frame errors are not analysed by CodonMuse.

## Example Dataset

An example dataset is provided. To run CodonMuSe on this dataset execute the following command

`CodonMuSe.py -f M_pneumoniae_CDS.fasta -tscan M_pneumoniae_tRNAscan.txt -tc 4 -ind -Mb gw -par -m Mb_Sc_St`

