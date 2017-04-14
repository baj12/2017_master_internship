# Cassandra Lab Meeting
_Date : 03/10/2017_

## Task 01 : 
- Showing the phylogenetic distribution of Oskar in insect orders
- We also want to see if the evolution of oskar is in the same than the evolution of insect orders


### Sub-task 01: Update our data
__From Tamsin work: We have 181 genomes__

With our script (based on Genomes data from NCBI), we have found __247 genomes__ 

$\rightarrow$ Result : Stock on a CSV file wirh species name/id, order name/id, family name/id, number acession

_NB: If the Genome data in NCBI database, the script still works for updating_


### Sub-task 02: Creating EVM / Augustus pipeline

#### EVM
Build a genome consensus sequence with
- a genome sequence (Fasta file)
- gene prediction (GFF3 file)
- protein / transcript alignments (GFF3 file)
- a weight values to be applied to each type of evidence

#### Augustus
- Create an hidden markov model for a given sequence (here our gene consensus sequence)
- Predict all gene for a specific species beased on the HMM model 


### Sub-task 03 : Transcriptomes
Leo have done the work
- Have to be curated
- Maybe do a better analysis

## TO DO
$\Rightarrow$ Test the pipeline on one genome

## Questions
About HMM model for each order
- Choose one organism randomly?
- Choose 4/5 organisms per order to build a meta HMM? (maybe use a representative organism for each family per order)
