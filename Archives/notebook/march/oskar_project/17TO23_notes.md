# Friday, March 17th

## Done 
- [x] Understand all the pipeline and know all the commands to run

## To do
- Choose a model of each order (for determinating if we have to use EVM for all the orders)
- Finally testing EVM !!!!

## New idea
- Do a first filter to check which species has the oskar gene or by not by genome mapping with the oskar gene list
- After do only EMV_Augustus pipeline on species which don't have oskar sequence
    * If the species have already annotations, apply only Augustus
    * Else, apply all the pipeline

## _For the week-end_
- Begin Exonerate execution code


# Monday, March 20th

## To do
- Apply a first filter on genomes to skip which are already known to have oskar sequence

## Done
- [x] Test EVM (but something wrong after execute commands)
- [x] Update of genome\_search (some id are not updated according to eucrayotes to eukaryotes.txt)

## In progress
- Exonerate pipeline for finding Oskar (have to modify Exonerate command to be faster)


# Tuesday, March 21st

## Done
- [x] Have launched Exonerate on 200 genomes 

## In progress
- Create a file which contain infos from all annotations genomes 


# Wednesday, March 22nd

## In progress
- Reboot Exonerate pipeline (still running on the first 8 processus)
- Work on `infoTracker.py`

## To do
- Finish `infoTracker.py`
- Test EVM with data from NCBI


# Thursday, March 23th

## Done
- [x] pipelineGuide created

## Notes
- See two genomes with errors (have to check by hand - fixed if we check only genomic.fna)
    - GCA\_000002195.1 : No genomic.fna.gz (but we have genomic.fna) and h annotations 
    - $\rightarrow$ Use Augustus
    - GCA\_001014675.1 : No genomic.fna.gz (but we have genomic.fna) and no gff annnotations 
    - $\rightarrow$ GeneMark + Exonerate + EVM + Augustus
    
## TO DO                                                               
> - Convert into GeneBank file and use Augustus! $\rightarrow$ __1__  
> - Find Oskar $\rightarrow$ __37__ 

- Already known, probably which one are in Genome\_search 
- Figure out that finally, the Accesion number of the initial Genome_search where good, have to upload the 10 right genomes

> - Use Augustus $\rightarrow$ __207__ 
> - Use GeneMark + Exonerate + Evidence Modeller + Augutus $\rightarrow$ __1__

- Have to create a pipeline only for Augustus
    - Choose model for training (17 for each order at least)
    
## Notes
- We can shut down Exonerate Pipeline (useless)
