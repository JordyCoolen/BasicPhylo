#BASICPHYLO V0.2
#New implementation with mash

## Table of contents
* [Background](#GENERAL-INFO)
* [Output example](#Output example)
* [Authors](#Authors)
* [Database](#Database)
* [Literature](#Literature)

## Background
BASICPHYLO is a tool that is able to identify raw Illumina WGS sequence data from mycobacteria
based on its kmer content using mash.
By combining mash with a database of Mycobacteria phylogeny a UPGMA tree
can be calculated to see batch matching Mycobacteria species.

## Output example
![Alt text](images/Sample.png?raw=true "Output example")

## Authors
J.P.M. Coolen 

Pieter Koopman

Jodie Schildkraut (v0.1)

Paul Verhoeven (v0.1)

## Database
Tool currently supports use of 3 databases:

| Name       | Date | Content      | Notes                                        |
|------------|----|--------------|----------------------------------------------|
| Tortolli_10K.msh | x  | 149 strains  | https://doi.org/10.1016/j.meegid.2017.10.013 |
| 20230303_78_strains_10K.msh   | 20230303 | 78 strains   | expension                                    |
| 20230303_database_10K.msh   | 20230303 | 227 strains  | Complete database                            |

## Literature
Enrico Tortoli, Tarcisio Fedrizzi, Conor J. Meehan, Alberto Trovato, Antonella Grottola, Elisabetta Giacobazzi, Giulia Fregni Serpini, Sara Tagliazucchi, Anna Fabio, Clotilde Bettua, Roberto Bertorelli, Francesca Frascaro, Veronica De Sanctis, Monica Pecorari, Olivier Jousson, Nicola Segata, Daniela M. Cirillo,
The new phylogeny of the genus Mycobacterium: The old and the news,
Infection, Genetics and Evolution, Volume 56, 2017, Pages 19-25, ISSN 1567-1348,
https://doi.org/10.1016/j.meegid.2017.10.013.

mash
Ondov, B.D., Treangen, T.J., Melsted, P. et al.
Mash: fast genome and metagenome distance estimation using MinHash. 
Genome Biol 17, 132 (2016). https://doi.org/10.1186/s13059-016-0997-x