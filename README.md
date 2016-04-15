# Skyline-Parent-Mass-List-Filtering
This R-script filters parent mass lists created by Skyline. It can be used to discard charge states, unlikely to detect depend on peptide sequence and/or can discard unmodified peptides when looking for certain modifications


## Introduction

Targeted MS-aquisition by using parent mass lists (PML) can boost your sensitivity tremendously. With The Skyline software it's easy to create these PML for your protein(s) of interest. A possible drawback of this method is the possible length of these lists as they can extend the speed capabilities of your LC-MS setup. This tool reduces the PML dramatically by discarting peptide charge states, which are unlikely to detict according to the peptide sequence. Large repositries of unmodified as well as phosphorylated peptides have been scanned for most sequence dependend charge states and the results serve as input for the filtering process. This filtering can be applied to cover > 95 % of all peptide charge states, which already greatly reduces complexity of the PML, or to just filter for the most likely charge state per peptide.
A second option allows to specifically filter for modified peptides (rightnow phosphorylation and acetylation are supported). In addition charge state filtering will be applied based on empirical data (see above) and peptide bearing modified cleavage sites are discarted as well (e.g. acetylated lysine no longer serves as tryptic cleavage site).
If the PML is used for parrallel reaction monitoring (PRM) in an Q-exactive mass spectrometer optimal MSX combinations can be introduced by retention times predicted by Skyline or more roughly by m/z values.

## Prerequisites

TBD

## Running the script

TBD
