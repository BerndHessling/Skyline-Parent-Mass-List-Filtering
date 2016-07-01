# Skyline-Parent-Mass-List-Filtering
Reduce complexity of parent mass lists by e.g. filtering out unlikely peptide charge states.

## Introduction

Targeted MS-aquisition by using parent mass lists (PML) can boost your sensitivity tremendously. With The Skyline software it's easy to create these PML for your protein(s) of interest. A possible drawback of the method is length of these lists as they can extend the speed capabilities of your LC-MS setup. This tool reduces the PML size dramatically by discarting peptide charge states, unlikely to be detected according to peptide sequences. Large repositries of unmodified as well as phosphorylated peptides have been scanned for most sequence dependend charge states and the results serve as input for the filtering process. This filtering can be applied to cover > 95 % of all peptide charge states, which already greatly reduces complexity of the PML, or to just filter for the one most likely detectable charge state per peptide.
A second option allows to specifically filter for modified peptides (right now STY-phosphorylation and K-acetylation are supported). Additionally charge state filtering will be applied based on empiric data (see above) and peptides bearing modified cleavage sites are discarted as well (e.g. acetylated lysine no longer serves as tryptic cleavage site).
If the PML is used for parrallel reaction monitoring (PRM) in an Q-exactive mass spectrometer optimal MSX combinations can be introduced by retention times predicted by Skyline or more roughly by m/z values.

## Prerequisites

 - Skyline is recommended to generate raw PML
 - [R](https://cran.r-project.org/bin/windows/base/) needs to be installed (package was developed for version 3.2.4).
 - csv file exported from skyline or formatted in the same way (see [pdf](Creating-parent-mass-lists-with-skyline.pdf) for more details)

## Running the script

When starting the script three successive input fields will pop-up:
 1. "Input file": Browse to your Skyline exported csv file (see [pdf](Creating-parent-mass-lists-with-skyline.pdf) for more details)

 2. "Which peptides to keep:": Choose the filter applied to your peptide list:
   - "all"         : don't discard peptides
   - "STY-phospho" : keep only peptides with phosphorylation on S,T, or Y
   - "K-acetyl"    : keep only peptides with acetylation on K, but no acetylation on c-terminal K

 3. "Choose Charge state filter:" How restrictively should be filtered for charge states:
   - "optimal" : only charge states are allowed that are likely to be detected according to the number of basic amino acids and phosphorylations (empiricly derived by analysing large dataset of about 15k peptides)
   - "oneMax"  : only the one most likely charge State is accepted
   - "noFilter": all charge states are accepted

 4. "Choose MSX count:" Should optimal MSX IDs be included (only for PRM analysis):
   - "0": disabled
   - "1-10": number of MSX groups in your experiments

MSX IDs are selected by having optimal spread of peptides for each MSX group according to either retention time (if present in the input file) or m/z values.

## Output

The resulting csv file will be filterd with chosen settings and stored at the same place the original file was located. It can directly be imported into QExactive instruments from thermo, or m/z lists must be copied manually into the MS-software.

