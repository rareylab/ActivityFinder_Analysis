# ActivityFinderAnalysis

This repository is designed to analyse the ActivityFinder results and compare them to the sota, BioChemGraph.
We utilize both the built ActivityDB as well as output files of both the BioChemGraph and the ActivityFinder.
The structure of the project is as following:

## Environment

### pyproject.toml

Here all the project information including the required packages are listed

## Folder structure

### input_files

Input files for both ActivityFinder and BioChemgraph. 
The combined output files for ActivityFinder will be available with the StrAcTable publication, which should be published immediately after.

### state_presubmission
This is where all output files go. 
You could change the name in the definitions.py file as I did multiple times while this project went on.
Within that folder there are multiple other folders

#### output_files

Most output files.

### automated_analysis_outputs

The outputs of fully automated analysis like those in the create_sankey_diagram.py


### molecule_difference_data

All the output files and plots regarding differences in molecule data between biochemgraph and ActivityFinder

### plots

All the plots that are relevant

# Code

The code has historically grown and I now need to investigate remains of past maybe intelligent life (me).
All code is in the src folder and there some code is split up dependent on the specific parts of the problem it looks at.

# Querying the rest api of ActivityFinder at PP2

The file ActivityFinderRestApi.ipynb details how to query the rest api of the ActivityFinder and get the results for single PDB entries.

# General stuff

### Definitions.py

All the definitions and paths I need, also the swith from home to zbh for Database access

### create_sankey_diagram.py

This is more or less the complete analysis and creates a sankey diagram as a more complex venn diagram.
It details all reasons or descriptors of the data BioChemGraph finds, we find and both find.

### compare_jsons.py

Compares the jsons created of the ActivityFinder and BioChemGraph.
Creates all the nice venn diagrams

### extract_json_from_activityfinder.py

Extracts all jsons from the ActivityFinder results

### extract_json_from_biochemgraph.py

Extracts all jsons from BioChemgraph input files

### cheatsheet.py

Creates basic stats and venn diagrams for the comparison of BioChemGraph and ActivityFinder

### utils

Is a collection of functionality used by multiple scripts. Mostly functionality to split data depending on a property.

### utils/db.py

Collection of utility querying the different relevant databases.

### utils/getter.py

Collections of functions getting specific data.
Most functions are designed in a way that they only calculate or query the data if it is not present for the user.

### utils_/isualisation.py

Visualisation functions like venn diagrams and molecule drawing

### analyse_notfound_targets.py

Analyses the notfound targets that were not found in targetmode.
Tries to find the reasons why it does so, including sequence identity.

## Mutation analysis

There are some exploratory python scripts for mutation analysis in this directory.
They were mainly used to find cases with interesting mutations to showcase in the Paper.

# Sequence Stuff

### sequence_identity_plots.py

Creates additional plots about sequence identity data

### sequence_matching.py

Constructs the match between two sequences and gets the percent identity.

### find_slightly_different_target_from_biochem.py

This script is used to find perfect or almost perfect targets that the BioChem does not find.

# Molecule Stuff

### find_moldiff_examples_db.py

This is the complete analysis on the molmatching between chembio and ActivityFinder.

### generate_moldifference_example_plots.py

Makes single examples of two inchis, also outdated but still useful I think.
Precedes the BioChemGraph analysis.

### paper_molmatching_plots.py

Current paper molmatch plot script with all the interesting single cases.

# Execution

The first files you have to execute are the ones extracting the jsons from biochemgraph and ActivityFinder.
The latter takes a long time, the former is pretty fast.
Afterwards, execute compare_jsons.py to compare the jsons generated in the prior runs.

# Data Availability and ZBH Particularities in the Code

## ZBH Particularities
All scripts were written with the given zbh infrastructure in mind, so we use our db server whakaari, our ActivityDb instance, etc. 
If you want to get this code running on your machine, you will have to change these particularities

## Data Availability
Some of the used data to analyze why the ActivityFinder did not work will be made available shortly after the publiction
through the follow-up Publication about StrAcTable. I will update this repo with a link as soon as the data is public.