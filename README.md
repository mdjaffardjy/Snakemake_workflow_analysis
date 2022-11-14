
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![Version 1.0.1](https://img.shields.io/badge/version-v1.0.1-blue)]()


# Snakemake_workflow_analysis

In this repository is presented a Snakemake workflow analyser. It takes one or several snakefiles and it has three major libraries :

-  The snakemake parser that  allows users to extract information (info on processes or and the structure) of a single or multiple Snakemake workflows. It is found in the folder [src](/src/parsing_snkmk.py). It extracts the snakemake components of a workflow (notably processors) and helps extracting info on the workflow such as the number of processors (here, rules), as well as their code, name, inputs and outputs. It doesn't extract the structure. It also finds the tools present in each processor.
- The second one, [__search_biotools__](/src/search_biotools_dump.py) finds tools matches and informations in bio.tools using a dump found in the inputs.
- The last one, [__get_characteristics__](/src/get_characteristics.py) compiles in a dataframethe extracted information about the processors (rules) previously parsed with parsing_snkmk. 

 

## Contribute
Please submit GitHub issues to provide feedback or ask for new features, and contact us for any related question.


## Install and Run

### Conda environment

The conda environment can be found [here](/data/env_conda.txt).

### Source data

The source data can be found in the [crawl file](/data/wf_crawl_snakemake.json). This file contains git metadata as well as links for the workflows. One can download the Snakefiles from this json. It is heavily recommended to store the data in the form /<name owner>/<name project>/<number>.snakefile

### To run
The analyzer needs as input the address of a folder containing all the snakefiles to analyze.

In a script, import [parsing_snkmk](/src/parsing_snkmk.py) and use these two functions:

```
import parsing_snkmk as ps

path_wf = "path/to/snakefiles/folder"

# workflow parsing
wf_list = ps.parse_wf_folder(path_wf)
dict_tool, total_extr = ps.make_dict_wf_with_biotools_dump(wf)

```
Your results will be in the object wf_list. You will have additional information about the tools you extracted in dict_tool and total_extr.



