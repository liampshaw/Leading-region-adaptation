# The leading region of many plasmids is adapted for translational efficiency

A repository of analysis related to the paper:

**The leading region of many plasmids is adapted for translational efficiency**  
Liam P. Shaw, Arancha Peñil Celis, Manuel Ares Arroyo, Lijuan Luo, Tatiana Dimitriu, Fernando de la Cruz  
doi: [short_text](url) 

Associated datasets can be downloaded at Zenodo: [short_text](url)

**A note on versions.** This repository supersedes an earlier repository associated with a previous iteration of this work, whose code is archived [here](https://github.com/liampshaw/RM-and-leading-region) (Zenodo doi: [10.5281/zenodo.15235227](https://doi.org/10.5281/zenodo.15235227)). That repository corresponds to v1 of the preprint and is preserved as a record of the previous analysis for clarity.

## Detail

**`scripts/`**
* `Fig1_` makes Figure 1 etc.
* `seq_stats/` contains Rust scripts in `src/bin/` for sequence statistics on batched plasmid sequences: GC content (`gc_content.rs`) and palindrome density (`palindrome_density.rs`). These were created using generative AI coding tools then validated for accuracy, including against older results generated from Python scripts in the previous repository.
* `cluster_runDefenseFinder.sh` shows how we ran DefenseFinder 

**`data/`**
* `TableS1_` corresponds to Table S1 of the paper etc. 
* A few other miscellaneous small files are included
* Larger datasets that are needed to remake the figures can be downloaded from zenodo [here]() and should then be placed in a `data/zenodo` directory before running the figure scripts. 

