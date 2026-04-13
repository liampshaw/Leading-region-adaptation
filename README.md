# The leading region of many plasmids is adapted for translational efficiency

A repository of analysis related to the paper:

**The leading region of many plasmids is adapted for translational efficiency**  
Liam P. Shaw, Arancha Peñil Celis, Manuel Ares Arroyo, Lijuan Luo, Tatiana Dimitriu, Fernando de la Cruz  
doi: [short_text](url) 

Associated datasets can be downloaded at Zenodo: [short_text](url)

**A note on versions.** This repository supersedes an earlier repository associated with a previous iteration of this work, whose code is archived [here](https://github.com/liampshaw/RM-and-leading-region) (Zenodo doi: [10.5281/zenodo.15235227](https://doi.org/10.5281/zenodo.15235227)). That repository corresponds to v1 of the preprint and is preserved as a record of the previous analysis for clarity.

## Details

To be added.

Scripts:
`Fig1_` makes Figure 1 etc.

`seq_stats/` contains Rust scripts in `src/bin/` for sequence statistics on batched plasmid sequences: GC content (`gc_content.rs`) and palindrome density (`palindrome_density.rs`). These were created using generative AI coding tools then validated for accuracy, including against older results generated from Python scripts. 


## Data 

`TableSX_defensefinder_results_prodigal_meta.tsv` - results of running DefenseFinder vX.X on plasmid proteins identified with Prodigal vX.X. Complete systems only.
`seq_stats`: 
* `avg_gc.tsv` - results from running Rust script on plasmids for GC content (give command)


`targets_leading_lagging`: all the comparisons of leading/lagging regions for within-range RM targets for the 13 PTUs analysed

To be added/clarified.

`1751_plasmids_oriented.tar.gz` - n=1751 plasmids oriented so that leading strand comes first
`1751_defense_systems.tsv` - results of `defense-finder -a` (with antidefense systems detected) on n=1751 plasmids
*to be added* - CDS of all n=1751 plasmids



We then count palindromes using `scripts/rust/src/palindrome_density.rs`. This script can work for any value of k - we focus on k=4,6,8,10 because Type II RM targets are typically 4-8bp (and most are 5-6bp). We can use a sliding window or not, but we want to avoid autocorrelation for some statistics so in that case we don't use a sliding window. 

