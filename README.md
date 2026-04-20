# Nanopore metagenomics pipeline (Nextflow)

Here is a Nextflow pipeline to assemble Nanopore reads from isolates containing bacterial species. It is currently being adapted from shell scripts based on (this repository). 

 Main tasks are QC, assembly, annotation and identification. It is configured for both local (desktop) and HPC, showing tractable scaling potential. 

I have modified it to run samples in parallel, employ both conda and apptainer software sources, and use higher throughput tools. I have also simplified its result and artifact products. 

## Pipeline 

This graph is updated as the workflow expands to incorporate original pipeline elements. 

<img src="fig/dag.png" height=300 alt="Nextflow DAG">

## Performance 

![Nextflow timeline analysis screenshot](fig/timeline.png)

## Data 
<!-- This markdown code only seems to work on github -->
* Development: Human oral swabs from (SRA)[^1][^2][^3] downsampled for scaling analysis only, not the data from the original pipeline (unpublished)[^4].
* Production: ...


## Links and References

<!-- This markdown code only seems to work on github -->
[^1]: [bioproject/PRJNA624185](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA624185)
[^2]: [Baker JL *et al.*](https://www.ncbi.nlm.nih.gov/pubmed/33239396), "Deep metagenomics examines the oral microbiome during dental caries, revealing novel taxa and co-occurrences with host molecules.", Genome Res, 2021 Jan;31(1):64-74
[^3]: [Original study github](https://github.com/jonbakerlab/nanopore-oral-genomes). Not implemented here (see[^4]).
[^4]: [Original pipeline](https://github.com/karla-vasco/nanopore_bacterial_isolates) built for different dataset.
