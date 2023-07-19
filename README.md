# GSEA2
Module for GSEA Julia Implementation


# Algorithm Version and Language
  - Original source code: https://github.com/KwatMDPhD/GSEA.jl
  - GSEA.jl 2022.11.14 RC1 Build
  - Julia version 1.8.2
  - Current Docker image used: acastanza/gsea2:2022.11.14

# How to run via commandline
`python3 run.gsea2.py $PWD <expression.dataset> <gene.sets.database> <number.of.permutations> <phenotype.labels> <reverse.phenotypes> <permutation.type> <collapse.dataset> <chip.platform.file> <metric.for.ranking.genes> <enrichment.algorithm> <weighting.exponent> <max.gene.set.size> <min.gene.set.size> <seed.for.permutation> <override.gene.list.length.validation> <plot.graphs>`
  - Refer to parameters below

# Parameters
**\* --> Required Parameter**
- **Basic Parameters**
  - **expression dataset\***: Input expression dataset in gct or tsv format
  - **gene sets database\***: Gene sets database from GSEA website. Upload a gene set if your gene set is not listed as a choice from MSigDB.
  - **number of permutations\***: Number of permutations to perform, default = 1000
  - **phenotype labels\***: Cls file - .cls
  - **reverse phenotypes\***: By default GSEA computes differential expression internally for the first phenotype in the cls vs. the second. Setting "reverse.phenotypes" to "True" computes for the second phenotype vs. the first. This will have no effect on "numeric" format CLS files. Default = False.
  - **permutation type\***: Type of permutations to perform, default = phenotype
  - **collapse dataset\***: Select a mathematical option to collapse a dataset from Gene IDs or Microarray probe IDs to Gene Symbols as used in MSigDB. For Human RNA-seq datasets, "Sum_of_probes" is recommended. For most other cases "Max_probe" should be used. If an option other than "none" is selected, a CHIP file containing the ID to symbol mappings must be provided to the "chip platform file" parameter. Default = none.
  - chip platform file: A Gene ID annotation file from GSEA website. Upload your own chip file if the one corresponding to your platform is not listed in the drop-down menu. A chip file is only required if collapse dataset is set to something other than "none".
- **Algorithmic**
  - **metric for ranking genes\*** : Class separation metric - gene markers are ranked using this metric to produce the gene list.
  - **enrichment algorithm\***: Which algorithm to use for running GSEA the old Kolmogorov-Smirnov enrichment statistic (Classic GSEA), an "area under the curve" variant of the Kolmogorov-Smirnov enrichment statistic, or "CIDAC" (cumulative information divergence with antisymmetricity and complementation) the next-generation GSEA method.
  - **weighting exponent\***: Exponent for weighting the Kolmogorov-Smirnov or Kolmogorov-Smirnov Area enrichment method.
Setting 1.0 is equivalent to the old "scoring scheme=weighted" option.
Setting 0 when running the Kolmogorov-Smirnov enrichment method is equivalent to the unweighted Kolmogorov-Smirnov statistic ("scoring scheme=classic" option in prior versions of GSEA).
  - **max gene set size\***: Gene sets larger than this are excluded from the analysis
  - **min gene set size\***: Gene sets smaller than this are excluded from the analysis
  - **seed for permutation\***: Numerical seed used to initiate the random permutation matrix.
  - **override gene list length validation\***: GSEA will check the length of the input dataset to ensure a reasonable number of genes are included.
- **Reporting**
  - **plot graphs\***: Number of top gene sets to produce enrichment plots for.


## Contact Information
  - Anthony Castanza (acastanza@cloud.ucsd.edu)
  - Edwin Huang (edh021@cloud.ucsd.edu)
