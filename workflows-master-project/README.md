# workflows: Example workflows that utilize datarail pipelines and/or external resources.

## Installation instructions

While individual workflows can be run independently, it is strongly recommended to install the R package associated with this repository. The package will automatically install all the necessary dependencies and provide a `getData()` function that can be used to quickly retrieve the example datasets used by individual workflows. To install the package, simply run the following commands in your R shell:

    if( !require(devtools) ) install.packages( "devtools" )
    devtools::install_github( "datarail/workflows" )

## Digital Gene Expression

The following workflows are centered around Digital Gene Expression data analysis:

* Aligning FASTQ files to the reference genome is detailed in a separate datarail repo: https://github.com/datarail/rnaseq
* After aligning raw reads, the next step is typically converting the counts table to a proper genes-by-samples format.
* One of the first things to do a dataset is to perform Principal Components Analysis (PCA) and plot the resulting projection(s).
* [In development] Differential expression
* [In development] Gene Set Enrichment Analysis

