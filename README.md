# Bayesian Hierarchical Clustering

## The Purpose of the Algorithm

Bayesian Hierarchical Clustering (BHC) is a simple, fast clustering algorithm capable of partitioning datasets made up of multiple static measurements or time courses whilst requiring minimal parameter input.

## Algorithm Availability outside iPlant

The app is a wrapper for the R package, which can be downloaded from Bioconductor and ran locally if more in-depth parameterisation is desired.

## Basic Input/Output

BHC accepts a CSV expression file on input, with the first column being gene identifiers and the first row being numerical time points (if the dataset is a time course) or sample names (if the dataset is a collection of multiple static measurements).

The output comes in the form of a basic cluster list (`clusters.txt`), and is automatically parsed into BiNGO and MEME friendly forms for further investigation, which are provided in the `functional_analysis_inputs/` folder.

## How Does It Work?

The concept behind BHC is an enhancement of the basic bottom-up hierarchical clustering framework. Every gene starts out in its own separate cluster, and the algorithm iteratively merges the best pair of unmerged clusters until a full hierarchical tree is obtained. The hierarchy is then split into a final cluster configuration based on the posterior probabilities of pairs of clusters being merged being the correct thing to do.

This app also features an enhancement to the basic BHC algorithm which allows it to process time course data. The time courses can be treated with two different modes of operation, a squared exponential covariance function or a cubic spline (to make it akin to [SplineCluster][heard2005]).

For details, consult [Cooke et al., 2011][cooke2011].

## Test Run

If you want to get a feel for how BHC operates and what the output looks like without using your own data, you can find a demonstration input file at `ktpolanski/bhc_testdata/botrytis17.csv` under Community Data. Switch the Run Mode field to squared exponential or cubic spline as desired, as this is a time course dataset.

## Input in Detail

### Expression CSV

**Mandatory input.** The data you're going to analyse. The first column is to be gene identifiers, while the first row is to contain information on the samples (if the dataset is comprised of multiple static measurements) or time points (if the dataset is a time course). In case formatting reference is needed, consult `ktpolanski/bhc_testdata/botrytis17.csv` under Community Data for a time course demonstration file.

### Run Mode

Select an operation mode based on the data you're providing. **Multinomial (Static Data)** is appropriate if your dataset is a collection of multiple static measurements, whilst time course datasets have the choice between **Squared Exponential Covariance (Time Course)** or, if emulating spline-based methods is desired, **Cubic Spline (Time Course)**.

### Generate a Heatmap

If checked, a heatmap visualisation of the obtained clustering will be created. If the data is a collection of multiple static measurements, those will be clustered as well, just for heatmap purposes. The shown data is the one used to perform clustering, with the expression data being binned into three categories - red for low expression, white for neutral expression, blue for high expression. The heatmap will likely be very crowded and illegible when large amounts of data are provided on input.

## Output in Detail

### `clusters.txt`

The basic cluster export. Each cluster is denoted by a header line, followed by one line per gene in the cluster.

### `functional_analysis_inputs/`

A folder featuring the clusters reformatted into input formats accepted by BiNGO and MEME.

### `heatmap.png`

If the app was instructed to produce a heatmap, then an additional heatmap image will be generated. The expression data provided on input is discretised by the script, and this is reflected on the heatmap - red means the expression value got put in the low expression bin, white is the middle expression bin, while blue is the high expression bin. If the data is a collection of static measurements, the samples are also clustered for heatmap purposes, but this clustering is not exported in any form.

### `FullOutput.tar`

The complete output of the analysis, archived into a single file for ease of downloading to your computer.

[heard2005]: http://www.pnas.org/content/102/47/16939.short
[cooke2011]: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-399