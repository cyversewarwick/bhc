#!/usr/bin/Rscript

library(argparse)
library(RColorBrewer)
library(BHC)

#argument parsing
parser = ArgumentParser()
parser$add_argument('file', nargs=1, help='CSV with gene expression, first column gene names, first row time points')
parser$add_argument('--Pool', dest='pool', type='integer', default=1, help='Number of threads for BHC parallelisation. Default: 1')
parser$add_argument('--Mode', dest='mode', default='multinomial', help='Mode of operation (multinomial/time-course). Default: multinomial')
parser$add_argument('--MakeHeatmap', dest='heatmap', action='store_true', help='Flag. If specified, a heatmap of the resulting clustering will be produced')
args = parser$parse_args()

#basic data prep
data = read.csv(args$file,header=TRUE,row.names=1,check.names=FALSE)
genes = rownames(data)
samples = colnames(data)
data = data.matrix(data)

#BHC runs proper
if (args$mode == 'multinomial')
{
	#data discretisation
	percentiles = FindOptimalBinning(data, genes, transposeData=TRUE, verbose=TRUE)
	discreteData = DiscretiseData(t(data), percentiles=percentiles)
	discreteData = t(discreteData)
	hc = bhc(discreteData, genes, dataType=args$mode, numThreads=args$pool, verbose=TRUE)
	if (args$heatmap)
	{
		#cluster samples too, for the heatmap
		percentiles2 = FindOptimalBinning(data, genes, transposeData=FALSE, verbose=TRUE)
		discreteData2 = DiscretiseData(data, percentiles=percentiles2)
		discreteData2 = t(discreteData2)
		hc2 = bhc(discreteData2, samples, dataType='multinomial', numThreads=args$pool, verbose=TRUE)
		png('heatmap.png')
		heatmap(discreteData, Colv=hc2, Rowv=hc, scale="none", col=brewer.pal(11,'RdBu'))
		dev.off()
	}
} else {
	#standardise data
	standardisedData = (data-mean(data))/sd(data)
	#get them time points as time points
	samples = as.numeric(samples)
	#clustering proper
	hc = bhc(standardisedData, genes, timePoints=samples, dataType=args$mode, numThreads=args$pool, verbose=TRUE)
	if (args$heatmap)
	{
		png('heatmap.png')
		heatmap(standardisedData, Colv=NA, Rowv=hc, scale="none", col=brewer.pal(11,'RdBu'))
		dev.off()
	}
}

#write out clusters
WriteOutClusterLabels(hc, 'clusters.txt', verbose=TRUE)