#!/usr/bin/Rscript

library(argparse)
library(RColorBrewer)
library(BHC)
library(parallel)
library(gplots)

#argument parsing
parser = ArgumentParser()
parser$add_argument('file', nargs=1, help='CSV with gene expression, first column gene names, first row time points')
parser$add_argument('--Pool', dest='pool', type='integer', default=0, help='Number of threads for BHC parallelisation. Default: 0 (automatic parallelisation)')
parser$add_argument('--Mode', dest='mode', default='multinomial', help='Mode of operation (multinomial/time-course). Default: multinomial')
parser$add_argument('--MakeHeatmap', dest='heatmap', action='store_true', help='Flag. If specified, a heatmap of the resulting clustering will be produced')
args = parser$parse_args()

#toy around with pool
if (args$pool==0)
{
	args$pool = detectCores()
}
#basic data prep
data = read.csv(args$file,header=TRUE,row.names=1,check.names=FALSE)
genes = rownames(data)
samples = colnames(data)
data = data.matrix(data)
#standardise data
standardisedData = (data-mean(data))/sd(data)

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
		png('heatmap.png',width=8,height=6,units='in',res=300,family='mono')
		holder = par()
		rowmar = max(nchar(genes))*(0.2+1/(log10(length(genes))))*holder$cra[1]/holder$cra[2] * 3/4
		colmar = max(nchar(samples))*(0.2+1/(log10(length(samples))))*holder$cra[1]/holder$cra[2] * 2/3
		heatmap.2(standardisedData, Colv=hc2, Rowv=hc, tracecol=NA, scale="none", col=brewer.pal(11,'RdBu'), margins=1.5+c(colmar,rowmar))
		dev.off()
	}
} else {
	#get them time points as time points
	samples2 = as.numeric(samples)
	#clustering proper
	hc = bhc(standardisedData, genes, timePoints=samples2, dataType=args$mode, numThreads=args$pool, verbose=TRUE)
	if (args$heatmap)
	{
		png('heatmap.png',width=8,height=6,units='in',res=300,family='mono')
		holder = par()
		rowmar = max(nchar(genes))*(0.2+1/(log10(length(genes))))*holder$cra[1]/holder$cra[2] * 3/4
		colmar = max(nchar(samples))*(0.2+1/(log10(length(samples))))*holder$cra[1]/holder$cra[2] * 2/3
		heatmap.2(standardisedData, Colv=NA, Rowv=hc, tracecol=NA, scale="none", col=brewer.pal(11,'RdBu'), margins=1.5+c(colmar,rowmar))
		dev.off()
	}
}

#write out clusters
WriteOutClusterLabels(hc, 'clusters.txt', verbose=TRUE)