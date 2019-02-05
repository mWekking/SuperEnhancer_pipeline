# Introduction
Going from raw ChIP-seq reads to super-enhancer regions and comparing these regions within two groups can be a very tedious job. This pipeline was created to do this in a relative easy to use way. The pipeline features 8 clear steps (and 2 optinonal steps) that go from raw reads to super-enhancer regions and their associated genes. 
# Installation
This pipeline makes use of different libraries within python and Linux. Needed to run the full pipeline is:
python related:
```
- Python (>= 2.6 and < 3)
- argparse module
- os module
- sys module
- subprocess module
- time module
```

Linux related:
```
- Macs2
- samtools
- bedtools
``` 

and the most recent version of R

# Usage

```

superEnhancer_Pipeline_split.py [-h] {callpeaks,filterChromosomes,filterPeaks,sortGFF,ROSE,peakConsensus,compareSamples,DESeq2,filterResults,findGenes}

```
The pipeline has 8 command to find super-enhancers.

| Sub command | Discription |
| --- | --- |
| callPeaks | Calling peaks by making use of the macs peakcaller |
| filterChromosomes | **Optional**: filter certain chromosomes out I.E chromosome X and Y |
| filterPeaks | stretch peaks, remove peaks to close to TSS, remove MT dna and convert to gff |
| sortGFF | **Optional**: Sort the created gff file on chromosome and starting coordinates |
| ROSE | Finding the super-enhancer regions by making use of ROSE |
| peakConnsensus | Creating one consensus file containing only regions that are found in a certain amount of the samples |
| compareSamples | Find for each super-enhancer how many reads are found in each sample |
| DESeq2 | Use DESeq2 to find the log2foldchange per group | 
| filterResults | **Optional**: Filter out results that aren't significant or have a to low difference in either group |
| findGenes | Finding which genes are associated with the significant regions | 

# Required Steps

## callPeaks ##
The first step makes use of the mac peak calling algorithm. This is possible to do with or without input/control files. All other settings are the macs2 default setting.

**required positional flags**

**-i**

A list of files peaks are called for. on every row should be the relative path from the pipeline to the file. If peaks should be called with a control file, the row should contain the relative path to the sample and the relative path to the control file seperated by comma.

**optional flags**

**-o**

The name of the output directory. Default = Peakcalling_Output

**Example input**: ```python superEnhancerPipeline.py callpeaks peakCall_filelist.txt ```

## filterPeaks ##
Convert the narrowPeak file from the peakcalling into a gff file (needed for ROSE). optional: stretch peaks that are to small, filter out peaks that are to close to a TSS and remove peaks that are on mitochondrial DNA.

**required positional flags**

**-i**

A list of files that are converted. Every row should contain the relative path from the pipeline to the narrowPeak file

**-r**

The relative path to the refSEQ txt file that contains all the regulatory elements

**optional flags**

**-o/-output**

The name of the output directory. Default = filteredPeaks

**-p/-minPeakSize**

If peaks are smaller than the given size, they will be stretched. default = 2000 bp

**-tr/-treshold**

If peaks are closer than the given distance to a TSS, they will be filtered out. default = 5000


**-mt/-removeMT**

Remove all peaks that are from the mitochondria (True/False). default = True

**Example input**: ```python superEnhancerPipeline.py filterPeaks filterPeaks_filelist.txt ncbiREFSEQ.txt -o SEpeaks_filtered -p 5000 -tr 10000 -mt False



# Optional steps

