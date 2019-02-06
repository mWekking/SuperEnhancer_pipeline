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

-i

A list of files peaks are called for. on every row should be the relative path from the pipeline to the file. If peaks should be called with a control file, the row should contain the relative path to the sample and the relative path to the control file seperated by comma.

**optional flags**

-o

The name of the output directory. Default = Peakcalling_Output

**Example input**: ```python superEnhancerPipeline.py callpeaks peakCall_filelist.txt -o Macs_peakcalling ```

## filterPeaks ##
Convert the narrowPeak file from the peakcalling into a gff file (needed for ROSE). optional: stretch peaks that are to small, filter out peaks that are to close to a TSS and remove peaks that are on mitochondrial DNA.

**required positional flags**

-i

A list of files that are converted. Every row should contain the relative path from the pipeline to the narrowPeak file

-r

The relative path to the refSEQ txt file that contains all the regulatory elements

**optional flags**

-o/-output

The name of the output directory. Default = filteredPeaks

-p/-minPeakSize

If peaks are smaller than the given size, they will be stretched. default = 2000 bp

-tr/-treshold

If peaks are closer than the given distance to a TSS, they will be filtered out. default = 5000

**-mt/-removeMT**

Remove all peaks that are from the mitochondria (True/False). default = True

**Example input**: ```python superEnhancerPipeline.py filterPeaks filterPeaks_filelist.txt ncbiREFSEQ.txt -o SEpeaks_filtered -p 5000 -tr 10000 -mt False```

## ROSE ##

Find super enhancers by calling the ROSE_main() program

**required positional flags**

-i

A list of files that are used for ROSE. Each line should contain the relative path to each gff file and the bam file it is associated with, comma seperated.

**optional flags**

-g/-genomeBuild

Which genome build is used to map against (MM8, MM9, MM10, HG18, HG19). default = HG19

**Example input**: ```python superEnhancerPipeline.py ROSE ROSE_filelist.txt -g HG18```

## peakConsensus ##

Creating one consensus file containing only regions that are found in a certain amount of the samples

**required positional flags**

-i

A list of the relative paths to every _SuperEnhancers.table.txt file created by ROSE.

-n

The name of the bed file created (_consensus will be added. I.E. NAME_consensus.bed

-t

how many files should contain a certain peak to be considered true

**Example input**: ``` python superEnhancerPipeline.py peakConsensus peakConsensus_filelist.txt Analysis_1 3```

## compareSamples ##

Find per super-enhancer how many reads are from every seperate sample

**required positional flags**

-i

A list of files. Every row should contain the original bam file, the consensus file and in which group the sample belongs comma seperated (I.E. PBT1.bam,PBT_SFT_consensus.bed,A_PBT)

-c

The relative location of the consensus file

**Example input**: ```python superEnhancerPipeline.py compareSamples compareSamples.txt consensusFiles/Analysis_1_consensus.bed```

## DESeq2 ##

Compare two sample groups and find if certain super-enhancers appear to have a higher activity in either one of the groups
There are no flags but in the same folder as the pipeline and the DESeq2 R script should be three files. superEnhancer_counts.txt, superEnhancer_names.txt & results.txt. These files are generated in the comparareSamples step.

**Example input**: ```python superEnhancerPipeline.py DESeq2```

## filterResults ##

filter the DESeq2 results that aren't significant or below a certain fold change. The DESeq2 files should be in the same folder as the pipeline.

**Optional flags**

-s/-signficane

Filter out peaks that have a higher P-value than the given value

-a/adjustedP

Filter out the results peak higher adjusted P-value than the given value. If this flag is used the -s flag will be ignored

-f/-foldChange

filter out peaks that have a fold change higher than the given value

**Example input**: ```python superEnhancerPipeline.py filterResults -s 0.05 -f 0.4```

## findGenes

find genes that are associated with each super-enhancer

**Required positional flags**

-i 

The relative location of the results peak file, generated by DESeq2 and possibly filtered by the filterResults step.

-r

The relative location of the refrence genome file

**Optional flags**

-t/-threshold

The maximal distrance between a TSS and the super-enhancer to be associated with each other. default = 50.000 bp

-n/-geneName

A translation file in TSV format with the columns: RefSeq mRNA ID and Gene name, tab delimited. 

**Example input** ```python superEnhancerPipeline.py findGenes significant_adjFoldChange_results.txt ncbiREFSEQ.txt -t 100000 -n geneTranslateFile.txt```


# Optional steps

## filterChromosomes

Filter certain chromosomes out of your narrowPeak file

**Required positional flags**

-l

The relative path to the directory with narrowPeak files

-c

The chromosomes you want to keep comma seperated

**Example input**: ```python superEnhancerPipeline.py filterChromosomes Macs_peakcalling```

## sortGFF

Sort the .gff file on chromosome and location

**Required positional flags**

-i

A file list with on each row the relative path to the GFF files

-o

The name of output directory

**Example input**: ```python superEnhancerPipeline.py sortGFF sort_filelist.txt sorted_GFFs```

## filterResults

filter our DESeq2 results that aren't significant or below a certain log2foldchange

-s/-significance

Filter out peaks that have a higher P-value than the given value

-a/-adjustedP

 Filter out the result peaks that have a higher adjusted P-value than the given value. If this flag is used, the -s flag will be ignored
 
 -f/-foldChange
 
 Filter out the peaks that have a foldchange lower than the given value
 
 **Example input**: python superEnhancerPipeline.py filterResults -a 0.1 -f 0.4
