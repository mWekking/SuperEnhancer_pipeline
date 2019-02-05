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

| Sub command | Discription |
| --- | --- |
| callPeaks | Calling peaks by making use of the macs peakcaller |
| filterChromosomes | **Optional**: filter certain chromosomes out I.E chromosome X and Y |
| filterPeaks | stretch peaks, remove peaks to close to TSS, remove MT dna and convert to gff |
| sortGFF | **optional**: Sort the created gff file on chromosome and starting coordinates |
