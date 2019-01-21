import argparse
import os
import sys
import subprocess
import time

#PYTHON 2

#The location of the pipeline.py program, this is used troughout the progam to locate directories and files
programDirectory = str(os.path.abspath(__file__))[:-len(os.path.basename(__file__))]

#----------------------------------------------------------- call peaks ---------------------------
def locateFileTypes(fileType, location):
    locatedFiles = []
    filesInDir = os.listdir(programDirectory + location)
    for index, files in enumerate(filesInDir):
        if files[-len(fileType):] == fileType:
            locatedFiles.append(filesInDir[index])
    return locatedFiles

def checkDir(dirName):
    if os.path.isdir(dirName):
        pass
    else:
        subprocess.call(["mkdir", dirName])

def checkFile(fileName):
    if os.path.isfile(fileName):
        return True
    else:
        return False

def peakcalling(fileDict, outputDir):
    for sample in fileDict:
        if fileDict[sample] != None:
            print("---  Calling peaks for {}    ---\n".format(sample))
            subprocess.call(["macs2", "callpeak", "-t", sample \
                                                , "-c", fileDict[sample] \
                                                , "-n", sample[:-4] \
                                                , "--outdir", outputDir \
                                                , "-f", sample[-3:].upper() \
                                                , "--verbose", "2"])
            print("---  Finished calling peaks for {} ---".format(files))
        else:
            print("---  Calling peaks for {}    ---\n".format(files))
            subprocess.call(["macs2", "callpeak", "-t", sample \
                                                , "-n", sample[:-4] \
                                                , "--outdir", outputDir \
                                                , "-f", sample[-3:].upper() \
                                                , "--verbose", "2"])
            print("---  Finished calling peaks for {} ---".format(sample))
    return outputDir


#--------------------------------------------------- filter chromosome -----------------------------
def filterChromosome(chromosomeList, peakCalldir):
    #This step is used to Filter out certain chromosomes from the .narrowPeak file that is created with the peakcalling step.
    #This is nothing more than a simple AWK command with first creating the filter string and after this performing the awk step on every .narrowPeak file in the given directory
    #You could do this by hand but this saves alot of tedious typing
    try: 
        int(chromosomeList[0])
        filterString = "$1 == {}".format(chromosomeList[0]) 
    except ValueError:
        filterString = "$1 == \"{}\"".format(chromosomeList[0].upper()) 

    if len(chromosomeList) > 1:
        for chromosome in chromosomeList[1:]:
            try: 
                int(chromosome)
                filterString += "||"
                filterString += "$1 == {}".format(chromosome)
            except ValueError:
                filterString += "||"
                filterString += "$1 == \"{}\"".format(chromosome.upper())
        
        
    for files in locateFileTypes(".narrowPeak", peakCalldir):
        if "awked" not in files:
            with open(peakCalldir + "/" + files[:-11] + "_awked.narrowPeak", "wb") as f:
                subprocess.call(["awk", filterString, peakCalldir + "/" + files], stdout = f)

#-------------------------------------- filter peaks ---------------------------------

def createTranslateDict(trnsFile):
    #Create a dictionary in the form of {NR_xxxx.x : ENSTxxxxx.x, NM_xxxx.x : ENSTxxxxxx.x}
    #This is used to translate the gene names in the refrence genome file
    translateDict = {}
    with open(trnsFile, "r") as f:
        for line in f:
            splitLine     = line.strip().split()
            value = splitLine.pop(0)
            for key in splitLine:
                translateDict[key] = value
    return translateDict

def createGeneDict(refSeqFile):
    #create a Dictionary with the form of: {"X" : {ENSTxxxx.x : {bin : x, exonEnd : x, exonFramse : x}, ENSTxxxx.x : {bin : x, exonEnd : x, exonFramse : x} }
    # This is used to calculate the distance of each gene to the peaks that are found
    refSeq = {}
    refSeq     = {"X" : {}, "Y" : {}, "MT" : {}}
    for key in range(1,23):
        refSeq[str(key)] = {}
    refHeaders = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"]
    
    with open(refSeqFile, "r") as f:
        for line in f:
            tempDict = {}
            splitLine= line.rstrip().split()

            key = splitLine[1]
            
            for header, value in zip(refHeaders, splitLine):
                tempDict[header] = value

            #Removing al genes on chromosome Un (chrUn) and on locations like chr6_apd_hap1
            if "_" in tempDict["chrom"] or tempDict["chrom"][0:5] == "chrUn":    
                break                  
            else:
                refSeq[str(tempDict["chrom"][3:])][key] = tempDict
    return refSeq

def createPeakDict(peakFile):
    #Creates a dictionary from the peak file with  the form of: {File_Peakx : {"chromEnd" : x, "peak" : x}, File_Peakx : {"chromEnd" : x, "peak" : x} }
    peakDict    = {}
    peakHeaders = ["chrom", "chromStart", "chromEnd", "score", "strand", "signalValue", "pValue", "qValue", "peak"]
    with open(peakFile) as f:
        for line in f:
            tempDict = {}
            splitLine= line.rstrip().split()
            key = splitLine.pop(3)
            for header, value in zip(peakHeaders, splitLine):
                tempDict[header] = value

            peakDict[key] = tempDict
    return peakDict

def stretchPeaks(peakDict, minPeakSize):
    #stretch Peaks that are smaller than a certain treshold (minPeakSize) to the treshold size
    for pKeys in peakDict:
        peakStart = int(peakDict[pKeys]["chromStart"])
        peakEnd   = int(peakDict[pKeys]["chromEnd"])
        peakSize = peakEnd - peakStart
        if peakSize < minPeakSize:
            dif = (minPeakSize - peakSize) / 2
            if peakStart - dif > 0:
                peakDict[pKeys]["chromEnd"] = int(peakDict[pKeys]["chromEnd"]) + dif
                peakDict[pKeys]["chromStart"] = int(peakDict[pKeys]["chromStart"]) - dif
            else:
                peakDict[pKeys]["chromEnd"] = int(peakDict[pKeys]["chromEnd"]) + 2 * dif
                peakDict[pKeys]["chromStart"] = int(peakDict[pKeys]["chromStart"])
    return peakDict

def median(lst):
    #used to calculate the median without importing numpy as there are many reasons why not to use numpy.. 
    #Thats not true I just had some problems with installing numpy and writing this median program was easier than solving the errors
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0

def filterPeakDistance(peakDict, refSeqDict, treshold, removeMTdna):
    #filters out any peak that is too close (treshold) to any promoter sequence

    #endPeaks will be a dictionary that only contains the peaks that are further than the treshold distance from every gene
    endPeaks = {}
    counter = 0
    for pKey in peakDict:
        # When running this prints some progress information
        counter += 1
        if float(counter) / 200 == counter / 200:
            os.system("clear")
            perc = (float(counter) / float(len(peakDict))) * 100
            print ("%.1f" % perc) + "%\n" , str(counter) + "/" + str(len(peakDict)), "Peaks tested", "\nFar enough peaks", len(endPeaks)
            print("\n--- %s seconds ---" % (time.time() - start_time))

        #here the comparrison starts
        #the Enddistance is the closest distance from each peak to a gene in the refrence genome
        #The EndDistance always starts at a ridiculous high number so that the first distance calculated is always lower.  
        #endDistance = 10 ** 20
        endDistance = float("inf")
        for rKey in refSeqDict[peakDict[pKey]["chrom"]]:
            if endDistance > treshold and treshold != 0:
                #If the Peak is withing treshold distance of a gene in the refrence genome the whole peak gets skipped. It's to close to a gene so no need to test futher it gets thrown out
                    distance  = None
                    
                    peakStart = int(peakDict[pKey]["chromStart"])
                    peakEnd   = int(peakDict[pKey]["chromEnd"])
                    
                    if refSeqDict[peakDict[pKey]["chrom"]][rKey]["strand"] == "+":
                        refStart  = int(refSeqDict[peakDict[pKey]["chrom"]][rKey]["txStart"])
                        refEnd    = int(refSeqDict[peakDict[pKey]["chrom"]][rKey]["txEnd"]) 
                    elif refSeqDict[peakDict[pKey]["chrom"]][rKey]["strand"] == "-":
                        refStart  = refEnd    = int(refSeqDict[peakDict[pKey]["chrom"]][rKey]["txEnd"]) 
                        refEnd    = refEnd    = int(refSeqDict[peakDict[pKey]["chrom"]][rKey]["txStart"]) 

                    if refStart < peakStart < refEnd:
                        distance = peakStart - refStart
                    elif refStart < peakEnd < refEnd:
                        endDistance = 500
                        break
                    elif peakStart < refStart and peakEnd > refEnd:
                        peakMedian = median(list(range(peakStart, peakEnd + 1)))
                        if refStart > peakMedian:
                            distance = refStart - peakMedian
                        else:
                            distance = peakMedian - refStart
                    elif refStart > peakEnd:
                        distance = refStart - peakEnd    
                    elif peakStart > refEnd:
                        distance = peakStart - refStart
                        
                    if distance < endDistance:
                        endDistance = distance
            else:
                break

        if removeMTdna == True:
            if endDistance > treshold:
                if peakDict[pKey]["chrom"] == "MT":
                    pass 
                else:    
                    endPeaks[pKey] = peakDict[pKey]
        elif removeMTdna == False:
            if endDistance > treshold:
                endPeaks[pKey] = peakDict[pKey]
        
    return endPeaks

def writeGFF(endPeaksDict, outputDir, fileName):
    #write filtered peaks into .gff file
    with open("{}/{}/{}_filtered.gff".format(os.path.dirname(os.path.realpath(__file__)), outputDir, fileName), "wb") as f:
        for endData in endPeaksDict:
            seqName   = endPeaksDict[endData]["chrom"]
            source    = endData
            feature   = "."
            start     = endPeaksDict[endData]["chromStart"]
            end       = endPeaksDict[endData]["chromEnd"]
            score     = "."
            strand    = endPeaksDict[endData]["strand"]
            frame     = "."
            attribute = "."
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seqName, source, feature, start, end, score, strand, frame, attribute))

#----------------------------------------- sort .gff File ---------------------
def sortGFF(gffFilesList, outputFolder):
    #This simply sorts the .gff files, again this can be done by hand but this saves alot of typing
    for files in gffFilesList:
        with open(str(os.path.abspath(os.curdir)) + "/" + files[:-4] + "_sorted.gff", "wb") as f:
            subprocess.call(["sort", "-k1,1", "-k4,4", "-n", files], stdout = f)

#-------------------------------------------- call ROSE ----------------------------

def ROSEcall(rankByDict, genome , outputDir = "ROSE_files"):
    #The most timeconsuming step of the pipeline, here the actual Super enhancers are found.
    #This is done by using the ROSE program, programmed by the Young lab 
    #http://younglab.wi.mit.edu/super_enhancer_code.html
    #I added this step because it's is very time consuming.. so you can simply start this step and let the program run for the night

    for bamFile in rankByDict:
        subprocess.call(["python", "ROSE_main.py", "-g", genome, "-i", rankByDict[bamFile], "-r", bamFile, "-o", outputDir ])

#--------------------------------------------find consensus Peaks -------------------------------
def getFileList(fileListLocation):
    fileList = []
    with open(fileListLocation, "r") as f:
        for fileLocation in f:
            fileList.append(fileLocation.rstrip().replace("\\", r"/"))
    return fileList

def bedParser(inputFile): 
    outputDict = {}
    with open(inputFile, "r") as f:
        for lines in f:
            if lines[0] == "#":
                pass
            else:
                if lines[:9] == "REGION_ID":
                    bedHeader = lines.split()
                else:
                    tempDict = {}
                    splitLine = lines.split()
                    for header, value in zip(bedHeader, splitLine):
                        tempDict[header] = value
                    try:
                        outputDict[tempDict["CHROM"]][tempDict["REGION_ID"]] = tempDict
                    except KeyError:
                        outputDict[tempDict["CHROM"]] = {}
                        outputDict[tempDict["CHROM"]][tempDict["REGION_ID"]] = tempDict
    return outputDict

def bedOverlap(bedFileList):
    genomeList = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
    chromosomes = None
    for bedfile in bedFileList:
        if len(bedfile.keys()) > chromosomes:
            chromosomes = bedfile.keys()
    peakList = [[] for _ in range(len(chromosomes))]
    for bedfile in bedFileList:
        for x in genomeList:
            if x not in bedfile:
                bedfile[x] = {}
        for chromosome in chromosomes:
            for peak in bedfile[chromosome]:
                start = int(bedfile[chromosome][peak]["START"])
                stop = int(bedfile[chromosome][peak]["STOP"])
                peakListIndex = chromosomes.index(bedfile[chromosome][peak]["CHROM"])
                match_idxs = []
                for regionidx, region in enumerate(peakList[peakListIndex]):
                    if (    region[0] < start <= region[-1] or
                            region[0] < stop <= region[-1] or
                            start < region[0] < stop or
                            start < region[-1] < stop):
                        match_idxs.append(regionidx)

                if len(match_idxs) == 1:
                    region = peakList[peakListIndex][match_idxs[0]]
                    region.extend([start, stop])
                    peakList[peakListIndex][match_idxs[0]] = sorted(region)
                elif len(match_idxs) > 1:
                    target_idx = match_idxs.pop(0)
                    region = peakList[peakListIndex][target_idx]
                    for source_idx in reversed(match_idxs):
                        region.extend(peakList[peakListIndex].pop(source_idx))
                    region.extend([start, stop])
                    peakList[peakListIndex][target_idx] = sorted(region)
                else:
                    peakList[peakListIndex].append([start, stop])
    return chromosomes, peakList

def fileWriter(peakList, chromosomeList, filename, treshold):
    with open("{}/{}_consensus.bed".format(os.path.abspath(os.curdir), filename), "wb") as f:
        for chromosomeNumber, peak_per_chr in zip(chromosomeList, peakList):
                for peak_region in peak_per_chr:
                    if len(peak_region) >= 2 * treshold:
                        start = min(peak_region)
                        stop  = max(peak_region)
                        toWrite= [str(chromosomeNumber), str(start), str(stop)]
                        f.write("\t".join(toWrite) + '\n')

#--------------------------------------------- compareGroups --------
def compareBamBed(bamBedFilesList):
    createdBedFiles = []
    for bamFile in bamBedFilesList:
        fileName       = bamFile[:-4] + "_Compared.bam"
        bedFile = bamBedFilesList[bamFile][0]
        relativePath = os.path.dirname(os.path.realpath(__file__)) + "/"
        with open(relativePath + fileName, "wb") as f1:
            print "Working on comparing {} with {} ...".format(bamFile, bedFile)
            subprocess.call(["samtools", "view", relativePath + bamFile, "-L", relativePath + bedFile, "-b"], stdout= f1)
            print "created {}".format(fileName)
        with open(relativePath + fileName[:-4] + ".bed", "wb") as f2:
            createdBamfile = [relativePath + fileName, bamBedFilesList[bamFile][1]]
            subprocess.call(["bedtools", "bamtobed", "-i", createdBamfile[0]], stdout = f2)
            print "converted it to {}".format(fileName[:-4] + ".bed")
            createdBedFiles.append([relativePath + fileName[:-4] + ".bed", createdBamfile[1]])
        print "Done\n"
    return createdBedFiles
        
def superEnhancerReadOverlap(consensusFile, readsFiles, pValueTreshold, foldChangeTreshold):
    relativePath = os.path.dirname(os.path.realpath(__file__)) + "/"
    consensusDict = {}
    with open(relativePath + consensusFile, "r") as f:
        for line in f:
            chromosome, start, stop = line.rstrip().split()
            if chromosome not in consensusDict:
                consensusDict[chromosome] = []
            consensusDict[chromosome].append((int(start), int(stop)))
    
    header = ["Chr", "Start", "Stop"]
    countHeader = ["Enhancer"]
    writeNameLines  = []
    writeNameLines.append("Name\tSample")

    samplesDicts = []
    sampleName    = {}
    for files in readsFiles:
        tempDict = {}
        with open(files[0], "r") as f:
            if files[1] not in sampleName:
                sampleName[files[1]]    = 1
                header.append(files[1] + str(1))
                countHeader.append(files[1] + str(1))
                writeNameLines.append("{}\t{}".format(files[1] + str(1), files[1]))
            else:
                sampleName[files[1]] += 1
                header.append(files[1] + str(sampleName[files[1]]))
                countHeader.append(files[1] + str(sampleName[files[1]]))
                writeNameLines.append("{}\t{}".format(files[1] + str(sampleName[files[1]]), files[1]))       
            for line in f:
                chromosome, start, stop, cigar, score, strand = line.rstrip().split()
                if chromosome not in tempDict:
                    tempDict[chromosome] = []
                tempDict[chromosome].append((strand, int(start), int(stop)))
        samplesDicts.append(tempDict)

    writeLines      = []
    writeCountLines = []
    counter = 1
    for chromosomes in consensusDict:
        for superEnhancers in consensusDict[chromosomes]:
            superEnhancerStart = superEnhancers[0]
            superEnhancerStop = superEnhancers[1]
            scores = []
            for samples in samplesDicts:
                score = 0
                for reads in samples[chromosomes]:
                    readStrand             = reads[0]

                    if readStrand == "+" or readStrand == ".":
                        readStart          = reads[1]
                        readStop           = reads[2]
                    else:
                        readStart          = reads[2]
                        readStop           = reads[1]

                    if (superEnhancerStart <= readStart < superEnhancerStop or
                        superEnhancerStart < readStop  <= superEnhancerStop or
                        readStart <= superEnhancerStart < readStop          or
                        readStart <  superEnhancerStop <= readStop):
                        score += 1

                scores.append(score)

            writeLines.append("{}\t{}\t{}\t{}".format(chromosomes, superEnhancerStart, superEnhancerStop, "\t".join(map(str, scores))))
            writeCountLines.append("Chr{}_Peak{}\t{}".format(chromosomes, counter, "\t".join(map(str, scores))))
            counter += 1
    
    with open(relativePath + "reads_per_SuperEnhancers.bed", "wb") as f:
        f.write("\t".join(header) + '\n')
        for lines in writeLines:
            f.write(lines + "\n")

    with open(relativePath + "superEnhancer_names.txt", "wb") as f:
        for lines in writeNameLines:
            f.write(lines + "\n")               
    
    with open(relativePath + "superEnhancer_counts.txt", "wb") as f:
        f.write("\t".join(countHeader) + "\n")
        for lines in writeCountLines:
            f.write(lines + "\n")
                
    subprocess.call(["Rscript", "--save", "DESeq2_analysis.R"])
    
    with open(relativePath + "results.txt", "r") as f1, open(relativePath + "reads_per_SuperEnhancers.bed", "r") as f2:
        with open("results_withCoords.txt", "wb") as wf:
            #remove the Headers
            f1.readline()
            f2.readline()
            wf.write( "\t".join(header) + "\n")
            for lines1, lines2 in zip(f1, f2):
                writeLine       = lines1.rstrip().split()
                coordinatesLine = lines2.rstrip().split()
                start           = coordinatesLine[1]
                stop            = coordinatesLine[2]
                writeLine.extend([start, stop])
                wf.write( "\t".join(writeLine) + "\n")

    if pValueTreshold[0] == "a":
        treshold     = float(pValueTreshold[1])
        filterString = "$7<{}".format(treshold)
        with open("Significant_results.txt", "wb") as f:
            subprocess.call(["awk", filterString, relativePath + "results_withCoords.txt"], stdout= f)
    elif pValueTreshold[0] == "s":
        treshold     = float(pValueTreshold[1])
        filterString = "$6<{}".format(treshold)
        with open("Significant_results.txt", "wb") as f:
            subprocess.call(["awk", filterString, relativePath + "results_withCoords.txt"], stdout= f)

    if foldChangeTreshold[0] == True:
        if pValueTreshold[0] == True:
            fileName = "Significant_results.txt"
        else:
            fileName = "results_withCoords.txt"
            
        treshold     = float(foldChangeTreshold[1])
        filterString = "$3>{}||$3<{}".format(treshold,-treshold)
        with open(fileName[:-4] + "_adjFoldChange.txt", "wb") as f:
            subprocess.call(["awk", filterString, relativePath + fileName], stdout= f) 

#----------------------------------------------------------find genes

def readConsensus(consensusFile):
    endDict = {}
    with open(consensusFile, "r") as f:
        for lines in f:
            tempDict = {}
            name, baseMean, logChange, lfcSE, stat, pvalue, padj, start, stop = lines.rstrip().split()
            tempDict = {"name" : name, "baseMean" : baseMean, "logChange" : logChange, \
                        "lfcSE": lfcSE, "stat" : stat, "pvalue" : pvalue, "padj" : padj, \
                        "start": start, "stop" : stop}
            chromosome = name[3:5].strip("_")
            if chromosome not in endDict:
                endDict[chromosome] = []
            endDict[chromosome].append(tempDict)
    return endDict

def geneFinder(consensusDict, geneDict, geneDistance = 50000):
    #find genes that are on the same chromosome and have their TSS withing the geneDistance of a super-enhancer
    with open("Significant_results_adjFoldChange_withGenes.txt", "wb") as f, open("Gene_centered_file.txt", "wb") as f2:
        f.write("ID\tchrom\tstart\tstop\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tgenes\n")
        f2.write("gene\tdistance\tchrom\tstart\tstop\tID\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n")
        for chromosomes in consensusDict:
            for peaks in consensusDict[chromosomes]:
                geneList   = []
                endGeneDict   = {}

                peakStart   = int(peaks["start"])
                peakEnd     = int(peaks["stop"])
                peakMiddle  = median(list(range(peakStart,peakEnd)))

                for genes in geneDict[chromosomes]:                                        
                    if geneDict[chromosomes][genes]["strand"] == "+":
                        TSS = int(geneDict[chromosomes][genes]["txStart"])
                    if geneDict[chromosomes][genes]["strand"] == "-":
                        TSS = int(geneDict[chromosomes][genes]["txEnd"])
                    
                    if (peakStart <= TSS < peakEnd or \
                        (peakStart - TSS <= geneDistance and peakStart - TSS > 0) or \
                        (TSS - peakEnd <= geneDistance and TSS - peakEnd > 0)):  
                        if TSS < peakMiddle:
                            distance = peakMiddle - TSS
                        elif peakMiddle < TSS:
                            distance = TSS - peakMiddle
                        endGeneDict[genes] = distance
                        
                for gene, distance in sorted(endGeneDict.iteritems(), key=lambda (k,v): (v,k)):
                    geneList.append(gene) 
                    f2.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene, distance, chromosomes, peaks["start"], peaks["stop"], peaks["name"], peaks["baseMean"], peaks["logChange"], peaks["lfcSE"], peaks["stat"], peaks["pvalue"], peaks["padj"]))

                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(peaks["name"], chromosomes, peaks["start"], peaks["stop"], peaks["baseMean"], peaks["logChange"], peaks["lfcSE"], peaks["stat"], peaks["pvalue"], peaks["padj"], ",".join(geneList)))

#--------------------------------------------------------- input ------------------------------
#Here the input from the ubuntu console gets parsed
parser = argparse.ArgumentParser(description= "This is a easy to use pipeline, to call peaks, find significant super-enhancers and the genes they are associated with.")

subparsers = parser.add_subparsers()
parser_peakCalling = subparsers.add_parser("callpeaks", help= "Use Macs2 to call peaks on bamfiles")
#required flags
parser_peakCalling.add_argument("i" , help= ".txt filelist file. each row must be a sample and a control file, seperated by a comma")
parser_peakCalling.add_argument("-o", help= "The name of the output folder, This will be located in the folder where this .py file is located" )

parser_filterChromosomes = subparsers.add_parser("filterChromosomes", help= "Filter certain chromosomes out of your .narrowPeak file")
#required flags
parser_filterChromosomes.add_argument("l", help= "The directory name where the .narrowPeaks files are located")
parser_filterChromosomes.add_argument("c", help= "which chromosomes do you want to extract? i.e 1,2,4")

parser_filterPeaks = subparsers.add_parser("filterPeaks", help= "Convert .narrowPeak file to .gff file with the possibilitie to filter out peaks that are too close to a TSS or stretch peaks to a certain size")
#Required flags
parser_filterPeaks.add_argument("i", nargs=1, type=str, help = "The .txt file containing the relative paths from this program and containing the names of each .narrowPeaks file")
parser_filterPeaks.add_argument("r", nargs=1, type=str, help = "This is the refSeq .txt file that contains all regulatory elements")
#Optional flags
parser_filterPeaks.add_argument("-o", "-output",      nargs=1, type=str, help   = "The location of the output folder")
parser_filterPeaks.add_argument("-p", "-minPeakSize", nargs=1, type=int, help   = "If you want to stretch the peaks smaller than a certain size to other size. default = 2000" )
parser_filterPeaks.add_argument("-tr", "-treshold",    nargs=1, type=int, help  = "The minimal distance each enhancer should be of any promoter region. default = 5000")
parser_filterPeaks.add_argument("-mt", "-removeMT", nargs=1, type=str, help    = "Remove Mitochondrial DNA from the peaks. (True/False) default= True")
#input example: python argpars.py SFTregD2_peaks.narrowPeak ncbiREFSEQ.txt -o results/filteredPeaks -p 2500 -tr 5100'

parser_sortGFF = subparsers.add_parser("sortGFF", help = "Sort the .gff file on chromosome and on position")
#required flags
parser_sortGFF.add_argument("i", nargs=1, type=str, help = "The .txt file containing the relative paths from this program and containing the names of each .gff file")
parser_sortGFF.add_argument("o", nargs=1, type=str, help = "The name of the outputFolder")

parser_ROSE = subparsers.add_parser("ROSE", help = "Find super enhancers by calling the ROSE_main() program")
#required flags
parser_ROSE.add_argument("i", nargs=1, type=str, help = "The .txt file containing the relative paths from this program and containing the names of each sorted .gff file you want to find super enhancers in. on the next line after eachr .gff file should be the name of the .bam file that is assiociated with it")
parser_ROSE.add_argument("-g", "-genomeBuild", type=str, help = "Which genome build is used to map to. (choose: MM8,MM9,MM10,HG18,HG19) default= HG19")

parser_peakConsensus = subparsers.add_parser("peakConsensus", help = "create one consensus .bed file out of multiple superEnhancer.table.txt")
#required flags
parser_peakConsensus.add_argument("i", nargs=1, type=str, help = "The location of the files list. This list should contain the relative paths leading to the .bed Super enhancer files from the location of this program")
parser_peakConsensus.add_argument("n", nargs=1, type=str, help = "This wil be the name of the end .bed file (to this name '_consensus' will be added)")
parser_peakConsensus.add_argument("t", nargs=1, type=int, help = "In how many files should a certain peak overlap to be added to the consensus .bed file")

parser_compareGroups = subparsers.add_parser("compareGroups", help = "Compare two sample groups and find if certain super-enhancers appear to have a higher activity in either one of the groups")
#required flags
parser_compareGroups.add_argument("i", nargs=1, type=str, help = "A file containing all the bam files, bed files and sample groups. all bam files, bed files and sample groups should be comma seperated and each couple should be on a newline. The sample groups are compared alphabeticly, so a positive foldchange means that there is more activity in group B. two example lines could be: PBT1.bam,PBT_SFT_consensus.bed,A_PBT & SFT1.bam,PBT_SFT_consensus.bed,B_SFT")
parser_compareGroups.add_argument("c", nargs=1, type=str, help = "The relative location of the consensus super enhancer file")
#optional flags
parser_compareGroups.add_argument("-s", "-significane", nargs=1, type=str, help= "Filter out peaks that have a higher P-value than the given value")
parser_compareGroups.add_argument("-a", "-adjustedP",   nargs=1, type=str, help= "Filter out the result peaks that have a higher adjusted P-value than the given value. If this flag is used, the -s flag will be ignored")
parser_compareGroups.add_argument("-f", "-foldChange",  nargs=1, type=str, help= "Filter out the peaks that have a foldchange lower than the given value")

parser_findGenes = subparsers.add_parser("findGenes", help = "Read the end consensus file and compare it to a refrence genome to find which genes are regulated by the super enhancers that are found")
#required flags
parser_findGenes.add_argument("i", nargs=1, type=str, help = "The relative location of the consensus peak file")
parser_findGenes.add_argument("r", nargs=1, type=str, help = "The location of the refrence genome file")
#optional flags
parser_findGenes.add_argument("-t", "-threshold", nargs=1, type=int, help="The maximal distance between the Gene and the super enhancer. default = 50.000")
args = parser.parse_args()

if sys.argv[1] == "callpeaks": 
    bamFileDir = {}
    fileList = args.i
    with open(fileList) as f:
        for lines in f:
            files = lines.rstrip().split(",")
            if len(files) > 1:
                #This step checks if there is a input/control file added, if not the peaks will be called withoud a input file.
                bamFileDir[files[0]] = files[1]
            else:
                bamFileDir[files[0]] = None

    if args.o:
        peakcallDir = args.o
    else:
        peakcallDir = "Peakcalling_Output"    
    checkDir(peakcallDir)
    
    
    if len(bamFileDir) > 0:
        peakcalling(fileDict= bamFileDir, outputDir= peakcallDir)   
    else:
        print("Error: There where no files found to peform peakcalling on ..")

if sys.argv[1] == "filterChromosomes":
    peakcallDir = args.l
    chromeList  = args.c.split(",")
    filterChromosome(chromosomeList= chromeList, peakCalldir= peakcallDir)
                                                              
if sys.argv[1] == "filterPeaks":
    #getting input
    narrowPeakFilesList = []
    with open(args.i[0], "r") as f:
        for lines in f:
            narrowPeakFilesList.append(lines.rstrip().replace("\\", r"/"))

    referenceGenome = args.r[0]

    if args.o:
        outPut = args.o[0]
    else:
        outPut = "filteredPeaks"
    checkDir(outPut)

    if args.p:
        inputMinPeakSize  = int(args.p[0])
    else:
        inputMinPeakSize = 2000

    if args.tr:
        inputTreshold = int(args.tr[0])
    else:
        inputTreshold = 5000

    if args.mt:
        removeMT = args.mt[0]
        if removeMT.lower() == "false":
            removeMT = False
        else:
            removeMT = True
    else:
        removeMT = True

    start_time = time.time()
    for peakFile in narrowPeakFilesList:
        fileName = peakFile[-peakFile[::-1].index("/"):][:-11]
        refSeq        = createGeneDict(referenceGenome)
        peakDict      = createPeakDict(peakFile)
        if inputMinPeakSize > 0:
            peakDict      = stretchPeaks(peakDict, minPeakSize = inputMinPeakSize)
        endPeaks      = filterPeakDistance(peakDict, refSeq, treshold = inputTreshold, removeMTdna= removeMT)
        writeGFF(endPeaks, outPut, fileName)
        
    
    print("\n--- %s seconds ---" % (time.time() - start_time))

if sys.argv[1] == "sortGFF":
    gffFiles = []
    with open(args.i[0], "r") as f:
        for lines in f:
            gffFiles.append(lines.rstrip().replace("\\", r"/"))

    sortGFF(gffFiles, args.o[0])

if sys.argv[1] == "ROSE":
    allFiles = []
    with open(args.i[0], "r") as f:
        for lines in f:
            allFiles.append(lines.rstrip().replace("\\", r"/"))
    rankByDict = {}
    for gffFiles, bamFiles in zip(allFiles[::2], allFiles[1::2]):
        rankByDict[bamFiles] = gffFiles
    
    if args.g:
        genomebuild = args.g[0].upper()
    else:
        genomebuild = "HG19"
    
    ROSEcall(rankByDict, genome= genomebuild)

if sys.argv[1] == "peakConsensus":
    fileListDir     = args.i[0]
    outputName      = args.n[0]
    treshold        = args.t[0]
    fileList = getFileList(fileListLocation = fileListDir)
    inputPeaks = []
    for files in fileList:
        bedFile = bedParser(files)
        inputPeaks.append(bedFile)
    chroms, endPeaks = bedOverlap(inputPeaks)

    fileWriter(endPeaks, chroms, outputName, treshold)

if sys.argv[1] == "compareGroups":
    bambedFiles = {}
    with open(args.i[0]) as f:   
        for lines in f:
            bamFile, bedFile, group= lines.rstrip().split(",")
            bambedFiles[bamFile] = [bedFile, group]

    createdBedFiles = compareBamBed(bambedFiles)
    consensusFile   = args.c[0]
    
    if args.a:
        adjustedPValue = ["a", args.a[0]]
    elif args.s:
        adjustedPValue = ["s", args.s[0]]
    else:
        adjustedPValue = [False]

    if args.f:
        foldChangeValue = [True, args.f[0]]
    else:
        foldChangeValue = [False]

    start_time = time.time()
    superEnhancerReadOverlap(consensusFile, createdBedFiles, adjustedPValue, foldChangeValue)
    print("\n--- %s seconds ---" % (time.time() - start_time))

if sys.argv[1] == "findGenes":
    consensusFile     = args.i[0]
    consensusDict     = readConsensus(consensusFile)

    referenceGenome   = args.r[0]
    referenceGenomeDict = createGeneDict(referenceGenome)
    
    geneFinder(consensusDict, referenceGenomeDict)