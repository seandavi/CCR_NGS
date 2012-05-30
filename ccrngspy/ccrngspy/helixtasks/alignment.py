import helix.jobs
import os


class BwaAlignment(helix.jobs.Job):
    threads = 24
    nodes = '1:c24:gpfs'
    
    def __init__(self,config,fastqFiles,bamFile,options=None):
        """
        config: a config instance
        fastqFiles: a list of length 1 or 2 (single- or paired-end)
        bamFile: output bam file name
        options: not used, but could be for changing default BWA behavior"""
        inputs = fastqFiles
        outputs = [bamFile,bamFile + ".finished"]
        hash = {'bwaExe':config['bwaexe'],
                'bwaPrefix':config['bwaprefix'],
                'threads':self.threads,
                'fastqFile':fastqFiles[0],
                'bamFile':bamFile,
                'bamPrefix':bamFile.replace('.bam',''),
                'samtoolsExe':config['samtoolsExe'],
                'saiFile':fastqFiles[0]+'.sai'}
        if(len(fastqFiles)==1):
            cmdString = """#!/bin/bash
            %(bwaExe)s aln -q 2 -t %(threads)d %(bwaPrefix)s %(fastqFile)s > %(saiFile)s 
            %(bwaExe)s samse -n 1 %(bwaPrefix)s %(saiFile)s %(fastqFile)s | %(samtoolsExe)s view -bS - | samtools sort -m 5000000000 - %(bamPrefix)s
            touch %(bamFile)s.finished
            """ % hash
            super(BwaAlignment,self).__init__(command=cmdString,
                                              inputs=fastqFiles,
                                              outputs=[bamFile,bamFile+'.finished'],
                                              nodes=self.nodes)
                                              
        if(len(fastqFiles)==2):
            hash['saiFile2']=fastqFiles[1]+".sai"
            hash['fastqFile2']=fastqFiles[1]
            hash['threads']=int(hash['threads']/2)
            cmdString = """#!/bin/bash
            %(bwaExe)s aln -q 2 -t %(threads)d %(bwaPrefix)s %(fastqFile)s > %(saiFile)s &
            bg
            %(bwaExe)s aln -q 2 -t %(threads)d %(bwaPrefix)s %(fastqFile2)s > %(saiFile2)s &
            bg
            wait
            %(bwaExe)s sampe -n 1 -N 1 %(bwaPrefix)s %(saiFile)s %(saiFile2)s %(fastqFile)s %(fastqFile2)s | %(samtoolsExe)s view -bS - | samtools sort -m 5000000000 - %(bamPrefix)s
            touch %(bamFile)s.finished
            """ % hash
            super(BwaAlignment,self).__init__(command=cmdString,
                                              inputs=fastqFiles,
                                              nodes=self.nodes,
                                              outputs=[bamFile,bamFile+'.finished'])


class MergeBam(helix.jobs.Job):
    threads = 2
    nodes="1:c16:gpfs"

    def __init__(self,config,bamFiles,outBam):
        """
        config: a config instance
        bamFiles: a list of input bam files
        outBam: an output bam file name"""
        inputs = bamFiles
        outputs = [outBam,outBam+'.bai']
        hash = {'picardPrefix':config['picardPrefix'],
                'outputBam':outBam,
                'samtoolsExe':config['samtoolsExe'],
                'inputs':" ".join(["INPUT=" + x for x in bamFiles])}
        cmdString = """#!/bin/bash
        java -jar %(picardPrefix)s/MergeSamFiles.jar ASSUME_SORTED=True OUTPUT=%(outputBam)s %(inputs)s
        %(samtoolsExe)s index %(outputBam)s """ % hash
        super(MergeBam,self).__init__(command=cmdString,
                                      inputs=bamFiles,
                                      nodes=self.nodes,
                                      outputs=outputs)

class MarkDuplicates(helix.jobs.Job):
    threads = 2
    nodes="1:c16:gpfs"

    def __init__(self,config,bamFile,outBam):
        """
        config: a config instance
        bamFile: a SINGLE input bam files
        outBam: an output bam file name"""
        inputs = bamFile
        outputs = [outBam,outBam+'.bai',outBam+".dupmetrics"]
        hash = {'picardPrefix':config['picardPrefix'],
                'outputBam':outBam,
                'metrics':outBam + ".dupmetrics",
                'samtoolsExe':config['samtoolsExe'],
                'inputs':bamFile}
        cmdString = """#!/bin/bash
        java64 -Xmx4g -jar %(picardPrefix)s/MarkDuplicates.jar ASSUME_SORTED=True METRICS_FILE=%(metrics)s OUTPUT=%(outputBam)s INPUT=%(inputs)s CREATE_INDEX=true""" % hash
        super(MarkDuplicates,self).__init__(command=cmdString,
                                            inputs=bamFile,
                                            nodes=self.nodes,
                                            outputs=outputs)


class IndelRealign(helix.jobs.Job):
    pass

class QualityRecalibration(helix.jobs.Job):
    pass


    
