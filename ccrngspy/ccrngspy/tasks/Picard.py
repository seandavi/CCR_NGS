#!/usr/bin/env python
"""
Originally written by Kelly Vincent
pretty output and additional picard wrappers by Ross Lazarus for rgenetics
Runs all available wrapped Picard tools.
usage: picard_wrapper.py [options]
code Ross wrote licensed under the LGPL
see http://www.gnu.org/copyleft/lesser.html

Re-working of the code by Kenny Daily, Feb. 2012
Converted to using argparse and sub-parsers, have kept PicardBase intact so far.

"""

import argparse
import os
import sys
import subprocess
import tempfile
import shutil
import time
import logging

from ccrngspy import utils

logger = utils.make_local_logger("Picard logging", level="debug", color="green")

galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy %s tool output - see http://getgalaxy.org/" />
<title></title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
"""
galhtmlattr = """Galaxy tool %s run at %s</b><br/>"""
galhtmlpostfix = """</div></body></html>\n"""


def stop_err(msg):
    sys.stderr.write('%s\n' % msg)
    sys.exit()
    

def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))


class PicardBase():
    """
    simple base class with some utilities for Picard
    adapted and merged with Kelly Vincent's code april 2011 Ross
    lots of changes...
    """
    
    def __init__(self, opts=None, arg0=None):
        """ common stuff needed at init for a picard tool
        """

        # Only if opts is passed in do we want to set up
        # Without opts, we can use the argparse(...) method to set up the necessary arguments
        if opts:
            self.opts = opts

            if self.opts.outdir == None:
                self.opts.outdir = os.getcwd() # fixmate has no html file eg so use temp dir
            assert self.opts.outdir <> None,'## PicardBase needs a temp directory if no output directory passed in'

            self.picname = self.baseName(opts.jar)
            self.picname = self.baseName(opts.jar)

            if self.picname.startswith('picard'):
                self.picname = opts.picard_cmd # special case for some tools like replaceheader?

            self.progname = self.baseName(arg0)
            self.version = '0.002'
            self.delme = [] # list of files to destroy
            self.title = opts.title
            self.inputfile = opts.input
            try:
                os.makedirs(opts.outdir)
            except:
                pass
            try:
                os.makedirs(opts.tmpdir)
            except:
                pass
            self.log_filename = os.path.join(self.opts.outdir,'%s.log' % self.picname)
            self.metricsOut =  os.path.join(opts.outdir,'%s.metrics.txt' % self.picname)

            # removed to use our logging code - logger is set at the module level
            # self.setLogging(logfname=self.log_filename)
 
    def baseName(self,name=None):
        return os.path.splitext(os.path.basename(name))[0]

    # Removed to use our logging code - logger is set at the module level
    # def setLogging(self,logfname="picard_wrapper.log"):
    #     """setup a logger
    #     """
    #     logging.basicConfig(level=logging.INFO,
    #                 filename=logfname,
    #                 filemode='a')


    def readLarge(self,fname=None):
        """ read a potentially huge file.
        """
        try:
            # get stderr, allowing for case where it's very large
            tmp = open(fname, 'rb')
            s = ''
            buffsize = 1048576
            try:
                while True:
                    more = tmp.read(buffsize)
                    if len(more) > 0:
                        s += more
                    else:
                        break
            except OverflowError:
                pass
            tmp.close()
        except Exception, e:
            stop_err('Read Large Exception : %s' % str(e))   
        return s

    def constructCL(self, cl=None, output_dir=None):
        """Constrct the command line

        """

        assert cl <> None, 'PicardBase runCL needs a command line as cl'

        if output_dir == None:
            output_dir = self.opts.outdir
        if type(cl) == type([]):
            cl = ' '.join(cl)

        return cl
    
    def runCL(self, cl=None, output_dir=None):
        """Run a command line.
        we have galaxy's temp path as opt.temp_dir so don't really need isolation
        sometimes stdout is needed as the output - ugly hacks to deal with potentially vast artifacts

        """

        # stdout and stderr redirected to these
        # fd,templog = tempfile.mkstemp(dir=output_dir,suffix='rgtempRun.txt')
        # tlf = open(templog,'wb')
        # fd,temperr = tempfile.mkstemp(dir=output_dir,suffix='rgtempErr.txt')
        # tef = open(temperr,'wb')

        cl = self.constructCL(cl=cl, output_dir=output_dir)
        
        process = subprocess.Popen(cl, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, cwd=output_dir)
        rval = process.wait()
        # tlf.close()
        # tef.close()
        # stderrs = self.readLarge(temperr)
        # stdouts = self.readLarge(templog)        

        if rval > 0:
            s = '## executing %s returned status %d and stderr: \n%s\n' % (cl,rval,stderrs)
            stdouts = '%s\n%s' % (stdouts,stderrs)
        else:
            s = '## executing %s returned status %d and nothing on stderr\n' % (cl,rval)

        # logging.info(s)
        # os.unlink(templog) # always
        # os.unlink(temperr) # always
        return s, stdouts, rval  # sometimes s is an output
    
    def runPic(self, jar, cl):
        """
        cl should be everything after the jar file name in the command
        """
        runme = ['java -Xmx%s' % self.opts.maxjheap]
        runme.append(" -Djava.io.tmpdir='%s' " % self.opts.tmpdir)
        runme.append('-jar %s' % jar)
        runme += cl
        s,stdouts,rval = self.runCL(cl=runme, output_dir=self.opts.outdir)
        return stdouts,rval

    def samToBam(self,infile=None,outdir=None):
        """
        use samtools view to convert sam to bam
        """
        fd,tempbam = tempfile.mkstemp(dir=outdir,suffix='rgutilsTemp.bam')
        cl = ['samtools view -h -b -S -o ',tempbam,infile]
        tlog,stdouts,rval = self.runCL(cl,outdir)
        return tlog,tempbam,rval

    def sortSam(self, infile=None,outfile=None,outdir=None):
        """
        """
        print '## sortSam got infile=%s,outfile=%s,outdir=%s' % (infile,outfile,outdir)
        cl = ['samtools sort',infile,outfile]
        tlog,stdouts,rval = self.runCL(cl,outdir)
        return tlog

    def cleanup(self):
        for fname in self.delme:
            try:
                os.unlink(fname)
            except:
                pass
                    
    def prettyPicout(self,transpose,maxrows):
        """organize picard outpouts into a report html page
        """
        res = []
        try:
            r = open(self.metricsOut,'r').readlines()
        except:
            r = []        
        if len(r) > 0:
            res.append('<b>Picard on line resources</b><ul>\n')
            res.append('<li><a href="http://picard.sourceforge.net/index.shtml">Click here for Picard Documentation</a></li>\n')
            res.append('<li><a href="http://picard.sourceforge.net/picard-metric-definitions.shtml">Click here for Picard Metrics definitions</a></li></ul><hr/>\n')
            if transpose:
                res.append('<b>Picard output (transposed to make it easier to see)</b><hr/>\n')       
            else:
                res.append('<b>Picard output</b><hr/>\n')  
            res.append('<table cellpadding="3" >\n')
            dat = []
            heads = []
            lastr = len(r) - 1
            # special case for estimate library complexity hist
            thist = False
            for i,row in enumerate(r):
                if row.strip() > '':
                    srow = row.split('\t')
                    if row.startswith('#'):
                        heads.append(row.strip()) # want strings
                    else:
                        dat.append(srow) # want lists
                    if row.startswith('## HISTOGRAM'):
                        thist = True
            if len(heads) > 0:
                hres = ['<tr class="d%d"><td colspan="2">%s</td></tr>' % (i % 2,x) for i,x in enumerate(heads)]
                res += hres
                heads = []
            if len(dat) > 0:
                if transpose and not thist:
                    tdat = map(None,*dat) # transpose an arbitrary list of lists
                    tdat = ['<tr class="d%d"><td>%s</td><td>%s&nbsp;</td></tr>\n' % ((i+len(heads)) % 2,x[0],x[1]) for i,x in enumerate(tdat)] 
                else:
                    tdat = ['\t'.join(x).strip() for x in dat] # back to strings :(
                    tdat = ['<tr class="d%d"><td colspan="2">%s</td></tr>\n' % ((i+len(heads)) % 2,x) for i,x in enumerate(tdat)]
                res += tdat
                dat = []
            res.append('</table>\n')   
        return res

    def fixPicardOutputs(self,transpose,maxloglines):
        """
        picard produces long hard to read tab header files
        make them available but present them transposed for readability
        """
        logging.shutdown()
        self.cleanup() # remove temp files stored in delme
        rstyle="""<style type="text/css">
        tr.d0 td {background-color: oldlace; color: black;}
        tr.d1 td {background-color: aliceblue; color: black;}
        </style>"""    
        res = [rstyle,]
        res.append(galhtmlprefix % self.progname)   
        res.append(galhtmlattr % (self.picname,timenow()))
        flist = [x for x in os.listdir(self.opts.outdir) if not x.startswith('.')] 
        pdflist = [x for x in flist if os.path.splitext(x)[-1].lower() == '.pdf']
        if len(pdflist) > 0: # assumes all pdfs come with thumbnail .jpgs
            for p in pdflist:
                pbase = os.path.splitext(p)[0] # removes .pdf
                imghref = '%s.jpg' % pbase
                mimghref = '%s-0.jpg' % pbase # multiple pages pdf -> multiple thumbnails without asking!
                if mimghref in flist:
                    imghref=mimghref # only one for thumbnail...it's a multi page pdf
                res.append('<table cellpadding="10"><tr><td>\n')
                res.append('<a href="%s"><img src="%s" title="Click image preview for a print quality PDF version" hspace="10" align="middle"></a>\n' % (p,imghref)) 
                res.append('</tr></td></table>\n')   
        if len(flist) > 0:
            res.append('<b>The following output files were created (click the filename to view/download a copy):</b><hr/>')
            res.append('<table>\n')
            for i,f in enumerate(flist):
                fn = os.path.split(f)[-1]
                res.append('<tr><td><a href="%s">%s</a></td></tr>\n' % (fn,fn))
            res.append('</table><p/>\n') 
        pres = self.prettyPicout(transpose,maxloglines)
        if len(pres) > 0:
            res += pres
        l = open(self.log_filename,'r').readlines()
        llen = len(l)
        if llen > 0: 
            res.append('<b>Picard Tool Run Log</b><hr/>\n') 
            rlog = ['<pre>',]
            if llen > maxloglines:
                n = min(50,int(maxloglines/2))
                rlog += l[:n]
                rlog.append('------------ ## %d rows deleted ## --------------\n' % (llen-maxloglines))
                rlog += l[-n:]
            else:
                rlog += l
            rlog.append('</pre>')
            if llen > maxloglines:
                rlog.append('\n<b>## WARNING - %d log lines truncated - <a href="%s">%s</a> contains entire output</b>' % (llen - maxloglines,self.log_filename,self.log_filename))
            res += rlog
        else:
            res.append("### Odd, Picard left no log file %s - must have really barfed badly?\n" % self.log_filename)
        res.append('<hr/>The freely available <a href="http://picard.sourceforge.net/command-line-overview.shtml">Picard software</a> \n') 
        res.append('generated all outputs reported here running as a <a href="http://getgalaxy.org">Galaxy</a> tool')   
        res.append(galhtmlpostfix) 
        outf = open(self.opts.htmlout,'w')
        outf.write(''.join(res))   
        outf.write('\n')
        outf.close()

    def makePicInterval(self,inbed=None,outf=None):
        """
        picard wants bait and target files to have the same header length as the incoming bam/sam 
        a meaningful (ie accurate) representation will fail because of this - so this hack
        it would be far better to be able to supply the original bed untouched
        Additional checking added Ross Lazarus Dec 2011 to deal with two 'bug' reports on the list
        """
        assert inbed <> None
        bed = open(inbed,'r').readlines()
        sbed = [x.split('\t') for x in bed] # lengths MUST be 5
        lens = [len(x) for x in sbed]
        strands = [x[3] for x in sbed if not x[3] in ['+','-']]
        maxl = max(lens)
        minl = min(lens)
        e = []
        if maxl <> minl:
            e.append("## Input error: Inconsistent field count in %s - please read the documentation on bait/target format requirements, fix and try again" % inbed)
        if maxl <> 5:
            e.append("## Input error: %d fields found in %s, 5 required - please read the warning and documentation on bait/target format requirements, fix and try again" % (maxl,inbed))
        if len(strands) > 0:
            e.append("## Input error: Fourth column in %s is not the required strand (+ or -) - please read the warning and documentation on bait/target format requirements, fix and try again" % (inbed))
        if len(e) > 0: # write to stderr and quit
            print >> sys.stderr, '\n'.join(e)
            sys.exit(1)
        thead = os.path.join(self.opts.outdir,'tempSamHead.txt')
        if self.opts.datatype == 'sam':
            cl = ['samtools view -H -S',self.opts.input,'>',thead]
        else:
            cl = ['samtools view -H',self.opts.input,'>',thead]
        self.runCL(cl=cl,output_dir=self.opts.outdir)
        head = open(thead,'r').readlines()
        s = '## got %d rows of header\n' % (len(head))
        logging.info(s)
        o = open(outf,'w')
        o.write(''.join(head))
        o.write(''.join(bed))
        o.close()
        return outf

    def cleanSam(self, insam=None, newsam=None, picardErrors=[],outformat=None):
        """
        interesting problem - if paired, must remove mate pair of errors too or we have a new set of errors after cleaning - missing mate pairs!
        Do the work of removing all the error sequences
        pysam is cool
        infile = pysam.Samfile("-", "r")
        outfile = pysam.Samfile("-", "w", template = infile)
        for s in infile: outfile.write(s)

        errors from ValidateSameFile.jar look like
        WARNING: Record 32, Read name SRR006041.1202260, NM tag (nucleotide differences) is missing
        ERROR: Record 33, Read name SRR006041.1042721, Empty sequence dictionary.
        ERROR: Record 33, Read name SRR006041.1042721, RG ID on SAMRecord not found in header: SRR006041

        """
        assert os.path.isfile(insam), 'rgPicardValidate cleansam needs an input sam file - cannot find %s' % insam
        assert newsam <> None, 'rgPicardValidate cleansam needs an output new sam file path'
        removeNames = [x.split(',')[1].replace(' Read name ','') for x in picardErrors if len(x.split(',')) > 2]
        remDict = dict(zip(removeNames,range(len(removeNames))))
        infile = pysam.Samfile(insam,'rb')
        info = 'found %d error sequences in picardErrors, %d unique' % (len(removeNames),len(remDict))
        if len(removeNames) > 0:
            outfile = pysam.Samfile(newsam,'wb',template=infile) # template must be an open file
            i = 0
            j = 0
            for row in infile:
                dropme = remDict.get(row.qname,None) # keep if None
                if not dropme:
                    outfile.write(row)
                    j += 1
                else: # discard
                    i += 1
            info = '%s\n%s' % (info, 'Discarded %d lines writing %d to %s from %s' % (i,j,newsam,insam))
            outfile.close()
            infile.close()
        else: # we really want a nullop or a simple pointer copy
            infile.close()
            if newsam:
                shutil.copy(insam,newsam)
        logging.info(info)

    def argparse(self, parser):
        # All tools
        parser.add_argument('-i', '--input', dest='input', help='Input SAM or BAM file', type=str)
        parser.add_argument('-e', '--inputext', default=None, type=str)
        parser.add_argument('-o', '--output', default=None, type=str)
        parser.add_argument('-n', '--title', default="Pick a Picard Tool", type=str)
        parser.add_argument('-t', '--htmlout', default=None)
        parser.add_argument('-d', '--outdir', default=None)
        parser.add_argument('-x', '--maxjheap', default='4g')
        parser.add_argument('-b', '--bisulphite', default='false')
        parser.add_argument('-s', '--sortorder', default='query')     
        parser.add_argument('--tmpdir', default='/tmp')
        parser.add_argument('-j', '--jar', default='')    
        # parser.add_argument('--picard-cmd', default=None)

        # Many tools
        parser.add_argument('--output_format', dest='output_format', help='Output format', default="sam")
        parser.add_argument('--bai_file', dest='bai_file', help='The path to the index file for the input bam file')
        parser.add_argument('--ref', dest='ref', help='Built-in reference with fasta and dict file', default=None)
        parser.add_argument('--assumesorted', default='True')
        parser.add_argument('--readregex', default="[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*")
        parser.add_argument('--ref_file', dest='ref_file', help='Fasta to use as reference', default=None)

        # Need subparsers for the various commands
        subparsers = parser.add_subparsers(help="Picard sub-command help", dest='subparser_name')

        # # CreateSequenceDictionary
        # createseqdictparser = subparsers.add_parser("CreateSequenceDictionary", help="CreateSequenceDictionary help")
        # createseqdictparser.add_argument('--species-name', dest='species_name', help='Species name to use in creating dict file from fasta file')
        # createseqdictparser.add_argument('--build-name', dest='build_name', help='Name of genome assembly to use in creating dict file from fasta file')
        # createseqdictparser.add_argument('--trunc-names', dest='trunc_names', help='Truncate sequence names at first whitespace from fasta file')
        # createseqdictparser.set_defaults(func=create_sequence_dictionary)

        # MarkDuplicates
        markdupparser = subparsers.add_parser("MarkDuplicates", help="MarkDuplicates help")
        markdupparser.add_argument('--remdups', default='true', help='Remove duplicates from output file')
        markdupparser.add_argument('--optdupdist', default="100", help='Maximum pixels between two identical sequences in order to consider them optical duplicates.')
        markdupparser.set_defaults(func=mark_duplicates)

        # CollectRnaSeqMetrics
        collectrnaseqmetricsparser = subparsers.add_parser("CollectRnaSeqMetrics", help="CollectRNASeqMetrics help")
        collectrnaseqmetricsparser.add_argument('--ribosomal_intervals', help='Location of ribosomal sequences in genome.', default=None)
        collectrnaseqmetricsparser.add_argument('--minimum_length', type=int, help='Minimum length [default: %(default)s]', default=500)
        collectrnaseqmetricsparser.add_argument('--chart_output', help='Output of PDF file', default=None)
        collectrnaseqmetricsparser.add_argument('--ignore_sequence', help='Ignore this sequence', default=None)
        collectrnaseqmetricsparser.add_argument('--rrna_fragment_precentage', type=float, help='rRNA fragment precentage. [default: %(default)s]', default=0.8)
        collectrnaseqmetricsparser.add_argument('--metric_accumulation_level', type=str, help='The level(s) at which to accumulate metrics. [default: %(default)s]',
                                                choices=["ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP"], default="SAMPLE")
        collectrnaseqmetricsparser.add_argument('--stop_after', type=int, help='Stop after N reads [default: %(default)s]', default=0)
        collectrnaseqmetricsparser.set_defaults(func=collect_rnaseq_metrics)

        return parser

    def set_options(self, args):
        """Use args from argparse.parse_args to populate the class.

        """

        self.__init__(opts=args)


def setup(args):
    """Do things that all functions may require.
    
    """

    assert args.input <> None

    args.tmp_dir = args.outdir
    args.haveTempout = False # we use this where sam output is an option

    args.rval = 0
    args.stdouts = 'Not run yet'

    pic = PicardBase(args, args.subparser_name)

    # define temporary output
    # if output is sam, it must have that extension, otherwise bam will be produced
    # specify sam or bam file with extension
    if args.output_format == 'sam':
        args.suff = '.sam'
    else:
        args.suff = ''

    args.tmp_fd, args.tempout = tempfile.mkstemp(dir=args.tmpdir, suffix=args.suff)
    
    cl = ['VALIDATION_STRINGENCY=LENIENT',]

    return (pic, args, cl)

def setup_and_cleanup(fxn):
    def wrapper(args):

        (pic, args, cl) = setup(args)
        args = fxn(args, pic, cl)
        cleanup(args, pic)

    return wrapper

def cleanup(args, pic):
    """Finish and clean up, see if the command failed.
    
    """
    
    if args.haveTempout:
        # Some Picard tools produced a potentially intermediate bam file. 
        # Either just move to final location or create sam
        if os.path.exists(args.tempout):
            shutil.move(args.tempout, os.path.abspath(args.output))

    # if (args.htmlout <> None):# or doFix: # return a pretty html page
    #     pic.fixPicardOutputs(transpose=doTranspose,maxloglines=maxloglines)
            
    if args.rval <> 0:
        print >> sys.stderr, '## exit code=%d; stdout=%s' % (args.rval, args.stdouts)
        # signal failure

@setup_and_cleanup
def mark_duplicates(args, pic, cl):

    # (pic, rval, stdouts, cl) = setup(args)
    
    # assume sorted even if header says otherwise
    cl.append('ASSUME_SORTED=%s' % (args.assumesorted))
    # input
    cl.append('INPUT=%s' % args.input)
    # outputs
    cl.append('OUTPUT=%s' % args.output) 
    cl.append('METRICS_FILE=%s' % pic.metricsOut)
    # remove or mark duplicates
    cl.append('REMOVE_DUPLICATES=%s' % args.remdups)
    # the regular expression to be used to parse reads in incoming SAM file
    cl.append('READ_NAME_REGEX="%s"' % args.readregex)
    # maximum offset between two duplicate clusters
    cl.append('OPTICAL_DUPLICATE_PIXEL_DISTANCE=%s' % args.optdupdist)
    args.stdouts, args.rval = pic.runPic(args.jar, cl)
    return args

    # cleanup(args, pic, rval, stdouts)

# @setup_and_cleanup
# def create_sequence_dictionary(args, pic, cl):

#     assert args.ref_file, "Did not specify a ref file"
#     assert args.ref, "Did not specify ref"

#     csd = 'CreateSequenceDictionary'
#     realjarpath = os.path.split(args.jar)[0]
#     jarpath = os.path.join(realjarpath,'%s.jar' % csd) # for refseq

#     tmp_ref_fd, tmp_ref_name = tempfile.mkstemp(dir=args.tmpdir , prefix = pic.picname)
#     args.ref = '%s.fasta' % tmp_ref_name

#     # build dict
#     dict_file_name = '%s.dict' % tmp_ref_name
#     os.symlink(args.ref_file, args.ref)

#     cl = ['REFERENCE=%s' % args.ref]
#     cl.append('OUTPUT=%s' % dict_file_name)
#     cl.append('URI=%s' % os.path.basename(args.ref_file))
#     cl.append('TRUNCATE_NAMES_AT_WHITESPACE=%s' % args.trunc_names)

#     if args.species_name:
#         cl.append('SPECIES=%s' % args.species_name)

#     if args.build_name:
#         cl.append('GENOME_ASSEMBLY=%s' % args.build_name)

#     pic.delme.append(dict_file_name)
#     pic.delme.append(args.ref)
#     pic.delme.append(tmp_ref_name)

#     stdouts,rval = pic.runPic(jarpath, cl)

#     return stdouts, rval

#     # run relevant command(s)

@setup_and_cleanup
def collect_rnaseq_metrics(args, pic, cl):
    """Run CollectRNASeqMetrics tool.

    http://picard.sourceforge.net/command-line-overview.shtml#CollectRnaSeqMetrics
    
    """

    assert args.ref_flat, "Did not specify a ref_flat file, which is required."
    assert args.output, "Did not specify an output file, which is required."

    # inputs
    cl.append('INPUT=%s' % args.input)

    # Gene annotations in refFlat format.
    cl.append('REF_FLAT=%s' % (args.ref_flat))

    # Gene annotations in refFlat format.
    cl.append('REFERENCE_SEQUENCE=%s' % (args.ref_file))

    # locations of rRNA seqs in genome in interval_list format.
    cl.append('RIBOSOMAL_INTERVALS=%s' % args.ribosomal_intervals)

    cl.append('MINIMUM_LENGTH=%i' % args.minimum_length)
    cl.append('CHART_OUTPUT=%s' % args.chart_output)
    cl.append('IGNORE_SEQUENCE=%s' % args.ignore_sequence)

    cl.append('RRNA_FRAGMENT_PERCENTAGE=%d' % args.rrna_fragment_percentage)
    cl.append('METRIC_ACCUMULATION_LEVEL=%s' % args.metric_accumulation_level)

    cl.append('ASSUME_SORTED=%s' % (args.assumesorted.lower()))

    cl.append('STOP_AFTER=%i' % (args.stop_after))
    
    # outputs
    cl.append('OUTPUT=%s' % args.output) 

    args.stdouts, args.rval = pic.runPic(args.jar, cl)
    return args


        
    
def __main__():

    doFix = False # tools returning htmlfile don't need this
    doTranspose = True # default
    maxloglines = 100 # default 

    parser = argparse.ArgumentParser()

    # Add Picard arguments
    picard = PicardBase()
    parser = picard.argparse(parser)

    args = parser.parse_args()
    args.func(args)

    # # CollectInsertSizeMetrics
    # parser.add_argument('', '--taillimit', default="0")
    # parser.add_argument('', '--histwidth', default="0")
    # parser.add_argument('', '--minpct', default="0.01")
    # parser.add_argument('', '--malevel', default='')
    # parser.add_argument('', '--deviations', default="0.0")

    # # CollectAlignmentSummaryMetrics
    # parser.add_argument('', '--maxinsert', default="20")
    # parser.add_argument('', '--adaptors', default='')

    # # FixMateInformation and validate
    # # CollectGcBiasMetrics
    # parser.add_argument('', '--windowsize', default='100')
    # parser.add_argument('', '--mingenomefrac', default='0.00001')    

    # # AddOrReplaceReadGroups
    # parser.add_argument('', '--rg-opts', dest='rg_opts', help='Specify extra (optional) arguments with full, otherwise preSet')
    # parser.add_argument('', '--rg-lb', dest='rg_library', help='Read Group Library')
    # parser.add_argument('', '--rg-pl', dest='rg_platform', help='Read Group platform (e.g. illumina, solid)')
    # parser.add_argument('', '--rg-pu', dest='rg_plat_unit', help='Read Group platform unit (eg. run barcode) ')
    # parser.add_argument('', '--rg-sm', dest='rg_sample', help='Read Group sample name')
    # parser.add_argument('', '--rg-id', dest='rg_id', help='Read Group ID')
    # parser.add_argument('', '--rg-cn', dest='rg_seq_center', help='Read Group sequencing center name')
    # parser.add_argument('', '--rg-ds', dest='rg_desc', help='Read Group description')

    # # ReorderSam
    # parser.add_argument('', '--allow-inc-dict-concord', dest='allow_inc_dict_concord', help='Allow incomplete dict concordance')
    # parser.add_argument('', '--allow-contig-len-discord', dest='allow_contig_len_discord', help='Allow contig length discordance')

    # # ReplaceSamHeader
    # repsamheaderparser = subparsers.add_parser("ReplaceSamHeader", help="ReplaceSamHeader help")
    # repsamheaderparser.add_argument('--header-file', dest='header_file', help='sam or bam file from which header will be read')

    # #estimatelibrarycomplexity
    # parser.add_argument('','--minid', default="5")
    # parser.add_argument('','--maxdiff', default="0.03")
    # parser.add_argument('','--minmeanq', default="20")

    # #hsmetrics
    # parser.add_argument('','--baitbed', default=None)
    # parser.add_argument('','--targetbed', default=None)

    # #validate
    # parser.add_argument('','--ignoreflags', action='append', type="string")
    # parser.add_argument('','--maxerrors', default=None)
    # parser.add_argument('','--datatype', default=None)
    # parser.add_argument('','--bamout', default=None)
    # parser.add_argument('','--samout', default=None)

    # args.sortme = (args.assumesorted == 'false')

    # need to add
    # instance that does all the work



    # if pic.picname == 'AddOrReplaceReadGroups':
    #     # sort order to match Galaxy's default
    #     cl.append('SORT_ORDER=coordinate')
    #     # input
    #     cl.append('INPUT=%s' % args.input)
    #     # outputs
    #     cl.append('OUTPUT=%s' % tempout)
    #     # required read groups
    #     cl.append('RGLB="%s"' % args.rg_library)
    #     cl.append('RGPL="%s"' % args.rg_platform)
    #     cl.append('RGPU="%s"' % args.rg_plat_unit)
    #     cl.append('RGSM="%s"' % args.rg_sample)
    #     if args.rg_id:
    #         cl.append('RGID="%s"' % args.rg_id)
    #     # optional read groups
    #     if args.rg_seq_center:
    #         cl.append('RGCN="%s"' % args.rg_seq_center)
    #     if args.rg_desc:
    #         cl.append('RGDS="%s"' % args.rg_desc)
    #     stdouts,rval = pic.runPic(args.jar, cl)
    #     haveTempout = True

  #   elif pic.picname == 'BamIndexStats':
  #       tmp_fd, tmp_name = tempfile.mkstemp(dir=tmp_dir)
  #       tmp_bam_name = '%s.bam' % tmp_name
  #       tmp_bai_name = '%s.bai' % tmp_bam_name
  #       os.symlink(opts.input, tmp_bam_name)
  #       os.symlink(opts.bai_file, tmp_bai_name)
  #       cl.append('INPUT=%s' % (tmp_bam_name))
  #       pic.delme.append(tmp_bam_name)
  #       pic.delme.append(tmp_bai_name)
  #       pic.delme.append(tmp_name)
  #       stdouts,rval = pic.runPic(opts.jar, cl)
  #       f = open(pic.metricsOut,'a')
  #       f.write(stdouts) # got this on stdout from runCl
  #       f.write('\n')
  #       f.close()
  #       doTranspose = False # but not transposed

  #   elif pic.picname == 'EstimateLibraryComplexity':
  #       cl.append('I=%s' % opts.input)
  #       cl.append('O=%s' % pic.metricsOut)
  #       if float(opts.minid) > 0:
  #           cl.append('MIN_IDENTICAL_BASES=%s' % opts.minid)
  #       if float(opts.maxdiff) > 0.0:
  #           cl.append('MAX_DIFF_RATE=%s' % opts.maxdiff)
  #       if float(opts.minmeanq) > 0:
  #           cl.append('MIN_MEAN_QUALITY=%s' % opts.minmeanq)
  #       if opts.readregex > '':
  #           cl.append('READ_NAME_REGEX="%s"' % opts.readregex)
  #       if float(opts.optdupdist) > 0:
  #           cl.append('OPTICAL_DUPLICATE_PIXEL_DISTANCE=%s' % opts.optdupdist)
  #       stdouts,rval = pic.runPic(opts.jar, cl)

  #   elif pic.picname == 'CollectAlignmentSummaryMetrics':
  #       # Why do we do this fakefasta thing? 
  #       # Because we need NO fai to be available or picard barfs unless it matches the input data.
  #       # why? Dunno Seems to work without complaining if the .bai file is AWOL....
  #       fakefasta = os.path.join(opts.outdir,'%s_fake.fasta' % os.path.basename(args.ref))
  #       try:
  #           os.symlink(args.ref,fakefasta)
  #       except:
  #           s = '## unable to symlink %s to %s - different devices? Will shutil.copy'
  #           info = s
  #           shutil.copy(args.ref,fakefasta)
  #       pic.delme.append(fakefasta)
  #       cl.append('ASSUME_SORTED=true')
  #       adaptlist = opts.adaptors.split(',')
  #       adaptorseqs = ['ADAPTER_SEQUENCE=%s' % x for x in adaptlist]
  #       cl += adaptorseqs
  #       cl.append('IS_BISULFITE_SEQUENCED=%s' % opts.bisulphite)
  #       cl.append('MAX_INSERT_SIZE=%s' % opts.maxinsert)
  #       cl.append('OUTPUT=%s' % pic.metricsOut)
  #       cl.append('R=%s' % fakefasta)
  #       cl.append('TMP_DIR=%s' % opts.tmpdir)
  #       if not opts.assumesorted.lower() == 'true': # we need to sort input
  #           sortedfile = '%s.sorted' % os.path.basename(opts.input)
  #           if opts.datatype == 'sam': # need to work with a bam 
  #               tlog,tempbam,trval = pic.samToBam(opts.input,opts.outdir)
  #               pic.delme.append(tempbam)
  #               try:
  #                   tlog = pic.sortSam(tempbam,sortedfile,opts.outdir)
  #               except:
  #                   print '## exception on sorting sam file %s' % opts.input
  #           else: # is already bam
  #               try:
  #                   tlog = pic.sortSam(opts.input,sortedfile,opts.outdir)
  #               except : # bug - [bam_sort_core] not being ignored - TODO fixme
  #                   print '## exception %s on sorting bam file %s' % (sys.exc_info()[0],opts.input)
  #           cl.append('INPUT=%s.bam' % os.path.abspath(os.path.join(opts.outdir,sortedfile)))
  #           pic.delme.append(os.path.join(opts.outdir,sortedfile))
  #       else:
  #           cl.append('INPUT=%s' % os.path.abspath(opts.input)) 
  #       stdouts,rval = pic.runPic(opts.jar, cl)
       
        
  #   elif pic.picname == 'CollectGcBiasMetrics':
  #       assert os.path.isfile(args.ref),'PicardGC needs a reference sequence - cannot read %s' % args.ref
  #       # sigh. Why do we do this fakefasta thing? Because we need NO fai to be available or picard barfs unless it has the same length as the input data.
  #       # why? Dunno 
  #       fakefasta = os.path.join(opts.outdir,'%s_fake.fasta' % os.path.basename(args.ref))
  #       try:
  #           os.symlink(args.ref,fakefasta)
  #       except:
  #           s = '## unable to symlink %s to %s - different devices? May need to replace with shutil.copy'
  #           info = s
  #           shutil.copy(args.ref,fakefasta)
  #       pic.delme.append(fakefasta)
  #       x = 'rgPicardGCBiasMetrics'
  #       pdfname = '%s.pdf' % x
  #       jpgname = '%s.jpg' % x
  #       tempout = os.path.join(opts.outdir,'rgPicardGCBiasMetrics.out')
  #       temppdf = os.path.join(opts.outdir,pdfname)
  #       cl.append('R=%s' % fakefasta)
  #       cl.append('WINDOW_SIZE=%s' % opts.windowsize)
  #       cl.append('MINIMUM_GENOME_FRACTION=%s' % opts.mingenomefrac)
  #       cl.append('INPUT=%s' % opts.input)
  #       cl.append('OUTPUT=%s' % tempout)
  #       cl.append('TMP_DIR=%s' % opts.tmpdir)
  #       cl.append('CHART_OUTPUT=%s' % temppdf)
  #       cl.append('SUMMARY_OUTPUT=%s' % pic.metricsOut)
  #       stdouts,rval = pic.runPic(opts.jar, cl)
  #       if os.path.isfile(temppdf):
  #           cl2 = ['convert','-resize x400',temppdf,os.path.join(opts.outdir,jpgname)] # make the jpg for fixPicardOutputs to find
  #           s,stdouts,rval = pic.runCL(cl=cl2,output_dir=opts.outdir)
  #       else:
  #           s='### runGC: Unable to find pdf %s - please check the log for the causal problem\n' % temppdf
  #       lf = open(pic.log_filename,'a')
  #       lf.write(s)
  #       lf.write('\n')
  #       lf.close()
        
  #   elif pic.picname == 'CollectInsertSizeMetrics':
  #       """ <command interpreter="python">
  #  picard_wrapper.py -i "$input_file" -n "$out_prefix" --tmpdir "${__new_file_path__}" --deviations "$deviations"
  #  --histwidth "$histWidth" --minpct "$minPct" --malevel "$malevel"
  #  -j "${GALAXY_DATA_INDEX_DIR}/shared/jars/picard/CollectInsertSizeMetrics.jar" -d "$html_file.files_path" -t "$html_file"
  # </command>
  #       """
  #       isPDF = 'InsertSizeHist.pdf'
  #       pdfpath = os.path.join(opts.outdir,isPDF)
  #       histpdf = 'InsertSizeHist.pdf'
  #       cl.append('I=%s' % opts.input)
  #       cl.append('O=%s' % pic.metricsOut)
  #       cl.append('HISTOGRAM_FILE=%s' % histpdf)
  #       #if opts.taillimit <> '0': # this was deprecated although still mentioned in the docs at 1.56
  #       #    cl.append('TAIL_LIMIT=%s' % opts.taillimit)
  #       if  opts.histwidth <> '0':
  #           cl.append('HISTOGRAM_WIDTH=%s' % opts.histwidth)
  #       if float(opts.minpct) > 0.0:
  #           cl.append('MINIMUM_PCT=%s' % opts.minpct)
  #       if float(opts.deviations) > 0.0:
  #           cl.append('DEVIATIONS=%s' % opts.deviations)
  #       if opts.malevel:
  #           malists = opts.malevel.split(',')
  #           malist = ['METRIC_ACCUMULATION_LEVEL=%s' % x for x in malists]
  #           cl += malist
  #       stdouts,rval = pic.runPic(opts.jar, cl)
  #       if os.path.exists(pdfpath): # automake thumbnail - will be added to html 
  #           cl2 = ['mogrify', '-format jpg -resize x400 %s' % pdfpath]
  #           pic.runCL(cl=cl2,output_dir=opts.outdir)
  #       else:
  #           s = 'Unable to find expected pdf file %s<br/>\n' % pdfpath
  #           s += 'This <b>always happens if single ended data was provided</b> to this tool,\n'
  #           s += 'so please double check that your input data really is paired-end NGS data.<br/>\n'
  #           s += 'If your input was paired data this may be a bug worth reporting to the galaxy-bugs list\n<br/>'
  #           logging.info(s)
  #       if len(stdouts) > 0:
  #          logging.info(stdouts)
        
  # elif pic.picname == 'MarkDuplicates':
  #     mark_duplicate(args)


  #   elif pic.picname == 'FixMateInformation':
  #       cl.append('I=%s' % opts.input)
  #       cl.append('O=%s' % tempout)
  #       cl.append('SORT_ORDER=%s' % opts.sortorder)
  #       stdouts,rval = pic.runPic(opts.jar,cl)
  #       haveTempout = True
        
  #   elif pic.picname == 'ReorderSam':
  #       # input
  #       cl.append('INPUT=%s' % opts.input)
  #       # output
  #       cl.append('OUTPUT=%s' % tempout)
  #       # reference
  #       cl.append('REFERENCE=%s' % args.ref)
  #       # incomplete dict concordance
  #       if opts.allow_inc_dict_concord == 'true':
  #           cl.append('ALLOW_INCOMPLETE_DICT_CONCORDANCE=true')
  #       # contig length discordance
  #       if opts.allow_contig_len_discord == 'true':
  #           cl.append('ALLOW_CONTIG_LENGTH_DISCORDANCE=true')
  #       stdouts,rval = pic.runPic(opts.jar, cl)
  #       haveTempout = True

  #   elif pic.picname == 'ReplaceSamHeader':
  #       cl.append('INPUT=%s' % opts.input)
  #       cl.append('OUTPUT=%s' % tempout)
  #       cl.append('HEADER=%s' % opts.header_file)
  #       stdouts,rval = pic.runPic(opts.jar, cl)
  #       haveTempout = True

  #   elif pic.picname == 'CalculateHsMetrics':
  #       maxloglines = 100
  #       baitfname = os.path.join(opts.outdir,'rgPicardHsMetrics.bait')
  #       targetfname = os.path.join(opts.outdir,'rgPicardHsMetrics.target')
  #       baitf = pic.makePicInterval(opts.baitbed,baitfname)
  #       if opts.targetbed == opts.baitbed: # same file sometimes
  #           targetf = baitf
  #       else:
  #           targetf = pic.makePicInterval(opts.targetbed,targetfname)   
  #       cl.append('BAIT_INTERVALS=%s' % baitf)
  #       cl.append('TARGET_INTERVALS=%s' % targetf)
  #       cl.append('INPUT=%s' % os.path.abspath(opts.input))
  #       cl.append('OUTPUT=%s' % pic.metricsOut)
  #       cl.append('TMP_DIR=%s' % opts.tmpdir)
  #       stdouts,rval = pic.runPic(opts.jar,cl)
           
  #   elif pic.picname == 'ValidateSamFile':
  #       import pysam
  #       doTranspose = False
  #       sortedfile = os.path.join(opts.outdir,'rgValidate.sorted')
  #       stf = open(pic.log_filename,'w')
  #       tlog = None
  #       if opts.datatype == 'sam': # need to work with a bam 
  #           tlog,tempbam,rval = pic.samToBam(opts.input,opts.outdir)
  #           try:
  #               tlog = pic.sortSam(tempbam,sortedfile,opts.outdir)
  #           except:
  #               print '## exception on sorting sam file %s' % opts.input
  #       else: # is already bam
  #           try:
  #               tlog = pic.sortSam(opts.input,sortedfile,opts.outdir)
  #           except: # bug - [bam_sort_core] not being ignored - TODO fixme
  #               print '## exception on sorting bam file %s' % opts.input
  #       if tlog:
  #           print '##tlog=',tlog
  #           stf.write(tlog)
  #           stf.write('\n')
  #       sortedfile = '%s.bam' % sortedfile # samtools does that      
  #       cl.append('O=%s' % pic.metricsOut)
  #       cl.append('TMP_DIR=%s' % opts.tmpdir)
  #       cl.append('I=%s' % sortedfile)
  #       opts.maxerrors = '99999999'
  #       cl.append('MAX_OUTPUT=%s' % opts.maxerrors)
  #       if opts.ignoreflags[0] <> 'None': # picard error values to ignore
  #           igs = ['IGNORE=%s' % x for x in opts.ignoreflags if x <> 'None']
  #           cl.append(' '.join(igs))
  #       if opts.bisulphite.lower() <> 'false':
  #           cl.append('IS_BISULFITE_SEQUENCED=true')
  #       if opts.ref <> None or opts.ref_file <> None:
  #           cl.append('R=%s' %  args.ref)
  #       stdouts,rval = pic.runPic(opts.jar,cl)
  #       if opts.datatype == 'sam':
  #           pic.delme.append(tempbam)
  #       newsam = opts.output
  #       outformat = 'bam'
  #       pe = open(pic.metricsOut,'r').readlines()
  #       pic.cleanSam(insam=sortedfile, newsam=newsam, picardErrors=pe,outformat=outformat)
  #       pic.delme.append(sortedfile) # not wanted
  #       stf.close()
  #       pic.cleanup()



if __name__=="__main__":
    __main__()

