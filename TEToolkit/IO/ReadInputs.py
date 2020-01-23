'''
Created on Oct 12, 2011

this code is for reading input parameters.

@author: Ying Jin
@status: 
@contact: yjin@cshl.edu
@version: 
'''
# python modules
import sys
import os
import re
import logging
import time
import gzip
import collections
from math import log

from TEToolkit.Constants import *
from TEToolkit.ShortRead.ParseBEDFile import BEDFile, BAMFile, SAMFile


# Taken from HTSeq
class FileOrSequence(object):
    """ The construcutor takes one argument, which may either be a string,
    which is interpreted as a file name (possibly with path), or a
    connection, by which we mean a text file opened for reading, or
    any other object that can provide an iterator over strings
    (lines of the file).

    The advantage of passing a file name instead of an already opened file
    is that if an iterator is requested several times, the file will be
    re-opened each time. If the file is already open, its lines can be read
    only once, and then, the iterator stays exhausted.

    Furthermore, if a file name is passed that end in ".gz" or ".gzip"
    (case insensitive), it is transparently gunzipped.
    """

    def __init__(self, filename_or_sequence):
        self.fos = filename_or_sequence
        self.line_no = None
      
    def __iter__(self ):
        self.line_no = 1
        if isinstance( self.fos, str ):
            if self.fos.lower().endswith( ( ".gz" , ".gzip" ) ):
                lines = gzip.open( self.fos )
            else:
                lines = open( self.fos )
        else:
            lines = self.fos
        for line in lines:
            yield line
            self.line_no += 1
        if isinstance( self.fos, str ):
            lines.close()
        self.line_no = None
         
    def __repr__( self ):
        if isinstance( self.fos, str ):
            return "<%s object, connected to file name '%s'>" % (
                self.__class__.__name__, self.fos )
        else:
            return "<%s object, connected to %s >" % (
                self.__class__.__name__, repr( self.fos ) )
            
    def get_line_number_string( self ):
        if self.line_no is None:
            if isinstance( self.fos, str ):
                return "file %s closed" % self.fos
            else:
                return "file closed"
        if isinstance( self.fos, str ):
            return "line %d of file %s" % ( self.line_no, self.fos )
        else:
            return "line %d" % self.line_no


# Taken from HTSeq
class SAM_Reader(FileOrSequence):
    """Parse a SAM File"""

    def __iter__(self):
        for line in FileOrSequence.__iter__( self ):
            if line.startswith( "@" ):
                continue
            try:
                (qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, trash) = line.split("\t", 10)
            except ValueError, e:
                e.args = e.args + ( self.get_line_number_string(), )
                raise
            yield (qname, int(flag), rname, int(pos)-1, mapq, cigar, rnext, pnext, seq.upper())


def read_opts3(parser):
    args = parser.parse_args()
    if not os.path.isfile(args.tefile):
        logging.error("No such file: %s !\n" % args.tefile)
        sys.exit(1)
    if not os.path.isfile(args.gtffile):
        logging.error("No such file: %s !\n" % args.gtffile)
        sys.exit(1)
    # Check for BAM file
    if not os.path.isfile(args.bam) :
        logging.error("No such file: %s !\n" % args.bam)
        sys.exit(1)
    # Identify file format for subsequent processing (parsing)
    if args.fileformat == "BAM":
        args.parser = "BAM"
    elif args.fileformat == "SAM":
        args.parser = "SAM"
    else:
        logging.error("Does not support such file format: %s !\n" % args.fileformat)
        sys.exit(1)
    # What type of RNA-Seq experiment (stranded or not)
    if args.stranded not in ['yes', 'no', 'reverse']:
        logging.error("Does not support such stranded value: %s !\n" % args.stranded)
        sys.exit(1)

    # Method of assigning reads to annotation (gene or TE)
    if args.te_mode not in ['uniq', 'multi']:
        logging.error("multi-mapper counting mode %s not supported\n" % args.te_mode)
        parser.print_help()
        sys.exit(1)
    
    if args.sortByPos:
        args.sortByPos = True
    else:
        args.sortByPos = False

    if args.numItr < 0:
        args.numItr = 0
    if args.fragLength < 0:
        logging.error("the fragment length cannot be negative. \n")
        sys.exit(1)
    if args.minL < 0:
        logging.error("the minimum fragment length cannot be negative. \n")
        sys.exit(1)
    if args.maxL < 0:
        logging.error("the maximum fragment length cannot be negative. \n")
        sys.exit(1)

    # Level of logging for tool
    logging.basicConfig(level=(4 - args.verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr, filemode="w")
    
    args.error = logging.critical        # function alias
    args.warn = logging.warning
    args.debug = logging.debug
    args.info = logging.info
    
    args.argtxt = "\n".join(("# ARGUMENTS LIST:", \
        "# name = %s" % (args.prefix), \
        "# BAM file = %s" % (args.bam), \
        "# GTF file = %s " % (args.gtffile), \
        "# TE file = %s " % (args.tefile), \
        "# multi-mapper mode = %s " % (args.te_mode), \
        "# stranded = %s " % (args.stranded), \
        "# number of iteration = %d" % (args.numItr), \
        "# Alignments grouped by read ID = %s\n" % (not args.sortByPos)
    ))
    return args 


def read_chrlen_tbl(chrfile, error, info):
    ''' read in chrom_size file '''
 
    if not os.path.isfile(chrfile) :
        error("No such file: %s !\n" % (chrfile))
        sys.exit(1)
    try:
        f = open(chrfile,'r')
    except IOError :
        error("open file %s error !\n" %(chrfile))
        sys.exit(1)
    else :
        chrlen_map = dict()    # hash table for chromosome length.
        cnt = 0
        for line in f :
            cnt += 1
            line = line.strip()
            items = line.split('\t')    # skip empty line.
            if len(items) < 2 :
                info("Insufficient chromosome information at % s, line: %s. Skip!\n" % (chrfile, line))
            if re.match('^(c|C)hr', items[0]) and re.match('^[0-9]+$', items[1]) :
                chrlen_map[items[0]] = int(items[1])
            else:
                info("Format error at %s, line %d: %s. Skip!\n" % (chrfile,cnt,line))
        f.close()
        
    return chrlen_map


def read_short_reads(samples,parser,TEmode):
    '''read short reads from single or multple samples and stored in short read objects '''
    
    shortReads = []
 #   chroms = chrlen_tbl.keys()
    for i in range(len(samples)) :
        s = samples[i]
        if not os.path.isfile(s) :
            logging.error("No such file %s !\n" %(s))
            sys.exit(1)
        logging.info("reading sample file %s ...\n" %(s))
        time.sleep(1)

        #sbed = __bam2bed(s,0,error)
        #b = BEDFile(sbed,chroms)
        #b = parser(s,chroms)
        b = parser(s)
        t = b.build_fwtrack(TEmode)
        shortReads.append(t)
   
    return shortReads

def read_short_reads_sameFam(samples,parser,teIdx):
    '''read short reads from single or multple samples and stored in short read objects '''
    
    shortReads = []

    for i in range(len(samples)) :
        s = samples[i]
        if not os.path.isfile(s) :
            logging.error("No such file %s !\n" %(s))
            sys.exit(1)
        logging.info("reading sample file %s ...\n" %(s))
        time.sleep(1)

        b = parser(s)
        t = b.build_fwtrack_v2(teIdx)
        shortReads.append(t)
   
    return shortReads

   
def __bam2bed(sample,pairend,error):
    res = sample + ".bed"
    if pairend == 0:  # single end
        try:
            os.system("bamToBED -ed -i sample >res")
            res = __assignWeight(sample,".bed",error)
        except :
            error("file format error %s !\n" %(sample))
            sys.exit(0)
        
    else :
        try:
            os.system("bamToBED -bedpe -i sample >res")
            res = __assignWeight(sample,".bed",error)
        except :
            error("file format error %s !\n" %(sample))
            sys.exit(0)

    return res 


def __assignWeight(sample, suffix, error):

    src = sample + suffix
    dest = sample + ".bal.bed"

    lines = []
    cur_seqid = "-1"
    multi_num = 0
    
    try:
        f = open(src, 'r')
        of = open(dest, 'w')
    except IOError:
        error("open file %s error !\n" % src)
        sys.exit(1)
    else:
        for line in f:
            line = line.strip()
            arr = line.split('\t')
            if cur_seqid == arr[3]:
                lines.append(line)
                multi_num += 1
            else:
        
                if multi_num > 0:
                    val = 1/multi_num
                    for record in lines:
                        of.write(record + "\t" + str(val) + "\n")
                lines.clear()
                lines.append(line)
                cur_seqid = arr[3]
                multi_num = 1
    f.close()

    if multi_num > 0:
        val = 1/multi_num
        for record in lines:
            of.write(record + "\t" + str(val) + "\n")

    of.close()
    return dest
