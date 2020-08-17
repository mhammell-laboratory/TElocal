TElocal
=============

Version: 1.1.1

*NOTE* TElocal relies on specially curated and indexed GTF files, which are not
packaged with this software due to their size. Please go to 
`our website <http://hammelllab.labsites.cshl.edu/software#TElocal>`_
for instructions to download the prebuilt curated indices.

TElocal takes RNA-seq (and similar data) and annotates reads to both
genes & transposable elements.


`Github Page <https://github.com/mhammell-laboratory/TElocal>`_

`MHammell Lab <http://hammelllab.labsites.cshl.edu/software>`_

Created by Ying Jin, Eric Paniagua, Oliver Tam & Molly Hammell, February 2014

Copyright (C) 2014-2020 Talitha Forcier, Ying Jin, Eric Paniagua, Oliver Tam & Molly Hammell

Contact: Talitha Forcier (talitha@cshl.edu)

Requirements
------------

Python:     2.6.x or 2.7.x or 3.x

pysam:      0.9.x or greater


Installation
------------

1. Download compressed tarball.
2. Unpack tarball.
3. Navigate into unpacked directory.
4. Run the following::

    $ python setup.py install

If you want to install locally (e.g. /local/home/usr),
run this command instead::

    $ python setup.py install --prefix /local/home/usr

*NOTE* In the above example, you must add

    /local/home/usr/bin

to the PATH variable, and

     /local/home/usr/lib/pythonX.Y/site-packages

to the PYTHONPATH variable, where X refers to the major python version, and Y refers to the minor python version. (e.g. python2.7 if using python version 2.7.x, and python3.6 if using python version 3.6.x)


TElocal
=======

Usage
-----

::

    usage: TElocal -b alignment-file
                   --GTF genic-annot-file
                   --TE TE-annot-file
                   [optional arguments]

    Required arguments:
      -b | --BAM alignment-file    RNAseq alignment file (BAM preferred)
      --GTF genic-annot-file       GTF or ind file for gene annotations
      --TE TE-annot-file           locInd file for transposable element annotations

    Optional arguments:

      *Input/Output options*
      --stranded [option]   Is this a stranded library? (no, forward, or reverse).
                 no      -  Library is unstranded   
                 forward -  "Second-strand cDNA library (e.g. QIAseq stranded)
                 reverse -  "First-strand" cDNA library (e.g. Illumina TruSeq stranded)
                            DEFAULT: no.
      --sortByPos           Input file is sorted by chromosome position.
      --project [name]      Prefix used for output files (e.g. project name)
                            DEFAULT: TElocal_out

      *Analysis options*
      --mode [TE counting mode]
         How to count TE:
            uniq        (unique mappers only)
            multi       (distribute among all alignments).
         DEFAULT: multi
      -L | --fragmentLength [fragLength]
         Average length of fragment used for single-end sequencing
         DEFAULT: For paired-end, estimated from the input alignment file. For single-end, ignored by default.
      -i | --iteration 
         maximum number of iterations used to optimize multi-reads assignment. DEFAULT: 100

      *Other options*
      -h | --help
         Show help message
      --verbose [number]
         Set verbose level.
           0: only show critical messages
           1: show additional warning messages
           2: show process information
           3: show debug messages
         DEFAULT: 2
      --version
         Show program's version and exit

*NOTE* BAM files must be either unsorted or sorted by queryname. If the BAM files are sorted by position, please use the '--sortByPos' option


Example Command Lines
---------------------

If BAM files are unsorted, or sorted by queryname:: 

    TElocal -b RNAseq.bam --GTF gene_annots.gtf --TE te_annots.locInd --project sample_nosort_test

If BAM files are sorted by coordinates/position::

    TElocal --sortByPos -b RNAseq.bam --GTF gene_annots.gtf --TE te_annots.locInd --project sample_sorted_test

Cluster Usage Recommendations
-----------------------------

In our experience, we recommend around 20-30Gb of memory for analyzing human samples (hg19) with around 20-30 million mapped reads when running on a cluster.


Recommendations for TElocal input files
=============================================

TElocal can perform transposable element quantification from alignment results (e.g. BAM files) generated from a variety of programs. 
Given the variety of experimental systems, we could not provide an optimal alignment strategy for every approach. Therefore,
we recommend that users identify the optimal parameters for their particular genome and alignment program in order to get the best
results.

When optimizing the alignment parameters, we recommend taking these points into consideration:

*Allowing sufficient number of multi-mappers during alignment*

Most alignment programs provide only 1 alignment per read by default. We recommend reporting multiple alignments per read. We have found 
that reporting a maximum of 100 alignments per read provides an optimal compromise between the size of the alignment file and recovery 
of multi-mappers in many genome builds. However, we highly suggest that users optimize this parameter for their particular experiment, 
as this could significantly improve the quality of transposable element quantification.

*Paired end sequencing input*

For paired-end libraries, it is recommended that only alignments from properly paired reads are present in the input BAM file. I.e., each read 1 alignment should only have a single read 2 alignment. For example, if read 1 matched 3 genomic locations (A, B, C), then if read 2 also match 3 genomic locations (A', B', C'), then all three pairs of alignments could be used (and should be in the BAM file). However, if alignment C of read 1 was matched with more than one alignment of read 2 (e.g. C' and C*), then alignment C should be discarded (as there are unmatched alignments between read 1 and read 2). `STAR <https://github.com/alexdobin/STAR>`_ only outputs properly paired alignments by default, while `Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ requires the :code:`--no-mixed` parameter to be used.

*Specific recommendations when using STAR*

`STAR <https://github.com/alexdobin/STAR>`_ utilizes two parameters for optimal identification of multi-mappers `--outFilterMultimapNmax` and `--outAnchorMultimapNmax`. 
The author of STAR recommends that `--winAnchorMultimapNmax` should be set at twice the value used in `--outFilterMultimapNmax`, 
but no less than 50. In our study, we used the same number for both parameters (100), and found negligible differences in identifying 
multi-mappers. Upon further discussion with the author of STAR, we recommend that setting the same value for `--winAnchorMultimapNmax`
and `--outFilterMultimapNmax`, though we highly suggest users test multiple values of `--winAnchorMultimapNmax` to identify the 
optimal value for their experiment.


Copying & distribution
======================

TElocal is part of `TEToolkit suite <http://hammelllab.labsites.cshl.edu/software/>`_.

TElocal is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but *WITHOUT ANY WARRANTY*; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE*.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TElocal.  If not, see `this website <http://www.gnu.org/licenses/>`_.


