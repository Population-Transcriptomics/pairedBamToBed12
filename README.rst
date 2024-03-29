###############
*pairedBamToBed12*
###############

.. image:: https://travis-ci.org/Population-Transcriptomics/pairedBamToBed12.svg?branch=pairedbamtobed12
    :target: https://travis-ci.org/Population-Transcriptomics/pairedBamToBed12

``pairedBamToBed12`` converts *properly paired* BAM alignments to
BED12 format.  Typical *proper pairs* will be represented by a 2 blocks BED12
entry.  Additional blocks are produced when an alignment contains long deletion
(CIGAR N-op).  Thickness indicates the first read of the pair.  The BAM input
file must be grouped/sorted by query name (not alignment position).

.. code-block::

    Read 1:   >>>>>>>>>>>>
    Read 2:                     <<<<<<<<<<<<<-----<<<<<<<
    The pair: >>>>>>>>>>>>------>>>>>>>>>>>>>----->>>>>>>

==========================================================================
Rationale
==========================================================================

We developped ``pairedBamToBed12`` in order to represent a paired alignment
on a sigle line, with splicing information.  The existing tools did not provide
both features at the same time: ``bedtools bamtobed -split`` represents each
read on separate lines and, ``bedtools bamtobed -bedpe`` does not represent
spliced alignments.

==========================================================================
Installation
==========================================================================

In brief: download the source code, unpack it and enter the main source
directory, type ``make`` and a ``pairedBamToBed12`` executable file will appear
the ``bin`` directory.  Test it with ``make test``.  See below for an example.

.. code-block::

    wget https://github.com/Population-Transcriptomics/pairedBamToBed12/archive/pairedbamtobed12.zip
    unzip pairedbamtobed12.zip
    cd pairedBamToBed12-pairedbamtobed12
    make
    make test
    ls bin/pairedBamToBed12

==========================================================================
Usage and option summary
==========================================================================
**Usage**:
::

    pairedBamToBed12 [OPTIONS] -i <BAM>



.. tabularcolumns:: |p{4.5cm}|p{8.5cm}|

=============   ================================================================
Option          Description
=============   ================================================================
**-dblock**     Triggers the creation of a new block when an alignment contains
                short deletion from reference (CIGAR D-op).
**-color**      An R,G,B string for the color used with BED12 format. Default 
                is (255,0,0).
**-extraG** 	Ignore G mismatches on first bases (experimental option for use
                on CAGE alignments with BWA aln).
**-nsep**       A string after which the read names are allowed to differ.
                Default is ___.  Give an improbable value like 'nothankyou' to turn off.
**-qual**       The minimum (inclusive) mapQ sum for reporting
                the paired BAM into a BED12. Default is 0.
**-x**          Optional filename where unprocessed mapped pairs can be stored.
=============   ================================================================


==========================================================================
Default behavior
==========================================================================
By default it processes a *properly paired* pair of reads into a single BED12
line, where the start and end positions are the 5′ end of Read 1 and the 3′ end
of Read 2.  The BED12 blocks are used to indicate positions where the reads
match, and the thick part indicates where is the contribution of Read 1.  The
relative orientation of the mate pairs must be forward/reverse (which is the
standard in most libraries prepared for the Illumina platform). 

.. note::
    
    The BAM file must be sorted by read name.

.. note::
    Reads that are not followed by their mate or not properly paired will be skipped.

.. code-block:: bash

  $ pairedBamToBed12 -i 1proper-pair.bam 
  chr1	50053297	50053480	M00528:19:000000000-A88YD:1:1101:2241:12366	0	+	50053297	50053324	255,0,0	2	27,21	0,162

==========================================================================
Usage with transcriptome libraries
==========================================================================

In transcriptome analysis, the BED12 entries produced by ``pairedBamToBed12``
represent the minimal information about a cDNA that was given by a read pair.

``pairedBamToBed12`` was created for the analysis of CAGEscan_ libraries, which
are paired-end directional libraries of random-primed 5′ cDNAs.  The BED12
files are used to assemble *CAGEscan clusters* that combine all the pairs where
the 5′ end is in the same transcript start site peak, thus providing approximate
rudimentary transcript models for each peak.  A typical analysis can be found in
`Kratz et al., 2014`_.

This BED12 format is also supported in RIKEN's Zenbu_ genome browser, where one
can load data in this format and visualise it either as genome intervals or as
expression histograms.
    
.. NOTE::
    BWA has a bug that will set the *properly paired* flag for reads where one
    mate is aligned very near the end of a chromosome and the other is aligned
    very near the beginning of the next chromosome, when the ``-a`` option of
    ``sampe`` is large.  However, for CAGEscan, large numbers are necessary to
    span whole gene loci.   It is therefore recommended to sanitise the output
    of BWA with SAMtools, using its ``fixmate`` command, that corrects the
    *properly paired* flag since version 1.0.

.. NOTE::
    CAGE methods sometimes add an extra G at the beginning of the cDNAs (see
    http://population-transcriptomics.org/nanoCAGE/#extra-G).  This leads to
    1-base shifts of some TSS peaks.  From version 1.2, ``pairedBamToBed12``
    provides an experimental option, ``-extraG`` to shift the start or end
    (according to the strand) of the output of one base when a G mismatch
    is detected on the first base of Read1.  A more detailed description of
    the problem may be found in the supplemental material of the FANTOM3_ paper.
    The implementation here is very naive and incomplete, and was tested only on
    data produced by BWA's `sampe` command. It relies on the `MD` flag and does
    not understand clipping.  Thus, the ``-extraG`` option available here is
    not entierly satisfactory and may be removed in the future.  A better
    approach for instance would be to post-process the BAM file instead of
    implementing a correction here.

.. _CAGEscan:               http://dx.doi.org/10.1038/nmeth.1470
.. _`Kratz et al., 2014`: http://dx.doi.org/10.1101/gr.164095.113
.. _Zenbu:                  http://fantom.gsc.riken.jp/zenbu/
.. _FANTOM3:                http://science.sciencemag.org/content/309/5740/1559

==========================================================================
Bugs and limitations
==========================================================================

``pairedbamtobed12`` only pertains to pairs mapped on the same chromosome
and is therefore unfit for representing gene fusions or interchromosomal
interactions.

At the moment, ``pairedbamtobed12`` will output broken BED12 files if Read1
and Read2 overlap and align to the same splice junction.

==========================================================================
Copyright, authorship and license
==========================================================================

``pairedBamToBed12`` is distrubuted under the `GNU General Public License version 2`_.

The tool ``pairedBamToBed12`` is copyright 2013~2015 RIKEN.  It was originally
written by Nicolas Bertin as an addition to `Bedtools`_ 2.11.1.  It was then
ported to Bedtools 2.21.0 by Mickaël Mendez, and then finally forked from the
Bedtools source as a stand-alone program by Charles Plessy.  The documentation
was written by NB, MM and CP, and the regression tests were implemented by MM
and CP.

Bedtools is copyright Aaron Quinlan and others.

.. _`GNU General Public License version 2`: LICENSE
.. _Bedtools: https://github.com/arq5x/bedtools2
