.. _pairedbamtobed12:

###############
*pairedbamtobed12*
###############
``bedtools pairedbamtobed12`` converts *properly paired* BAM alignments to BED12 format.
Typically producing a 2 blocks BED12 entry for each *properly paired* BAM pair
Additional blocks are produced when an alignment contains long deletion (CIGAR N-op).
The BAM input file must be grouped/sorted by query name (not alignment position). 

==========================================================================
Usage and option summary
==========================================================================
**Usage**:
::

  bedtools pairedbamtobed12 [OPTIONS] -i <BAM>

**(or)**:
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
**-qual**       The minimum (inclusive) mapQ sum for reporting
                the paired BAM into a BED12. Default is 0.
**-x**          Optional filename where unprocessed mapped pairs can be stored.
=============   ================================================================


==========================================================================
Default behavior
==========================================================================
By default it processes a *properly paired* pair of reads 
(mapped on the same chromosome) into a single BED12 line, where the start
position is the 5′ end of Read 1 and the stop position is the 3′ end of Read 2.
The BED12 blocks are used to indicate positions where the reads match, and the
thick part indicates where is the contribution of Read 1. 

.. note::
    
    When using this tool, it is required that the BAM
    file is sorted/grouped by the read name.

.. note::
    All the reads that are not followed by their mate
    or not properly paired will be automatically skipped by pairedBamToBed12. 
    
    
.. code-block:: bash

  $ bedtools pairedbamtobed12 -i 1proper-pair.bam 
  chr1	50053297	50053480	M00528:19:000000000-A88YD:1:1101:2241:12366	0	+	50053297	50053324	255,0,0	2	27,21	0,162

==========================================================================
Usage with transcriptome libraries
==========================================================================

In transcriptome analysis, the BED12 entries produced by ``pairedbamtobed12``
represent the minimal information about a cDNA that was given by a read pair.

``pairedbamtobed12`` was created for the analysis of CAGEscan_ libraries, which
are we paired-end directional libraries of random-primed 5′ cDNAs.  The BED12
files are used to assemble *CAGEscan clusters* that combine all the pairs where
the 5′ end is in the same transcript start site peak, thus providing approximative
rudimentary transcript models for each peak.  A typical analysis can be found in
`Kratz et al., 2014`_.

This BED12 format is also supported in RIKEN's Zenbu_ genome browser, where one
can load data in this format and visualise it either as genome intervals or as
expression histograms.

.. _CAGEscan:               http://dx.doi.org/10.1038/nmeth.1470
.. _`Kratz et al., 2014`: http://dx.doi.org/10.1101/gr.164095.113
.. _Zenbu:                  http://fantom.gsc.riken.jp/zenbu/

==========================================================================
Advantages and limitations in comparison with ``bamtobed``
==========================================================================

The advantage compared to ``bamtobed -split`` is that ``pairedbamtoned12``
reports the whole pair on a single line, and the advantage compared with
``bamtobed -bedpe``, is that it reports spliced alignments.

The limitation of ``pairedbamtoned12`` is that it only pertains to pairs mapped
on the same chromosome and is therefore unfit for representing gene fusions or
interchromosomal interactions.
