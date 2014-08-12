.. _pairedbamtobed1212:

###############
*pairedbamtobed12*
###############
``bedtools pairedbamtobed12`` Converts 'properly paired' BAM alignments to BED12 format.
Typically producing a 2 blocks BED12 entry for each 'properly paired' BAM pair
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
by default it processes a ‘properly paired’ pair of reads 
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
  
