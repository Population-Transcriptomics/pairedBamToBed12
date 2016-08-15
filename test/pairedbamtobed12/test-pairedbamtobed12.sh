BT=${BT-../../bin/pairedBamToBed12}

SUCCESS=0
FAILURES=0
check()
{
	if diff $1 $2; then
    		echo ok
		SUCCESS=$((SUCCESS + 1))
	else
    		echo fail
		FAILURES=$((FAILURES +1))
	fi
}

###########################################################
###########################################################
#                       BAM files                         #
###########################################################
###########################################################
samtools view -Sb pair-sorted.sam > pair-sorted.bam 2>/dev/null
samtools view -Sb pair-unsorted.sam > pair-unsorted.bam 2>/dev/null
samtools view -Sb pair-and-missingmate.sam > pair-and-missingmate.bam 2>/dev/null
samtools view -Sb missingmate-and-pair.sam > missingmate-and-pair.bam 2>/dev/null
samtools view -Sb pair-with-name-containing-slash.sam > pair-with-name-containing-slash.bam 2>/dev/null
samtools view -Sb proper_pair_plus_strand.sam > proper_pair_plus_strand.bam 2> /dev/null
samtools view -Sb proper_pair_minus_strand.sam > proper_pair_minus_strand.bam 2> /dev/null
samtools view -Sb proper_pair_spliced.sam > proper_pair_spliced.bam 2> /dev/null
samtools view -Sb proper_pair_overlap.sam > proper_pair_overlap.bam 2> /dev/null
samtools view -Sb proper_pair_bad_mapq.sam > proper_pair_bad_mapq.bam 2> /dev/null
samtools view -Sb not_proper_pair.sam > not_proper_pair.bam 2> /dev/null
samtools view -Sb bug_proper_pair_different_chrom.sam > bug_proper_pair_different_chrom.bam 2> /dev/null
samtools view -Sb proper_pair_nsep-option.sam > proper_pair_nsep-option.bam 2> /dev/null


##################################################################
#  Test paired reads when sorted by name and position
##################################################################
printf "Test  1: paired reads when sorted by name and position...\n"
echo \
"chr14	50053297	50053480	M00528:19:000000000-A88YD:1:1101:2241:12366	0	+	50053297	50053324	255,0,0	2	27,21	0,162" > exp
$BT -i pair-sorted.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test paired reads when sorted by name and not by position
##################################################################
echo "Test  2: paired reads when sorted by name and not by position..."
echo \
"chr14	50053297	50053480	M00528:19:000000000-A88YD:1:1101:2241:12366	0	+	50053297	50053324	255,0,0	2	27,21	0,162" > exp
$BT -i pair-unsorted.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test skipping missing mate first read
##################################################################
echo "Test  3: skipping missing mate first read..."
echo \
"chr4	8210431	8210761	M00528:19:000000000-A88YD:1:1101:2318:12845	120	+	8210431	8210458	255,0,0	2	27,21	0,309" > exp

$BT -i missingmate-and-pair.bam > obs 2> /dev/null
check obs exp
rm obs exp


##################################################################
#  Test skipping missing mate first read (test stderr)
##################################################################
echo "Test  4: skipping missing mate first read (test stderr)..."
echo \
"*****WARNING: Query M00528:19:000000000-A88YD:1:1101:2241:12366 is not followed by his mate in your BAM file. Skipping" > exp
$BT -i missingmate-and-pair.bam > /dev/null 2> obs
check obs exp
rm obs exp


##################################################################
#  Test skipping missing mate last read
##################################################################
echo "Test  5: skipping missing mate last read..."
echo \
"chr14	50053297	50053480	M00528:19:000000000-A88YD:1:1101:2241:12366	0	+	50053297	50053324	255,0,0	2	27,21	0,162" > exp
$BT -i pair-and-missingmate.bam > obs 2> /dev/null
check obs exp
rm obs exp


##################################################################
#  Test skipping missing mate last read (stderr)
##################################################################
echo "Test  6: skipping missing mate last read (stderr)..."
echo \
"*****WARNING: Query M00528:19:000000000-A88YD:1:1101:2318:12845 is the last read and has no mate. Skip and exit. " > exp
$BT -i pair-and-missingmate.bam > /dev/null 2> obs
check obs exp
rm obs exp
#exit


##################################################################
#  Test paired reads with name ending with '/*'
##################################################################
echo "Test  7: paired reads with name ending with '/*'..."
echo \
"chr14	50053297	50053480	M00528:19:000000000-A88YD:1:1101:2241:12366	0	+	50053297	50053324	255,0,0	2	27,21	0,162" > exp
$BT -i pair-with-name-containing-slash.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test proper pairs on plus or minus strands
##################################################################
echo "Test 8a: proper pair on the plus strand..."
echo \
"chr1	0	99	proper_pair_plus_strand_no_mismatch	80	+	0	30	255,0,0	2	30,30	0,69
chr2	0	99	proper_pair_plus_strand_with_extra_G	80	+	0	30	255,0,0	2	30,30	0,69
chr3	0	99	proper_pair_plus_strand_other_mismatch1	80	+	0	30	255,0,0	2	30,30	0,69
chr3	0	99	proper_pair_plus_strand_other_mismatch2	80	+	0	30	255,0,0	2	30,30	0,69" > exp
$BT -i proper_pair_plus_strand.bam > obs
check obs exp
rm obs exp

echo "Test 8b: proper pair on the minus strand..."
echo \
"chr1	0	99	proper_pair_minus_strand_no_mismatch	80	-	69	99	255,0,0	2	30,30	0,69
chr2	0	99	proper_pair_minus_strand_with_extra_G	80	-	69	99	255,0,0	2	30,30	0,69
chr3	0	99	proper_pair_minus_strand_with_extra_C	80	-	69	99	255,0,0	2	30,30	0,69
chr4	0	99	proper_pair_minus_strand_other_mismatch	80	-	69	99	255,0,0	2	30,30	0,69" > exp
$BT -i proper_pair_minus_strand.bam > obs
check obs exp
rm obs exp

##################################################################
# Test proper pairs with removal of extra Gs
##################################################################

echo "Test 9a: proper pair on the plus strand, with extraG handling..."
echo \
"chr1	0	99	proper_pair_plus_strand_no_mismatch	80	+	0	30	255,0,0	2	30,30	0,69
chr2	1	99	proper_pair_plus_strand_with_extra_G	80	+	1	30	255,0,0	2	29,30	0,68
chr3	0	99	proper_pair_plus_strand_other_mismatch1	80	+	0	30	255,0,0	2	30,30	0,69
chr3	0	99	proper_pair_plus_strand_other_mismatch2	80	+	0	30	255,0,0	2	30,30	0,69" > exp
$BT -extraG -i proper_pair_plus_strand.bam > obs
check obs exp
rm obs exp

echo "Test 9b: proper pair on the minus strand, with extraG handling..."
echo \
"chr1	0	99	proper_pair_minus_strand_no_mismatch	80	-	69	99	255,0,0	2	30,30	0,69
chr2	0	98	proper_pair_minus_strand_with_extra_G	80	-	69	98	255,0,0	2	30,29	0,69
chr3	0	99	proper_pair_minus_strand_with_extra_C	80	-	69	99	255,0,0	2	30,30	0,69
chr4	0	99	proper_pair_minus_strand_other_mismatch	80	-	69	99	255,0,0	2	30,30	0,69" > exp
$BT -extraG -i proper_pair_minus_strand.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test proper pair spliced
##################################################################
echo "Test 10: proper pair spliced..."
echo \
"chr1	0	99	proper_pair_spliced	80	+	0	40	255,0,0	3	20,10,30	0,30,69" > exp
$BT -i proper_pair_spliced.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test proper pair overlap
##################################################################
echo "Test 11: proper pair overlap..."
echo \
"chr1	0	49	proper_pair_overlap	80	+	0	30	255,0,0	2	30,30	0,19" > exp
$BT -i proper_pair_overlap.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test skipping proper pair with bad mapq when using -qual argument
##################################################################
echo "Test 12: skipping proper pair with bad mapq when using -qual argument..."
touch exp
$BT -i proper_pair_bad_mapq.bam -qual 11 > obs
check obs exp
rm obs exp


##################################################################
#  Test skip not proper pair
##################################################################
echo "Test 13: skip not proper pair..."
touch exp
$BT -i not_proper_pair.bam > obs 2> /dev/null
check obs exp
rm obs exp


##################################################################
#  Test skip not proper pair (stderr)
##################################################################
echo "Test 14: skip not proper pair (stderr)..."
echo \
"*****WARNING: Query not_proper_pair is not a proper pair. Skipping
*****WARNING: Query not_proper_pair is the last read and has no mate. Skip and exit. " > exp
$BT -i not_proper_pair.bam > /dev/null 2> obs 
check obs exp
rm obs exp


##################################################################
#  Test skip proper when reads are on different chromosomes
##################################################################
echo "Test 15: skip proper when reads are on different chromosomes..."
touch exp
$BT -i bug_proper_pair_different_chrom.bam > obs 2> /dev/null
check obs exp
rm obs exp


##################################################################
#  Test skip proper when reads are on different chromosomes (stderr)
##################################################################
echo "Test 16: skip proper when reads are on different chromosomes (stderr)..."
echo \
"*****WARNING: Query bug_proper_pair_different_chrom is not on the same chromosome than his mate. Skipping
*****WARNING: Query bug_proper_pair_different_chrom is the last read and has no mate. Skip and exit. " > exp
$BT -i bug_proper_pair_different_chrom.bam > /dev/null 2> obs
check obs exp
rm obs exp


##################################################################
#  Test -nsep option
##################################################################
echo "Test 17: -nsep option..."
echo \
"chr1	0	99	proper_pair_plus_strand%%%string1	80	+	0	30	255,0,0	2	30,30	0,69" > exp
$BT -nsep %%% -i proper_pair_nsep-option.bam > obs
check obs exp
rm obs exp

echo "Test 18: -nsep option (stderr)..."
echo \
"*****WARNING: Query proper_pair_plus_strand%%%string1 is not followed by his mate in your BAM file. Skipping
*****WARNING: Query proper_pair_plus_strand%%%string2 is the last read and has no mate. Skip and exit. " > exp
$BT -i proper_pair_nsep-option.bam 2> obs
check obs exp
rm obs exp

echo "$SUCCESS success(es) and $FAILURES failure(s)."

rm *.bam

if [ $FAILURES -eq 0 ]
then
	exit 0
else
	exit 1
fi
