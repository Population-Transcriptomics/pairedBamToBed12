/*****************************************************************************
  pairedBamToBed12.cpp

  2014      - Nicolas Bertin <nbertin@gsc.riken.jp> and
            - Mickaël Mendez <mickael.mendez@riken.jp>
  2015-2016 - Charles Plessy <plessy@riken.jp>
  CLST-DGT RIKEN Yokohama

  This work was supported by a research grant from the Japanese Ministry of
  Education, Culture, Sports, Science and Technology (MEXT) for the RIKEN Omics
  Science Center to Yoshihide Hayashizaki, and by a research grant from MEXT to
  the RIKEN Center for Life Science Technologies.

  directly inspired from 
     bamToBed.cpp
     (c) 2009 - Aaron Quinlan
     Hall Laboratory
     Department of Biochemistry and Molecular Genetics
     University of Virginia
     aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "BlockedIntervals.h"
#include "bedFile.h"
using namespace BamTools;

#include <vector>
#include <map>
#include <algorithm>    // for swap()
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

// define our program name
#define PROGRAM_NAME "pairedBamToBed12"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void pairedbamtobed12_help(void);

void ConvertPairedBamToBed12(const string &bamFile, int minMapQuality,
                             bool trackUnprocessed, const string &unprocessedBamFile,
                             bool delAsBlock, bool extraG, const string &color, const string &nameSeparator);

void PrintPairedBed12(const BamAlignment &bam1, const BamAlignment &bam2, 
                      const RefVector &refs, bool delAsBlock, bool extraG,
                      string color = "255,0,0");

void ParseCigarBed12(const vector<CigarOp> &cigar, bool delAsBlock, 
                     int &currStart, vector<int> &blockStarts, 
                     vector<int> &blockEnds, int &alignmentEnd);

void RenameRead(BamAlignment& bam);

void SimpleGCorrection(const BamAlignment& bam, string strand, unsigned int &alignmentStart, unsigned int &alignmentEnd);

string TruncateReadName(string name, string separator);

void SaveRead(const BamAlignment& bam, BamWriter& writer,
              bool trackUnprocessed);

int pairedbamtobed12_main(int argc, char* argv[]) {
    
    // our configuration variables
    bool showHelp = false;
  
    // input files
    string bamFile            = "stdin";
    string unprocessedBamFile = "unprocessedPair.bam";
    string color              = "255,0,0";
    string nameSeparator      = "___";
    
    bool haveBam              = true;
    bool delAsBlock           = false;
    bool extraG               = false;
    bool trackUnprocessed     = false;
    
    int minMapQuality         = 0;

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);
        
        if((PARAMETER_CHECK("-h", 2, parameterLength)) || 
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
        else if(PARAMETER_CHECK("-dblock", 7, parameterLength)) {
            delAsBlock = true;
        }
        else if(PARAMETER_CHECK("-extraG", 7, parameterLength)) {
            extraG = true;
        }
        else if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bamFile = argv[i + 1];
                i++;
          }
        }
        else if(PARAMETER_CHECK("-color", 6, parameterLength)) {
            if ((i+1) < argc) {
                color = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-nsep", 5, parameterLength)) {
            if ((i+1) < argc) {
                nameSeparator = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-qual", 5, parameterLength)) {
            if ((i+1) < argc) {
                minMapQuality = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-x", 2, parameterLength)) {
            trackUnprocessed = true;
            if ((i+1) < argc) {
                unprocessedBamFile = argv[i + 1];
                i++;
            }
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }
  
    // make sure we have an input files
    if (haveBam == false) {
      cerr << endl << "*****" << endl << "*****ERROR: Need -i (BAM) file. " << endl << "*****" << endl;
      showHelp = true;
    }

    // if there are no problems, let's convert BAM to BED12
    if (!showHelp) {
        ConvertPairedBamToBed12(bamFile, minMapQuality, 
                                trackUnprocessed, unprocessedBamFile, 
                                delAsBlock, extraG, color, nameSeparator);
    }
    else {
        pairedbamtobed12_help();
    }
    return 0;
}

void pairedbamtobed12_help(void) {

    cerr << "\nTool:    pairedBamToBed12 1.1" << endl;
    cerr << "Summary: Converts 'properly paired' BAM alignments to BED12 format." << endl;
    cerr << "         Typically producing a 2 blocks BED12 entry for each 'properly paired' BAM pair." << endl;
    cerr << "         Additional blocks are produced when an alignment contains long deletion (CIGAR N-op)." << endl;
    cerr << "         The BAM input file must be grouped/sorted by query name (not alignment position)." << endl << endl;
    
    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bam> " << endl << endl;
    
    cerr << "Options: " << endl;
        
    cerr << "\t-dblock\t"  << "Triggers the creation of a new block when an alignment contains short deletion from reference (CIGAR D-op)" << endl << endl;
    
    cerr << "\t-extraG\t"  << "Ignore G mismatches on first bases (for use on CAGE alignments with BWA aln)." << endl << endl;
    
    cerr << "\t-color\t"   << "An R,G,B string for the color used with BED12 format." << endl;
    cerr                   << "\t\tDefault is (255,0,0)." << endl << endl;

    cerr << "\t-nsep\t"    << "A string after which the read names are allowed to differ." << endl;
    cerr                   << "\t\tDefault is ___.  Give an improbable value like 'nothankyou' to turn off." << endl << endl;
    
    cerr << "\t-qual\t"    << "The minimum (inclusive) mapQ sum for reporting the paired BAM into a BED12." << endl;
    cerr                   << "\t\tDefault is (0)." << endl << endl;
    
    cerr << "\t-x\t"       << "Optional filename where unprocessed mapped pairs can be stored." << endl << endl;
    
    // end the program here
    exit(1);
}

 
/*
  Assumptions:
     1.  The BAM file is grouped/sorted by query name,
         not alignment position
     2.  The BAM file only contains properly paired reads
*/
void ConvertPairedBamToBed12(const string &bamFile, int minMapQuality,
                             bool trackUnprocessed, const string &unprocessedBamFile,
                             bool delAsBlock, bool extraG, const string &color, const string &nameSeparator) 
{
    bool lastReadHasMate;
    // open the BAM file
    BamReader reader;
    if (!reader.Open(bamFile)) {
        cerr << "Failed to open BAM file " << bamFile << endl;
        exit(1);
    }

    // get header & reference information
    string header = reader.GetHeaderText();
    RefVector refs = reader.GetReferenceData();

    
    BamWriter writer;
    if (trackUnprocessed == true){
        if (!writer.Open(unprocessedBamFile, header, refs)){
            cerr << "Failed to write into BAM file " << unprocessedBamFile << endl;
        }
    }
    BamAlignment bam1, bam2;
    while (reader.GetNextAlignment(bam1)) { //while1
        RenameRead(bam1);
        lastReadHasMate = false;
        
        while (reader.GetNextAlignment(bam2)) { //while2
            lastReadHasMate = true;
            RenameRead(bam2);
            
            // check if reads are properly paired
            if ((!bam1.IsProperPair()) || (!bam2.IsProperPair())){
                cerr << "*****WARNING: Query " << bam1.Name 
                     << " is not a proper pair. Skipping" << endl;
                SaveRead(bam1, writer, trackUnprocessed);
                lastReadHasMate = false;
                bam1 = bam2;
                continue;
            }
            // check if reads name are similar
            else if (TruncateReadName(bam1.Name, nameSeparator) != TruncateReadName(bam2.Name, nameSeparator)) {
                cerr << "*****WARNING: Query " << bam1.Name 
                     << " is not followed by his mate in your BAM file. Skipping" << endl;
                SaveRead(bam1, writer, trackUnprocessed);
                lastReadHasMate = false;
                bam1 = bam2;
                continue;
            } 
            // check if the reads are on the same chrom
            else if (refs.at(bam1.RefID).RefName != refs.at(bam2.RefID).RefName){
                cerr << "*****WARNING: Query " << bam1.Name 
                     << " is not on the same chromosome than his mate. Skipping" << endl;                
                SaveRead(bam1, writer, trackUnprocessed);
                lastReadHasMate = false;
                bam1 = bam2;
                continue;
            }
            // check the mapping quality
            else if (bam1.MapQuality + bam2.MapQuality < minMapQuality) {
                SaveRead(bam1, writer, trackUnprocessed);
                SaveRead(bam2, writer, trackUnprocessed);
                break;
            }
            
            // make sure that read1 is before read2
            if (bam1.Position > bam2.Position) {
                swap(bam1, bam2);
            }
            
            PrintPairedBed12(bam1, bam2, refs, delAsBlock, extraG, color);
            break;
        } //end while2

        // lastReadHasMate = false when the last read has no mate
        if (!lastReadHasMate) {
            cerr << "*****WARNING: Query " << bam1.Name
                 << " is the last read and has no mate. Skip and exit. " << endl;
            SaveRead(bam1, writer, trackUnprocessed);
        } 
    } //end while1
    
    reader.Close();

    if (trackUnprocessed){
        writer.Close();
    } 

    
}

void ParseCigarBed12(const vector<CigarOp> &cigar, bool delAsBlock, unsigned int &currStart, vector<int> &blockStarts, vector<int> &blockLengths, unsigned int &alignmentEnd) {

    int blockLength = 0;
    //  Rip through the CIGAR ops and figure out if there is more
    //  than one block for this alignment
    vector<CigarOp>::const_iterator cigItr = cigar.begin();
    vector<CigarOp>::const_iterator cigEnd = cigar.end();
    for (; cigItr != cigEnd; ++cigItr) {
        switch (cigItr->Type) {
            case ('M'):
                blockLength += cigItr->Length;
                currStart += cigItr->Length;
            case ('I'): break;
            case ('S'): break;
            case ('D'):
                if (delAsBlock == true) {
                    blockStarts.push_back(currStart + cigItr->Length);
                    blockLengths.push_back(blockLength);
                    currStart += cigItr->Length;
                    blockLength = 0;
                } else {
                    blockLength += cigItr->Length;
                    currStart += cigItr->Length;
                }
                break;
            case ('P'): break;
            case ('N'):
                blockStarts.push_back(currStart + cigItr->Length);
                blockLengths.push_back(blockLength);
                currStart += cigItr->Length;
                blockLength = 0;
            case ('H'): break; // for 'H' - do nothing, move to next op
            default:
                printf("%s *****ERROR: Invalid Cigar op type\n", PROGRAM_NAME); // shouldn't get here
                exit(1);
        }
    }
    // add the last block and set the alignment end 
    blockLengths.push_back(blockLength);
    alignmentEnd = currStart;
}

void SimpleGCorrection(const BamAlignment& bam, string strand, unsigned int &alignmentStart, unsigned int &alignmentEnd) {
    string md;
    bam.GetTag("MD", md);
    if ( (strand == "+") & (bam.QueryBases.substr(0,1) == "G") )  {
      md = md.substr(0,2);
      if (md == "0A" || md == "0C" || md == "0T") {
        alignmentStart = alignmentStart +1;
      }
    }
    if ((strand == "-") & (bam.QueryBases.substr(bam.QueryBases.length() -1, 1) == "C") ) {
      md = md.substr(md.length() -2, 2);
      if (md == "A0" || md == "G0" || md == "T0") {
        alignmentEnd = alignmentEnd -1;
      }
    }
}

void PrintPairedBed12(const BamAlignment &bam1, const BamAlignment &bam2, const RefVector &refs, bool delAsBlock, bool extraG, string color) {

    // set the chrom
    string chrom = refs.at(bam1.RefID).RefName;

    // Get the name and strand of the first mate.
    // Since bam1 and bam2 are proper pairs, and since they are sorted by
    // position (see the `swap` comand above), we are sure that bam1 is on
    // the plus strand and bam2 on the minus strand.

    string strand = "+";
    string name = "";
    if (bam1.IsFirstMate()) {
        name = bam1.Name;
    } else {
        strand = "-";
        name = bam2.Name;
    }

    // parse the CIGAR string and figure out the alignment blocks
    unsigned int bam1_alignmentEnd;
    unsigned int bam2_alignmentEnd;
    vector<int> blockLengths;
    vector<int> blockStarts;
    // very first block
    unsigned int currPosition1 = 0;
    blockStarts.push_back(currPosition1);
    // extract the additional block starts and lengths from the CIGAR string of bam1
    ParseCigarBed12(bam1.CigarData, delAsBlock, currPosition1, blockStarts, blockLengths, bam1_alignmentEnd);
    // new block with bam2.position start
    unsigned int currPosition2 = bam2.Position - bam1.Position;
    blockStarts.push_back(currPosition2);
    // extract the additional block starts and lengths from the CIGAR string of bam2
    ParseCigarBed12(bam2.CigarData, delAsBlock, currPosition2, blockStarts, blockLengths, bam2_alignmentEnd);

    // set the start and and position of the BED12
    unsigned int alignmentStart;
    unsigned int alignmentEnd;
    alignmentStart = bam1.Position;
    alignmentEnd = alignmentStart + bam2_alignmentEnd;

    // Shift TSS left or right if first base looks like G-addition.
    if (extraG) SimpleGCorrection(bam1, strand, alignmentStart, alignmentEnd);

    // write BED6 portion
    // the score is the sum of the MapQ
    printf("%s\t%d\t%d\t%s\t%d\t%s\t",
            chrom.c_str(),
            alignmentStart,
            alignmentEnd,
            name.c_str(),
            bam1.MapQuality + bam2.MapQuality,
            strand.c_str()
            );
    // write the remaining BED12 fields
    // write the txStart / txEnd mark the extend of the 5'read block(s)
    if (bam1.IsFirstMate()) {
        printf("%d\t%d\t", bam1.Position, bam1.Position + bam1_alignmentEnd);
    } else {
        printf("%d\t%d\t", bam2.Position, alignmentEnd);
    }
    printf("%s\t%d\t",
            color.c_str(),
            (int) blockStarts.size()
            );

    // write the comma delimited blockSizes
    unsigned int b;
    for (b = 0; b < blockLengths.size() - 1; ++b) {
        printf("%d,", blockLengths[b]);
    }
    printf("%d\t", blockLengths[b]);

    // write the comma delimited blockStarts
    for (b = 0; b < blockStarts.size() - 1; ++b) {
        printf("%d,", blockStarts[b]);
    }
    printf("%d\n", blockStarts[b]);
}

/* 
 * Remove '/' and what follows when '/' is found at the end of the read's name
 */
void RenameRead(BamAlignment& bam) {
    if (bam.Name.length() < 2)
        return;

    if (bam.Name.at(bam.Name.length() - 2) == '/')
        bam.Name = bam.Name.substr(0, bam.Name.length() - 2);
}
    
// Accessory function to make name comparison more readable.
string TruncateReadName(string name, string separator) {
    return name.substr(0, name.find(separator));
}

void SaveRead(const BamAlignment& bam, BamWriter& writer,
              bool trackUnprocessed){
    
    if (trackUnprocessed){
        writer.SaveAlignment(bam);
    }
}
