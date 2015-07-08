/*****************************************************************************
  bedtools.cpp

  Modified in 2015 by Charles Plessy <plessy@riken.jp> to run pairedBamToBed12
  unconditionally.
  
  bedtools is:
  (c) 2009-2011 - Aaron Quinlan
  Quinlan Laboratory
  Department of Public Health Sciences
  Center for Public Health genomics
  University of Virginia
  aaronquinlan@gmail.com
  
  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "pairedBamToBed12"

// colors for the term's menu 
#define RESET "\033[m"
#define GREEN "\033[1;32m"
#define BLUE "\033[1;34m"
#define RED "\033[1;31m"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

int pairedbamtobed12_main(int argc, char* argv[]); //

int main(int argc, char *argv[])
{
    return pairedbamtobed12_main(argc, argv);
}
