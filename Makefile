
# ==========================
# BEDTools Makefile
# (c) 2009 Aaron Quinlan
#
# Modified in 2015 by Charles Plessy <plessy@riken.jp> to
# only build pairedBamToBed12
#
# ==========================

SHELL := /bin/bash -e

# define our object and binary directories
export OBJ_DIR	= obj
export BIN_DIR	= bin
export SRC_DIR	= src
export UTIL_DIR	= src/utils
export CXX		= g++
#ifeq ($(DEBUG),1)
#export CXXFLAGS = -Wall -O0 -g -fno-inline -fkeep-inline-functions -D_FILE_OFFSET_BITS=64 -fPIC -DDEBUG -D_DEBUG
#else
export CXXFLAGS = -Wall -O2 -D_FILE_OFFSET_BITS=64 -fPIC $(INCLUDE)
#endif
export LIBS		= -lz
export BT_ROOT  = src/utils/BamTools/

SUBDIRS = $(SRC_DIR)/pairedBamToBed12

UTIL_SUBDIRS =	$(SRC_DIR)/utils/bedFile \
				$(SRC_DIR)/utils/BinTree \
				$(SRC_DIR)/utils/bedGraphFile \
				$(SRC_DIR)/utils/chromsweep \
				$(SRC_DIR)/utils/Contexts \
				$(SRC_DIR)/utils/FileRecordTools \
				$(SRC_DIR)/utils/general \
				$(SRC_DIR)/utils/gzstream \
				$(SRC_DIR)/utils/fileType \
				$(SRC_DIR)/utils/bedFilePE \
				$(SRC_DIR)/utils/KeyListOps \
				$(SRC_DIR)/utils/NewChromsweep \
				$(SRC_DIR)/utils/sequenceUtilities \
				$(SRC_DIR)/utils/tabFile \
				$(SRC_DIR)/utils/BamTools \
				$(SRC_DIR)/utils/BamTools-Ancillary \
				$(SRC_DIR)/utils/BlockedIntervals \
				$(SRC_DIR)/utils/Fasta \
				$(SRC_DIR)/utils/VectorOps \
				$(SRC_DIR)/utils/GenomeFile \
				$(SRC_DIR)/utils/RecordOutputMgr

BUILT_OBJECTS = $(OBJ_DIR)/*.o

all: $(OBJ_DIR) $(BIN_DIR) $(UTIL_SUBDIRS) $(SUBDIRS)
	@echo "- Building main pairedBamToBed12 binary."
	@$(CXX) $(CXXFLAGS) -c src/bedtools.cpp -o obj/bedtools.o -I$(UTIL_DIR)/version/
	@$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/pairedBamToBed12 $(BUILT_OBJECTS) -L$(UTIL_DIR)/BamTools/lib/ -lbamtools $(LIBS) $(LDFLAGS)
	@echo "done."

.PHONY: all

# make the "obj/" and "bin/" directories, if they don't exist
$(OBJ_DIR) $(BIN_DIR):
	@mkdir -p $@

# One special case: All (or almost all) programs requires the BamTools API files to be created first.
.PHONY: bamtools_api
bamtools_api:
	@$(MAKE) --no-print-directory --directory=$(BT_ROOT) api
$(UTIL_SUBDIRS) $(SUBDIRS): bamtools_api

# even though these are real directories, treat them as phony targets, forcing to always go in them are re-make.
# a future improvement would be the check for the compiled object, and rebuild only if the source code is newer.
.PHONY: $(UTIL_SUBDIRS) $(SUBDIRS)
$(UTIL_SUBDIRS) $(SUBDIRS): $(OBJ_DIR) $(BIN_DIR)
	@echo "- Building in $@"
	@$(MAKE) --no-print-directory --directory=$@

clean:
	@$(MAKE) --no-print-directory --directory=$(BT_ROOT) clean_api
	@echo " * Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
.PHONY: clean

test: all
	@cd test; bash test.sh

.PHONY: test
