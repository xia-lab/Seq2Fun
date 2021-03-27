CC = gcc
CXX = g++
CFLAGS = -O3 -DNDEBUG
CXXFLAGS = -O3 -pthread -std=c++11 -DNDEBUG
#LDLIBS = -ldl -lpthread -lz
LDLIBS = -lpthread -lz
INCLUDES	= -I./include -I./include/ncbi-blast+
#INCLUDES	= -I./include -I./include/ncbi-blast+ -I./zlib

BLASTOBJS = include/ncbi-blast+/algo/blast/core/pattern.o \
			 include/ncbi-blast+/algo/blast/core/blast_posit.o \
			 include/ncbi-blast+/algo/blast/composition_adjustment/matrix_frequency_data.o \
			 include/ncbi-blast+/algo/blast/core/blast_dynarray.o \
			 include/ncbi-blast+/algo/blast/core/matrix_freq_ratios.o \
			 include/ncbi-blast+/algo/blast/core/blast_encoding.o \
			 include/ncbi-blast+/algo/blast/core/blast_stat.o \
			 include/ncbi-blast+/algo/blast/core/blast_filter.o \
			 include/ncbi-blast+/algo/blast/core/blast_util.o \
			 include/ncbi-blast+/algo/blast/core/blast_message.o \
			 include/ncbi-blast+/algo/blast/core/ncbi_erf.o \
			 include/ncbi-blast+/algo/blast/core/blast_options.o \
			 include/ncbi-blast+/algo/blast/core/ncbi_math.o \
			 include/ncbi-blast+/algo/blast/core/blast_program.o \
			 include/ncbi-blast+/algo/blast/core/ncbi_std.o \
			 include/ncbi-blast+/algo/blast/core/blast_psi_priv.o \
			 include/ncbi-blast+/util/tables/raw_scoremat.o \
			 include/ncbi-blast+/algo/blast/core/blast_query_info.o \
			 include/ncbi-blast+/algo/blast/core/blast_seg.o 


BWTOBJS = bwt/bwt.o bwt/compactfmi.o bwt/sequence.o bwt/suffixArray.o

ifeq ($(uname -s), "Darwin")
LD_LIBS_STATIC = -Wl,-all_load -lpthread -lz -Wl,-noall_load
else
LD_LIBS_STATIC = -Wl,--whole-archive -lpthread -lz -Wl,--no-whole-archive
endif

all: makefile  seq2fun bwt/mkbwt
	mkdir -p ../bin
	cp seq2fun ../bin/
	cp bwt/mkbwt ../bin/mkbwt
	cp bwt/mkfmi ../bin/mkfmi

# use bwt/mkbwt as target for compiling everything in the bwt folder first
bwt/mkbwt:
	$(MAKE) -C bwt/ $(MAKECMDGOALS)

seq2fun: makefile bwt/mkbwt seq2fun.o transsearcher.o fragment.o bwtfmiDB.o \
	adaptertrimmer.o basecorrector.o \
	duplicate.o evaluator.o fastareader.o fastqreader.o filter.o filterresult.o htmlreporter.o htmlreporterall.o \
	jsonreporter.o  nucleotidetree.o options.o overlapanalysis.o peprocessor.o \
	polyx.o processor.o read.o seprocessor.o sequence.o stats.o threadconfig.o umiprocessor.o \
	unittest.o writer.o writerthread.o $(BLASTOBJS)
	$(CXX) $(LDFLAGS) -o seq2fun seq2fun.o transsearcher.o fragment.o bwtfmiDB.o \
	adaptertrimmer.o basecorrector.o \
	duplicate.o evaluator.o fastareader.o fastqreader.o filter.o filterresult.o htmlreporter.o htmlreporterall.o \
	jsonreporter.o  nucleotidetree.o options.o overlapanalysis.o peprocessor.o \
	polyx.o processor.o read.o seprocessor.o sequence.o stats.o threadconfig.o umiprocessor.o \
	unittest.o writer.o writerthread.o $(BWTOBJS) $(BLASTOBJS) $(LDLIBS)

#%.o : %.c makefile
%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
#%.o : %.cpp makefile
%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -f -v bwt/mkbwt bwt/mkfmi  seq2fun ../bin/*
	find . -name "*.o" -delete
	$(MAKE) -C bwt/ clean

static: LDFLAGS = -static
static: LDLIBS = $(LD_LIBS_STATIC)
static: all

debug: CXXFLAGS = -O3 -pthread -std=c++11 -g -Wall -Wpedantic -Wextra -Wconversion -fno-omit-frame-pointer
debug: CFLAGS = -g -O3 -Wall -Wno-uninitialized
debug: all

.PHONY: clean debug static
