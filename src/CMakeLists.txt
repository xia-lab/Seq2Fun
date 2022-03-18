cmake_minimum_required(VERSION 3.16)

set(ZLIB_NAME "zlib")
set(NCBI_BLAST_INC "include/ncbi-blast+")

#set(CMAKE_EXE_LINKER_FLAGS "-ldl -lpthread -lz")
set(CMAKE_C_COMPILER /usr/bin/gcc)
set(CMAKE_CXX_COMPILER /usr/bin/g++)

set(CMAKE_C_FLAGS "-O3 -g -Wno-unused-result -DNDEBUG")
set(CMAKE_CXX_FLAGS "-O3 -pthread -Wno-unused-result -std=c++11 -DNDEBUG")

#set(BWT_OBJS "bwt/bwt.o bwt/compactfmi.o bwt/sequence.o bwt/suffixArray.o")

#set(CMAKE_EXE_LINKER_FLAGS ${BWT_OBJS})

#project(Seq2Fun)

set(BWT_EXES mkbwt mkfmi)

add_custom_command(OUTPUT ${BWT_EXES}
  COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/bwt/
  COMMAND echo "running custom command"
)

#add_custom_target(run ALL
#  DEPENDS ${BWT_EXES}
#)


include_directories(
  include
  ${NCBI_BLAST_INC}
  ${ZLIB_NAME}
)

set(BLASTSRCS
  ./include/ncbi-blast+/algo/blast/core/pattern.c
  ./include/ncbi-blast+/algo/blast/core/blast_posit.c
  ./include/ncbi-blast+/algo/blast/composition_adjustment/matrix_frequency_data.c
  ./include/ncbi-blast+/algo/blast/core/blast_dynarray.c
  ./include/ncbi-blast+/algo/blast/core/matrix_freq_ratios.c
  ./include/ncbi-blast+/algo/blast/core/blast_encoding.c
  ./include/ncbi-blast+/algo/blast/core/blast_stat.c
  ./include/ncbi-blast+/algo/blast/core/blast_filter.c
  ./include/ncbi-blast+/algo/blast/core/blast_util.c
  ./include/ncbi-blast+/algo/blast/core/blast_message.c
  ./include/ncbi-blast+/algo/blast/core/ncbi_erf.c
  ./include/ncbi-blast+/algo/blast/core/blast_options.c
  ./include/ncbi-blast+/algo/blast/core/ncbi_math.c
  ./include/ncbi-blast+/algo/blast/core/blast_program.c
  ./include/ncbi-blast+/algo/blast/core/ncbi_std.c
  ./include/ncbi-blast+/algo/blast/core/blast_psi_priv.c
  ./include/ncbi-blast+/util/tables/raw_scoremat.c
  ./include/ncbi-blast+/algo/blast/core/blast_query_info.c
  ./include/ncbi-blast+/algo/blast/core/blast_seg.c
)

set(BWTSRCS
#  ./bwt/multikeyqsort.c
  ./bwt/bwt.c
#  ./bwt/suffixArray.c
  ./bwt/compactfmi.c
#  ./bwt/readFasta.c
  ./bwt/sequence.c
#  ./bwt/mkfmi.c
#  ./bwt/mkbwt.c
  ./bwt/suffixArray.c
)

set(SEQ2FUN_SRCS
  seq2fun.cpp
  transsearcher.cpp
  fragment.cpp
  bwtfmiDB.cpp
  adaptertrimmer.cpp
  basecorrector.cpp
  duplicate.cpp
  evaluator.cpp
  fastareader.cpp
  fastqreader.cpp
  filter.cpp
  filterresult.cpp
  htmlreporterall.cpp
  htmlreporter.cpp
  jsonreporter.cpp
  nucleotidetree.cpp
  options.cpp
  overlapanalysis.cpp
  peprocessor.cpp
  polyx.cpp
  processor.cpp
  read.cpp
  seprocessor.cpp
  sequence.cpp
  stats.cpp
  threadconfig.cpp
  umiprocessor.cpp
  unittest.cpp
  writer.cpp
  writerthread.cpp
)

set(SEQTRACT_SRCS
  seqtract.cpp
  transsearcher.cpp
  fragment.cpp
  bwtfmiDB.cpp
  adaptertrimmer.cpp
  basecorrector.cpp
  duplicate.cpp
  evaluator.cpp
  fastareader.cpp
  fastqreader.cpp
  filter.cpp
  filterresult.cpp
  htmlreporterall.cpp
  htmlreporter.cpp
  jsonreporter.cpp
  nucleotidetree.cpp
  options.cpp
  overlapanalysis.cpp
  peprocessor.cpp
  polyx.cpp
  processor.cpp
  read.cpp
  seprocessor.cpp
  sequence.cpp
  stats.cpp
  threadconfig.cpp
  umiprocessor.cpp
  unittest.cpp
  writer.cpp
  writerthread.cpp
  seqtractpeprocessor.cpp
  threadsconfig2.cpp
)
#message(STATUS "BLASTOBJS is: ${BLASTOBJS}")
message(STATUS "in src CMAKE_INSTALL_PREFIX is: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "in src EXECUTABLE_OUTPUT_PATH is ${EXECUTABLE_OUTPUT_PATH}")
message(STATUS "in src CMAKE_CURRENT_SOURCE_DIR is ${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(seq2fun
  ${BWTSRCS}
  ${SEQ2FUN_SRCS}
  ${BLASTSRCS}
  #${BWT_EXES}
)

#add_dependencies(seq2fun
#  ${BWT_OBJS}
#)

add_executable(seqtract
  ${BWTSRCS}
  ${SEQTRACT_SRCS}
  ${BLASTSRCS}
  #DEPENDS ${BWT_EXES}
)

target_link_libraries(seq2fun
  -lpthread 
  -lz
  -lm
)

target_link_libraries(seqtract
  -lpthread 
  -lz
  -lm
)

#add_subdirectory(include)

#ExternalProject_Add(${ZLIB_NAME} 
#  PREFIX zlib
#)

