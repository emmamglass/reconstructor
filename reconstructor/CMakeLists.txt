cmake_minimum_required (VERSION 3.15...3.26)
project (libdiamond)

include(CheckCXXCompilerFlag)
include(CheckSymbolExists)
include(CheckTypeSize)

option(BUILD_STATIC "BUILD_STATIC" OFF)
option(EXTRA "EXTRA" OFF)
option(STATIC_LIBGCC "STATIC_LIBGCC" OFF)
option(STATIC_LIBSTDC++ "STATIC_LIBSTDC++" OFF)
option(X86 "X86" ON)
option(ARM "ARM" OFF)
option(AARCH64 "AARCH64" OFF)
option(STRICT_BAND "STRICT_BAND" ON)
option(LEFTMOST_SEED_FILTER "LEFTMOST_SEED_FILTER" ON)
option(SEQ_MASK "SEQ_MASK" ON)
option(DP_STAT "DP_STAT" OFF)
option(SINGLE_THREADED "SINGLE_THREADED" OFF)
option(EIGEN_BLAS "EIGEN_BLAS" OFF)
option(WITH_ZSTD "WITH_ZSTD" OFF)
option(HIT_KEEP_TARGET_ID "HIT_KEEP_TARGET_ID" OFF)
option(LONG_SEEDS "LONG_SEEDS" OFF)
option(WITH_AVX512 "WITH_AVX512" OFF)
option(WITH_DNA "WITH_DNA" OFF)
option(WITH_MCL "WITH_MCL" OFF)
option(WITH_MIMALLOC "WITH_MIMALLOC" OFF)
option(USE_TLS "USE_TLS" OFF)
set(MAX_SHAPE_LEN 19)
set(BLAST_INCLUDE_DIR "" CACHE STRING "BLAST_INCLUDE_DIR")
set(BLAST_LIBRARY_DIR "" CACHE STRING "BLAST_LIBRARY_DIR")

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*|arm64.*)")
  set(X86 OFF)
  set(AARCH64 ON)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm.*|ARM.*)")
  set(X86 OFF)
  set(ARM ON)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "PPC64*|ppc64*|powerpc64*")
  set(X86 OFF)
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -maltivec")
  add_definitions(-DEIGEN_DONT_VECTORIZE)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^s390")
  set(X86 OFF)
endif()

if(STATIC_LIBSTDC++)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")
endif()

if(STRICT_BAND)
  add_definitions(-DSTRICT_BAND)
endif()

if(USE_TLS)
  add_definitions(-DUSE_TLS)
endif()

if(SINGLE_THREADED)
  add_definitions(-DSINGLE_THREADED)
endif()

if(SEQ_MASK)
  add_definitions(-DSEQ_MASK)
endif()

if(LEFTMOST_SEED_FILTER)
  add_definitions(-DLEFTMOST_SEED_FILTER)
endif()

if(DP_STAT)
  add_definitions(-DDP_STAT)
endif()

if(EXTRA)
  add_definitions(-DEXTRA)
endif()

if(HIT_KEEP_TARGET_ID)
  add_definitions(-DHIT_KEEP_TARGET_ID)
endif()

if(LONG_SEEDS)
  add_definitions(-DLONG_SEEDS)
endif()

if(WITH_AVX512)
  add_definitions(-DWITH_AVX512)
endif()

if(WITH_DNA)
  add_definitions(-DWITH_DNA)
endif()

if(WITH_MCL)
  add_definitions(-DWITH_MCL)
endif()

if(WITH_FAMSA)
  add_definitions(-DWITH_FAMSA)
endif()

add_definitions(-DMAX_SHAPE_LEN=${MAX_SHAPE_LEN})
add_definitions(-D_ITERATOR_DEBUG_LEVEL=0)

IF(STATIC_LIBGCC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc")
endif()

if(BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBRARIES OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

function(set_cxx_standard std flag)
  if(${CMAKE_VERSION} VERSION_LESS "3.1.0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
  else()
    set(CMAKE_CXX_STANDARD ${std} PARENT_SCOPE)
  endif()
endfunction(set_cxx_standard)

check_type_size(ptrdiff_t SIZEOF_PTRDIFF_T)
check_type_size(int SIZEOF_INT)

check_cxx_compiler_flag("-std=gnu++14" HAS_GNUPP14)
check_cxx_compiler_flag("-std=gnu++17" HAS_GNUPP17)
if(HAS_GNUPP17 OR ${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  set_cxx_standard(17 "-std=gnu++17")
elseif (HAS_GNUPP14)
  set_cxx_standard(14 "-std=gnu++14")
else()
  set_cxx_standard(11 "-std=gnu++11")
endif()

if(CMAKE_BUILD_MARCH)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${CMAKE_BUILD_MARCH}")
endif()

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_HAS_STD_BYTE=0)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-implicit-fallthrough -Wreturn-type -Wno-unused -Wno-unused-parameter -Wno-unused-variable -Wno-uninitialized -Wno-deprecated-copy -Wno-unknown-warning-option ")#-g -fsanitize=address -fno-omit-frame-pointer ")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-pragma-clang-attribute -Wno-overloaded-virtual -Wno-missing-braces") #-g -fsanitize=thread -fno-omit-frame-pointer" )
endif()

if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fabi-version=7")
  message("Setting -fabi-version for GCC 4.x")
endif()

set(DISPATCH_OBJECTS
        "reconstructor/diamond/src/dp/swipe/banded_3frame_swipe.cpp"
        "reconstructor/diamond/src/search/stage1_2.cpp"
        "reconstructor/diamond/src/tools/benchmark.cpp"
        "reconstructor/diamond/src/dp/swipe/swipe_wrapper.cpp"
        "reconstructor/diamond/src/masking/tantan.cpp"
        "reconstructor/diamond/src/dp/scan_diags.cpp"
        "reconstructor/diamond/src/dp/ungapped_simd.cpp"
        "reconstructor/diamond/src/dp/pfscan/pfscan.cpp"
        "reconstructor/diamond/src/dp/swipe/anchored_wrapper.cpp"
        "reconstructor/diamond/src/dp/score_profile.cpp"
        )

if(EXTRA)
  LIST(APPEND DISPATCH_OBJECTS "reconstructor/diamond/src/tools/benchmark_swipe.cpp")
endif()

add_library(arch_generic OBJECT ${DISPATCH_OBJECTS})
target_compile_options(arch_generic PUBLIC -DDISPATCH_ARCH=ARCH_GENERIC -DARCH_ID=0 -DEigen=Eigen_GENERIC)
target_include_directories(arch_generic PRIVATE "${CMAKE_SOURCE_DIR}/reconstructor/diamond/src/lib")
if(X86)
  add_library(arch_sse4_1 OBJECT ${DISPATCH_OBJECTS})
  target_include_directories(arch_sse4_1 PRIVATE "${CMAKE_SOURCE_DIR}/reconstructor/diamond/src/lib")
  add_library(arch_avx2 OBJECT ${DISPATCH_OBJECTS})
  target_include_directories(arch_avx2 PRIVATE "${CMAKE_SOURCE_DIR}/reconstructor/diamond/src/lib")
  if(WITH_AVX512)
    add_library(arch_avx512 OBJECT ${DISPATCH_OBJECTS})
    target_include_directories(arch_avx512 PRIVATE "${CMAKE_SOURCE_DIR}/reconstructor/diamond/src/lib")
  endif()
  add_definitions(-DWITH_AVX2)
  add_definitions(-DWITH_SSE4_1)
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    target_compile_options(arch_sse4_1 PUBLIC -DDISPATCH_ARCH=ARCH_SSE4_1 -DARCH_ID=1 -D__SSSE3__ -D__SSE4_1__ -D__POPCNT__ -DEigen=Eigen_SSE4_1)
    target_compile_options(arch_avx2 PUBLIC -DDISPATCH_ARCH=ARCH_AVX2 -DARCH_ID=2 /arch:AVX2 -D__SSSE3__ -D__SSE4_1__ -D__POPCNT__ -DEigen=Eigen_AVX2)
    if(WITH_AVX512)
      target_compile_options(arch_avx512 PUBLIC -DDISPATCH_ARCH=ARCH_AVX512 -DARCH_ID=3 /arch:AVX512 -D__SSSE3__ -D__SSE4_1__ -D__POPCNT__ -DEigen=Eigen_AVX512)
    endif()
  else()
    target_compile_options(arch_sse4_1 PUBLIC -DDISPATCH_ARCH=ARCH_SSE4_1 -DARCH_ID=1 -mssse3 -mpopcnt -msse4.1 -DEigen=Eigen_SSE4_1)
    target_compile_options(arch_avx2 PUBLIC -DDISPATCH_ARCH=ARCH_AVX2 -DARCH_ID=2 -mssse3 -mpopcnt -msse4.1 -msse4.2 -mavx -mavx2 -DEigen=Eigen_AVX2)
    if(WITH_AVX512)
      target_compile_options(arch_avx512 PUBLIC -DDISPATCH_ARCH=ARCH_AVX512 -DARCH_ID=3 -mssse3 -mpopcnt -msse4.1 -msse4.2 -mavx -mavx2 -mavx512f -mavx512bw -DEigen=Eigen_AVX512)
    endif()
  endif()
endif(X86)

# NEON is mandatory on Aarch64
if(AARCH64)
  add_definitions(-DWITH_NEON)
  add_library(arch_neon OBJECT ${DISPATCH_OBJECTS})
  target_include_directories(arch_neon PRIVATE "${CMAKE_SOURCE_DIR}/reconstructor/diamond/src/lib")
  target_compile_options(arch_neon PUBLIC -DDISPATCH_ARCH=ARCH_NEON -DARCH_ID=4 -D__ARM_NEON -D__aarch64__ -DEigen=Eigen_NEON)
endif(AARCH64)

# NEON is optional on Armv7, so we need to check for compiler support,
# and for the <sys/auxv.h> header used for runtime detection.
if(ARM)
check_symbol_exists(getauxval "sys/auxv.h" HAVE_GETAUXVAL)
check_cxx_compiler_flag("-mfpu=neon" HAVE_MFPU_NEON)
if(HAVE_MFPU_NEON)
    add_definitions(-DWITH_NEON)
    add_definitions(-DHAVE_MFPU_NEON)
    add_library(arch_neon OBJECT ${DISPATCH_OBJECTS})
    target_include_directories(arch_neon PRIVATE "${CMAKE_SOURCE_DIR}/reconstructor/diamond/src/lib")
    target_compile_options(arch_neon PUBLIC -DDISPATCH_ARCH=ARCH_NEON -DARCH_ID=4 -D__ARM_NEON -DEigen=Eigen_NEON -mfpu=neon)
  endif()
  if(HAVE_GETAUXVAL)
    add_definitions(-DHAVE_GETAUXVAL)
  endif()
endif(ARM)

set(OBJECTS
        reconstructor/diamond/src/libdiamond.cpp
        reconstructor/diamond/src/run/main.cpp
        reconstructor/diamond/src/basic/config.cpp
        reconstructor/diamond/src/stats/score_matrix.cpp
        reconstructor/diamond/src/data/queries.cpp
        reconstructor/diamond/src/data/seed_histogram.cpp
        reconstructor/diamond/src/output/daa/daa_record.cpp
        reconstructor/diamond/src/util/command_line_parser.cpp
        reconstructor/diamond/src/util/seq_file_format.cpp
        reconstructor/diamond/src/util/util.cpp
        reconstructor/diamond/src/basic/basic.cpp
        reconstructor/diamond/src/basic/hssp.cpp
        reconstructor/diamond/src/dp/ungapped_align.cpp
        reconstructor/diamond/src/run/tools.cpp
        reconstructor/diamond/src/chaining/greedy_align.cpp
        reconstructor/diamond/src/output/output_format.cpp
        reconstructor/diamond/src/output/clustering_variables.cpp
        reconstructor/diamond/src/output/clustering_format.cpp
        reconstructor/diamond/src/output/join_blocks.cpp
        reconstructor/diamond/src/data/frequent_seeds.cpp
        reconstructor/diamond/src/align/legacy/query_mapper.cpp
        reconstructor/diamond/src/output/blast_tab_format.cpp
        reconstructor/diamond/src/output/blast_pairwise_format.cpp
        reconstructor/diamond/src/run/double_indexed.cpp
        reconstructor/diamond/src/output/sam_format.cpp
        reconstructor/diamond/src/align/align.cpp
        reconstructor/diamond/src/search/setup.cpp
        reconstructor/diamond/src/data/taxonomy.cpp
        reconstructor/diamond/src/masking/masking.cpp
        reconstructor/diamond/src/data/seed_set.cpp
        reconstructor/diamond/src/util/simd.cpp
        reconstructor/diamond/src/output/taxon_format.cpp
        reconstructor/diamond/src/output/daa/view.cpp
        reconstructor/diamond/src/output/output_sink.cpp
        reconstructor/diamond/src/output/target_culling.cpp
        reconstructor/diamond/src/align/legacy/banded_swipe_pipeline.cpp
        reconstructor/diamond/src/util/io/compressed_stream.cpp
        reconstructor/diamond/src/util/io/deserializer.cpp
        reconstructor/diamond/src/util/io/file_sink.cpp
        reconstructor/diamond/src/util/io/file_source.cpp
        reconstructor/diamond/src/util/io/input_file.cpp
        reconstructor/diamond/src/util/io/input_stream_buffer.cpp
        reconstructor/diamond/src/util/io/output_file.cpp
        reconstructor/diamond/src/util/io/output_stream_buffer.cpp
        reconstructor/diamond/src/util/io/serializer.cpp
        reconstructor/diamond/src/util/io/temp_file.cpp
        reconstructor/diamond/src/util/io/text_input_file.cpp
        reconstructor/diamond/src/data/taxon_list.cpp
        reconstructor/diamond/src/data/taxonomy_nodes.cpp
        reconstructor/diamond/src/util/algo/MurmurHash3.cpp
        reconstructor/diamond/src/search/stage0.cpp
        reconstructor/diamond/src/data/seed_array.cpp
        reconstructor/diamond/src/output/paf_format.cpp
        reconstructor/diamond/src/util/system/system.cpp
        reconstructor/diamond/src/util/algo/greedy_vertex_cover.cpp
        reconstructor/diamond/src/util/sequence/sequence.cpp
        reconstructor/diamond/src/tools/tsv_record.cpp
        reconstructor/diamond/src/tools/tools.cpp
        reconstructor/diamond/src/util/system/getRSS.cpp
        reconstructor/diamond/src/lib/tantan/LambdaCalculator.cc
        reconstructor/diamond/src/util/algo/edge_vec.cpp
        reconstructor/diamond/src/util/string/string.cpp
        reconstructor/diamond/src/align/extend.cpp
        reconstructor/diamond/src/test/test.cpp
        reconstructor/diamond/src/align/ungapped.cpp
        reconstructor/diamond/src/align/gapped_score.cpp
        reconstructor/diamond/src/align/gapped_final.cpp
        reconstructor/diamond/src/align/full_db.cpp
        reconstructor/diamond/src/align/culling.cpp
        reconstructor/diamond/src/cluster/cluster_registry.cpp
        reconstructor/diamond/src/cluster/cascaded/cascaded.cpp
        reconstructor/diamond/src/align/output.cpp
        reconstructor/diamond/src/tools/roc.cpp
        reconstructor/diamond/src/test/data.cpp
        reconstructor/diamond/src/test/test_cases.cpp
        reconstructor/diamond/src/chaining/smith_waterman.cpp
        reconstructor/diamond/src/output/xml_format.cpp
        reconstructor/diamond/src/align/gapped_filter.cpp
        reconstructor/diamond/src/util/parallel/filestack.cpp
        reconstructor/diamond/src/util/parallel/parallelizer.cpp
        reconstructor/diamond/src/util/parallel/multiprocessing.cpp
        reconstructor/diamond/src/tools/benchmark_io.cpp
        reconstructor/diamond/src/lib/alp/njn_dynprogprob.cpp
        reconstructor/diamond/src/lib/alp/njn_dynprogproblim.cpp
        reconstructor/diamond/src/lib/alp/njn_dynprogprobproto.cpp
        reconstructor/diamond/src/lib/alp/njn_ioutil.cpp
        reconstructor/diamond/src/lib/alp/njn_localmaxstat.cpp
        reconstructor/diamond/src/lib/alp/njn_localmaxstatmatrix.cpp
        reconstructor/diamond/src/lib/alp/njn_localmaxstatutil.cpp
        reconstructor/diamond/src/lib/alp/njn_random.cpp
        reconstructor/diamond/src/lib/alp/sls_alignment_evaluer.cpp
        reconstructor/diamond/src/lib/alp/sls_alp.cpp
        reconstructor/diamond/src/lib/alp/sls_alp_data.cpp
        reconstructor/diamond/src/lib/alp/sls_alp_regression.cpp
        reconstructor/diamond/src/lib/alp/sls_alp_sim.cpp
        reconstructor/diamond/src/lib/alp/sls_basic.cpp
        reconstructor/diamond/src/lib/alp/sls_pvalues.cpp
        reconstructor/diamond/src/align/global_ranking/global_ranking.cpp
        reconstructor/diamond/src/align/global_ranking/extend.cpp
        reconstructor/diamond/src/tools/rocid.cpp
        reconstructor/diamond/src/lib/blast/blast_seg.cpp
        reconstructor/diamond/src/lib/blast/blast_filter.cpp
        reconstructor/diamond/src/lib/blast/nlm_linear_algebra.cpp
        reconstructor/diamond/src/stats/matrices/blosum45.cpp
        reconstructor/diamond/src/stats/matrices/blosum50.cpp
        reconstructor/diamond/src/stats/matrices/blosum62.cpp
        reconstructor/diamond/src/stats/matrices/blosum80.cpp
        reconstructor/diamond/src/stats/matrices/blosum90.cpp
        reconstructor/diamond/src/stats/matrices/pam250.cpp
        reconstructor/diamond/src/stats/matrices/pam30.cpp
        reconstructor/diamond/src/stats/matrices/pam70.cpp
        reconstructor/diamond/src/stats/stats.cpp
        reconstructor/diamond/src/stats/cbs.cpp
        reconstructor/diamond/src/stats/comp_based_stats.cpp
        reconstructor/diamond/src/stats/hauser_correction.cpp
        reconstructor/diamond/src/stats/matrix_adjust.cpp
        reconstructor/diamond/src/stats/matrix_adjust_eigen.cpp
        reconstructor/diamond/src/data/index.cpp
        reconstructor/diamond/src/data/dmnd/dmnd.cpp
        reconstructor/diamond/src/data/sequence_file.cpp
        reconstructor/diamond/src/tools/find_shapes.cpp
        reconstructor/diamond/src/data/block/block.cpp
        reconstructor/diamond/src/data/block/block_wrapper.cpp
        reconstructor/diamond/src/run/config.cpp
        reconstructor/diamond/src/data/sequence_set.cpp
        reconstructor/diamond/src/align/global_ranking/table.cpp
        reconstructor/diamond/src/output/daa/daa_write.cpp
        reconstructor/diamond/src/search/seed_complexity.cpp
        reconstructor/diamond/src/tools/view.cpp
        reconstructor/diamond/src/util/tsv/tsv.cpp
        reconstructor/diamond/src/tools/join.cpp
        reconstructor/diamond/src/basic/value.cpp
        reconstructor/diamond/src/masking/motifs.cpp
        reconstructor/diamond/src/align/alt_hsp.cpp
        reconstructor/diamond/src/data/fasta/fasta_file.cpp
        reconstructor/diamond/src/cluster/output.cpp
        reconstructor/diamond/src/cluster/realign.cpp
        reconstructor/diamond/src/cluster/reassign.cpp
        reconstructor/diamond/src/util/tsv/read_tsv.cpp
        reconstructor/diamond/src/tools/greedy_vertex_cover.cpp
        reconstructor/diamond/src/cluster/cascaded/recluster.cpp
        reconstructor/diamond/src/cluster/helpers.cpp
        reconstructor/diamond/src/search/kmer_ranking.cpp
        reconstructor/diamond/src/chaining/hamming_ext.cpp
        reconstructor/diamond/src/dna/smith_watermann.cpp
        reconstructor/diamond/src/stats/dna_scoring/build_score.cpp
        reconstructor/diamond/src/lib/blast/blast_message.cpp
        reconstructor/diamond/src/lib/blast/blast_stat.cpp
        reconstructor/diamond/src/lib/blast/blastn_score.cpp
        reconstructor/diamond/src/lib/blast/ncbi_std.cpp
        reconstructor/diamond/src/util/tsv/table.cpp
        reconstructor/diamond/src/util/tsv/file.cpp
        reconstructor/diamond/src/util/tsv/record.cpp
        reconstructor/diamond/src/cluster/cascaded/helpers.cpp
        reconstructor/diamond/src/cluster/cascaded/wrapper.cpp
        reconstructor/diamond/src/output/daa/merge.cpp
        reconstructor/diamond/src/chaining/backtrace.cpp
        reconstructor/diamond/src/util/tsv/merge.cpp
        reconstructor/diamond/src/util/tsv/join.cpp
        reconstructor/diamond/src/dp/scalar/smith_waterman.cpp
        reconstructor/diamond/src/align/short.cpp
        reconstructor/diamond/src/tools/tsv.cpp
        )

if(WITH_DNA)
  list(APPEND OBJECTS
          reconstructor/diamond/src/dna/smith_watermann.cpp
          reconstructor/diamond/src/dna/dna_index.cpp
          reconstructor/diamond/src/dna/seed_set_dna.cpp
          reconstructor/diamond/src/dna/extension.cpp
          reconstructor/diamond/src/dna/chain.cpp
          reconstructor/diamond/src/dna/extension_chain.cpp
          #reconstructor/diamond/src/dna/align.cpp
          reconstructor/diamond/src/lib/ksw2/ksw2_extz2_sse.c
          reconstructor/diamond/src/lib/ksw2/ksw2_extz.c
          reconstructor/diamond/src/lib/WFA2-lib.diamond/bindings/cpp/WFAligner.cpp
  )
endif()

if(WITH_MCL)
  list(APPEND OBJECTS reconstructor/diamond/src/contrib/mcl/mcl.cpp)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

if(BLAST_INCLUDE_DIR)
  set(BLAST_OBJ "reconstructor/diamond/src/data/blastdb/blastdb.cpp")
endif()

if(WITH_ZSTD OR BLAST_INCLUDE_DIR)
  set(ZSTD_OBJ "reconstructor/diamond/src/util/io/zstd_stream.cpp")
endif()

if(X86)
  if(WITH_AVX512)
    add_library(libdiamond MODULE $<TARGET_OBJECTS:arch_generic> $<TARGET_OBJECTS:arch_sse4_1> $<TARGET_OBJECTS:arch_avx2> $<TARGET_OBJECTS:arch_avx512> ${OBJECTS} ${BLAST_OBJ} ${ZSTD_OBJ})
  else()
    add_library(libdiamond MODULE $<TARGET_OBJECTS:arch_generic> $<TARGET_OBJECTS:arch_sse4_1> $<TARGET_OBJECTS:arch_avx2> ${OBJECTS} ${BLAST_OBJ} ${ZSTD_OBJ})
  endif()
elseif(ARM OR AARCH64)
  add_library(libdiamond MODULE $<TARGET_OBJECTS:arch_generic> $<TARGET_OBJECTS:arch_neon> ${OBJECTS} ${BLAST_OBJ} ${ZSTD_OBJ})
else()
  add_library(libdiamond MODULE $<TARGET_OBJECTS:arch_generic> ${OBJECTS} ${BLAST_OBJ} ${ZSTD_OBJ})
endif()

target_include_directories(libdiamond PRIVATE
        "${ZLIB_INCLUDE_DIR}"
        "${CMAKE_SOURCE_DIR}/reconstructor/diamond/src/lib")

if(BLAST_INCLUDE_DIR)
  function(find_blast_lib var library)
    find_library(${var} ${library} PATHS ${BLAST_LIBRARY_DIR} "/usr/lib/ncbi-blast+")
    if(${${var}} STREQUAL "${var}-NOTFOUND")
      message(FATAL_ERROR "Unable to find BLAST library: ${library}")
    endif()
    message("Found BLAST library: ${${var}}")
  endfunction(find_blast_lib)
  message("BLAST_INCLUDE_DIR: ${BLAST_INCLUDE_DIR}")
  message("BLAST_LIBRARY_DIR: ${BLAST_LIBRARY_DIR}")
  target_include_directories(libdiamond PRIVATE "${BLAST_INCLUDE_DIR}")
  find_blast_lib(SEQDB_LIBRARY seqdb)
  find_blast_lib(BIBLIO_LIBRARY biblio)
  find_blast_lib(BLASTDB_LIBRARY blastdb)
  find_blast_lib(BLASTDB_FORMAT_LIBRARY blastdb_format)
  find_blast_lib(GENERAL_LIBRARY general)
  find_blast_lib(GENOME_COLLECTION_LIBRARY genome_collection)
  find_blast_lib(MEDLINE_LIBRARY medline)
  find_blast_lib(PUB_LIBRARY pub)
  find_blast_lib(SEQ_LIBRARY seq)
  find_blast_lib(SEQCODE_LIBRARY seqcode)
  find_blast_lib(SEQSET_LIBRARY seqset)
  find_blast_lib(SEQUTIL_LIBRARY sequtil)
  find_blast_lib(SUBMIT_LIBRARY submit)
  find_blast_lib(XOBJMGR_LIBRARY xobjmgr)
  find_blast_lib(XOBJUTIL_LIBRARY xobjutil)
  find_blast_lib(XSER_LIBRARY xser)
  find_blast_lib(XNCBI_LIBRARY xncbi)
  find_blast_lib(XUTIL_LIBRARY xutil)
  find_blast_lib(LMDB_LIBRARY lmdb)
  target_link_libraries(libdiamond ${SEQDB_LIBRARY} ${BLASTDB_LIBRARY} ${BLASTDB_FORMAT_LIBRARY} ${GENERAL_LIBRARY}
          ${XOBJUTIL_LIBRARY} ${XOBJMGR_LIBRARY} ${GENOME_COLLECTION_LIBRARY} ${SEQ_LIBRARY} ${SEQCODE_LIBRARY}
          ${SEQSET_LIBRARY} ${SEQUTIL_LIBRARY} ${PUB_LIBRARY} ${MEDLINE_LIBRARY} ${BIBLIO_LIBRARY} ${SUBMIT_LIBRARY}
          ${XSER_LIBRARY} ${XNCBI_LIBRARY} ${XUTIL_LIBRARY} ${LMDB_LIBRARY})
  add_definitions(-DWITH_BLASTDB)
  set_target_properties(libdiamond PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    find_library(DBGHELP_LIBRARY Dbghelp)
    target_link_libraries(libdiamond ${DBGHELP_LIBRARY})
  else()
    target_link_libraries(libdiamond ${CMAKE_DL_LIBS})
  endif()
endif()

if(WITH_DNA)
  include_directories(reconstructor/diamond/src/lib/WFA2-lib.diamond)
  add_subdirectory(reconstructor/diamond/src/lib/WFA2-lib.diamond)
  target_link_libraries(libdiamond wfa2)
endif()

if(EIGEN_BLAS)
  add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE)
  find_package(LAPACKE)
  target_include_directories(libdiamond PRIVATE "${LAPACKE_INCLUDE_DIR}")
  target_link_libraries(libdiamond "${LAPACKE_LIBRARIES_DEP}")
endif()

if(WITH_ZSTD OR BLAST_INCLUDE_DIR)
  find_path(${ZSTD_INCLUDE_DIR} NAMES zstd.h PATHS ${ZSTD_INCLUDE_DIR})
  find_library(ZSTD_LIBRARY NAMES libzstd.a libzstd_static)
  target_include_directories(libdiamond PRIVATE "${ZSTD_INCLUDE_DIR}")
  target_link_libraries(libdiamond ${ZSTD_LIBRARY})
  add_definitions(-DWITH_ZSTD)
endif()

if(WITH_MIMALLOC)
  find_package(mimalloc 2.0 REQUIRED)
  add_definitions(-DWITH_MIMALLOC)
  target_link_libraries(libdiamond mimalloc)
  target_include_directories(libdiamond PRIVATE ${MIMALLOC_INCLUDE_DIR})
endif()

find_package(PythonExtensions REQUIRED)
python_extension_module(libdiamond)
set_target_properties(libdiamond PROPERTIES OUTPUT_NAME "libdiamond")
install(TARGETS libdiamond LIBRARY DESTINATION .)