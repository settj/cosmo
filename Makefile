# NOTE: needs boost, tclap, STXXL, KMC, and sdsl

CXX=g++
CPP_FLAGS=-pipe -m64 -std=c++14 -pedantic-errors -W -Wall -Wextra -Wpointer-arith \
					-Wunused -Wwrite-strings -g #-Wcast-qual #-Werror

CLANG_WARNINGS=-Wbool-conversions -Wshift-overflow -Wliteral-conversion # CLANG ONLY
BOOST_PATH=3rd_party_inst/boost

DEP_PATH=3rd_party_inst
KMC_PATH=./3rd_party_src/KMC
SDREADER_PATH=3rd_party_src/sdreader
INC=-isystem $(DEP_PATH)/include -isystem $(BOOST_PATH)/include -I$(SDREADER_PATH)
LIB=-L$(DEP_PATH)/lib -L./ -L$(BOOST_PATH)/lib -L$(SDREADER_PATH)
BOOST_FLAGS= -lboost_system -lboost_filesystem
DEP_FLAGS=$(INC) $(LIB) $(BOOST_FLAGS) -isystem $(KMC_PATH) -lsdsl -fopenmp #openmp is needed for logging
DEBUG_FLAGS= -g
NDEBUG_FLAGS=-DNDEBUG
OPT_FLAGS=-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing
NOPT_FLAGS=-O0

# Using Semantic Versioning: http://semver.org/
VERSION=0.7.0
BANNER='Copyright Alex Bowe (c) 2016'
CPP_FLAGS+=-DVERSION=\"$(VERSION)\" -DBANNER=\"$(BANNER)\"

k?=32
# TODO: quantize k
CPP_FLAGS+=-DK_LEN=$(k)

colors?=4000
CPP_FLAGS+=-DNUM_COLS=$(colors)

ifeq ($(asm),1)
CPP_FLAGS+=-S -fverbose-asm
endif

ifeq ($(optimise),0)
CPP_FLAGS+=$(NOPT_FLAGS)
else
CPP_FLAGS+=$(OPT_FLAGS)
endif

ifeq ($(debug),1)
CPP_FLAGS+=$(DEBUG_FLAGS)
else
CPP_FLAGS+=$(NDEBUG_FLAGS)
endif

ifeq ($(verbose),1)
CPP_FLAGS+=-DVERBOSE
endif

SDREADER_OBJS=$(SDREADER_PATH)/LowReader.o $(SDREADER_PATH)/HighReader.o $(SDREADER_PATH)/SDIter.o
KMC_OBJS=$(KMC_PATH)/kmc_api/kmc_file.o $(KMC_PATH)/kmc_api/kmer_api.o $(KMC_PATH)/kmc_api/mmer.o
BUILD_REQS=lut.hpp debug.hpp utility.hpp io.hpp sort.hpp kmer.hpp dummies.hpp debruijn_graph.hpp debruijn_graph_shifted.hpp 
COLOR_REQS=debruijn_graph_shifted.hpp cosmo-color.hpp kmer.hpp
DUMP_REQS=debruijn_graph_shifted.hpp
MERGE_REQS=debug.hpp utility.hpp io.hpp sort.hpp kmer.hpp dummies.hpp 
COLOR_MERGE_REQS=kmer.hpp
COUNT_REQS=io.hpp
VARI_MERGE_REQS=config.hpp debruijn_graph_shifted.hpp sort.hpp
TEST_VARI_DELETE_REQS=sort.hpp
DYN_TEST_REQS=kmer-counter.hpp dyn_debruijn_graph.hpp 
DYN_IN_OUT_TEST_REQS=debruijn_graph_shifted.hpp dyn_debruijn_graph.hpp 
UNI_REQS=debruijn_graph_shifted.hpp
COLOR_PD_REQS=debruijn_graph_shifted.hpp
READ_COLOR_REQS=debruijn_graph_shifted.hpp
TRANSPOSE_REQS=debruijn_graph_shifted.hpp
PACK_COLOR_REQS=kmer.hpp
PACK_COLOR_MERGE_REQS=''
BENCHMARK_REQS=debruijn_graph_shifted.hpp debruijn_hypergraph.hpp
TEST_REQS=kmer.hpp utility.hpp debug.hpp multi_bit_vector.hpp bgl_sdb_adapter.hpp 
BINARIES=cosmo-build cosmo-color cosmo-test cosmo-benchmark cosmo-benchmark-varord pack-color cosmo-color-pd cosmo-read-color transpose cosmo-merge vari-merge cosmo-dump cosmo-dump-full-edges color-merge vari-uni count_kmer vari-delete test-vari-delete dyn_test dyn_in_out_test dyn_add_test

default: all

# Stores last compiled options to know whether we need to recompile when comp variables change
.PHONY: force
compiler_flags: force
	@echo $(CPP_FLAGS) | cmp -s - $@ || echo $(CPP_FLAGS) > $@

lut.hpp: make_lut.py
	python make_lut.py > lut.hpp

cosmo-build: cosmo-build.cpp $(BUILD_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS) -lstxxl

cosmo-color: cosmo-color.cpp $(COLOR_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

cosmo-dump: cosmo-dump.cpp $(DUMP_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

cosmo-dump-full-edges: cosmo-dump-full-edges.cpp $(DUMP_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

cosmo-merge: cosmo-merge.cpp $(MERGE_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS) 

color-merge: color-merge.cpp $(COLOR_MERGE_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(SDREADER_OBJS) $(DEP_FLAGS) 

count_kmer: count_kmer.cpp $(COUNT_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

vari-merge: vari-merge.cpp $(VARI_MERGE_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS) -lstxxl

vari-delete: vari-delete.cpp $(VARI_MERGE_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(SDREADER_OBJS) $(DEP_FLAGS) -lstxxl

test-vari-delete: test-vari-delete.cpp $(TEST_VARI_DELETE_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(SDREADER_OBJS) $(DEP_FLAGS) -lstxxl

dyn_add_test: dyn_add_test.cpp $(DYN_TEST_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(SDREADER_OBJS) $(DEP_FLAGS) -lstxxl

dyn_test: dyn_test.cpp $(DYN_TEST_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(SDREADER_OBJS) $(DEP_FLAGS) -lstxxl

dyn_in_out_test: dyn_in_out_test.cpp $(DYN_IN_OUT_TEST_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(SDREADER_OBJS) $(DEP_FLAGS) -lstxxl

vari-uni: vari-uni.cpp $(UNI_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

cosmo-color-pd: cosmo-color-pd.cpp $(COLOR_PD_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

cosmo-read-color: cosmo-read-color.cpp $(READ_COLOR_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

transpose: transpose.cpp $(TRANSPOSE_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

pack-color: pack-color.cpp $(PACK_COLOR_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

pack-color-merge: pack-color-merge.cpp $(PACK_COLOR_MERGE_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

cosmo-benchmark: cosmo-benchmark.cpp $(BENCHMARK_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS)

cosmo-benchmark-varord: cosmo-benchmark.cpp $(BENCHMARK_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) -DVAR_ORDER

#cosmo-test: cosmo-test.cpp catch.hpp $(wildcard *_test.cpp) $(wildcard $(subst _test.cpp,.hpp,$(wildcard *_test.cpp)))
#	$(CXX) $(CPP_FLAGS) -o $@ $(filter-out %.hpp,$^) $(DEP_FLAGS) -lstxxl -fopenmp -lsdsl

cosmo-test: debruijn_graph_test.cpp $(TEST_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) -lboost_unit_test_framework

test: cosmo-test
	./cosmo-test

all: $(BINARIES)

clean:
		rm -rf $(BINARIES) *.o *.dSYM
