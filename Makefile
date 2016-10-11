CPPSOURCES := src/main.cpp src/FitsObject.cpp src/GenerateModel.cpp src/fetches_parameters.cpp
CCSOURCES := src/timer.C
BATMANDIR := $(shell pwd)/external/libbatman
OBJECTS := $(CPPSOURCES:.cpp=.o) $(CCSOURCES:.C=.o)
RUN := bin/synthetic_transits
COMMON := -Wno-write-strings -O2 -g -std=c++11 -Wall -Wextra
CFLAGS := -I/usr/local/cfitsio/include -I/usr/local/tclap/include -Iinclude -Imodelgen/include -Irgwtimer/include -I$(BATMANDIR)/c_src
LDFLAGS := -L/usr/local/cfitsio/lib -L$(BATMANDIR) -lcfitsio -lsqlite3 -lm -lbatman

TESTOBJECTS := src/compare_models.o src/fetches_parameters.o src/GenerateModel.o src/FitsObject.o external/libbatman/c_src/light_curve.o \
	external/libbatman/c_src/_rsky.o external/libbatman/c_src/_nonlinear_ld.o external/libbatman/c_src/_eclipse.o

testdatadir := $(shell python ./find_testdata_dir.py)

all: $(RUN)

$(RUN): bin $(BATMANDIR)/libbatman.a $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ${COMMON} ${LDFLAGS}

bin/compare_models: bin $(TESTOBJECTS)
	$(CXX) -o $@ $(TESTOBJECTS) -O2 -g -L/usr/local/cfitsio/lib -lcfitsio -lsqlite3 -lm

$(BATMANDIR)/libbatman.a:
	$(MAKE) -C $(BATMANDIR) libbatman.a

%.o: %.cpp
	$(CXX) -c $< -o $@ ${COMMON} ${CFLAGS}

%.o: %.C
	$(CXX) -c $< -o $@ ${COMMON} ${CFLAGS}

%.o: %.c
	$(CC) -c $< -o $@ -O2 -g ${CFLAGS}

bin:
	mkdir -p $@

clean:
	@-find . -name '*.o' -delete
	@-$(MAKE) -C $(BATMANDIR) clean
	@-rm $(RUN)

test: $(RUN)
	$(RUN) -o $(outputdir)/out-test.fits -c $(testdatadir)/MODELS_NG0522-2518_802_2016_TEST16.db \
		-i $(testdatadir)/NG0522-2518.fits

singleobjecttest: $(RUN)
	$(RUN) -o $(outputdir)/out-test.fits -c $(testdatadir)/single_object.db \
		-i $(testdatadir)/NG0522-2518.fits

duplicateobjecttest: $(RUN)
	$(RUN) -o $(outputdir)/out-test.fits -c $(testdatadir)/duplicate_objects.db \
		-i $(testdatadir)/NG0522-2518.fits

gdb: $(RUN)
	gdb --args $(RUN) -o $(outputdir)/out-gdb.fits -c $(testdatadir)/MODELS_NG0522-2518_802_2016_TEST16.db \
		-i $(testdatadir)/NG0522-2518.fits

valgrind: $(RUN)
	valgrind --leak-check=full $(RUN) -o $(outputdir)/out-valgrind.fits -c $(testdatadir)/MODELS_NG0522-2518_802_2016_TEST16.db \
		-i $(testdatadir)/NG0522-2518.fits 2>&1 | tee valgrind.log

define outputdir
$(if $(findstring ngts10,$(shell hostname)),/local/srw/synthetic-testdata,.)
endef

.PHONY: clean test gdb
