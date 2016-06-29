CPPSOURCES := $(wildcard src/*.cpp)
CCSOURCES := $(wildcard src/*.C)
BATMANDIR := $(shell pwd)/external/libbatman
OBJECTS := $(CPPSOURCES:.cpp=.o) $(CCSOURCES:.C=.o)
RUN := bin/synthetic_transits
COMMON := -Wno-write-strings -O2 -g -std=c++11
CFLAGS := -I/usr/local/cfitsio/include -I/usr/local/tclap/include -Iinclude -Imodelgen/include -Irgwtimer/include -I$(BATMANDIR)/c_src
LDFLAGS := -L/usr/local/cfitsio/lib -L$(BATMANDIR) -lcfitsio -lsqlite3 -lm -lbatman

all: $(RUN)

$(RUN): bin $(BATMANDIR)/libbatman.a $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ${COMMON} ${LDFLAGS}

$(BATMANDIR)/libbatman.a:
	$(MAKE) -C $(BATMANDIR) libbatman.a

%.o: %.cpp
	$(CXX) -c $< -o $@ ${COMMON} ${CFLAGS}

%.o: %.C
	$(CXX) -c $< -o $@ ${COMMON} ${CFLAGS}

%.o: %.c
	$(CC) -c $< -o $@ ${COMMON} ${CFLAGS}

bin:
	mkdir -p $@

clean:
	@-find . -name '*.o' -delete
	@-$(MAKE) -C $(BATMANDIR) clean
	@-rm $(RUN)

test: $(RUN)
	$(RUN) -o out-batman.fits -c ../testdata/MODELS_NG0522-2518_802_2016_TEST16.db \
		-i ../testdata/NG0522-2518.fits

gdb: $(RUN)
	gdb --args $(RUN) -o out.fits -c models.db -i NG0522-2518.fits

valgrind: $(RUN)
	valgrind $(RUN) -o out-batman.fits -c ../testdata/MODELS_NG0522-2518_802_2016_TEST16.db \
		-i ../testdata/NG0522-2518.fits


.PHONY: clean test gdb
