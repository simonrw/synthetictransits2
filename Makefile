CPPSOURCES := $(wildcard src/*.cpp)
CCSOURCES := $(wildcard src/*.C)
OBJECTS := $(CPPSOURCES:.cpp=.o) $(CCSOURCES:.C=.o)
RUN := bin/synthetic_transits
COMMON := -Wno-write-strings -O0 -g -std=c++11
CFLAGS := -I/usr/local/cfitsio/include -I/usr/local/tclap/include -Iinclude -Imodelgen/include -Irgwtimer/include
LDFLAGS := -L/usr/local/cfitsio/lib -lcfitsio -lsqlite3

all: $(RUN)

$(RUN): bin $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ ${COMMON} ${LDFLAGS}

%.o: %.cpp
	$(CXX) -c $< -o $@ ${COMMON} ${CFLAGS}

%.o: %.C
	$(CXX) -c $< -o $@ ${COMMON} ${CFLAGS}

bin:
	mkdir -p $@

clean:
	@-find . -name '*.o' -delete
	@-rm $(RUN)

test: $(RUN)
	$(RUN) -o out.fits -c models.db -i NG0522-2518.fits

gdb: $(RUN)
	gdb --args $(RUN) -o out.fits -c models.db -i NG0522-2518.fits

.PHONY: clean test gdb
