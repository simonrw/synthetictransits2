CPPSOURCES := $(wildcard src/*.cpp)
CCSOURCES := $(wildcard src/*.C)
OBJECTS := $(CPPSOURCES:.cpp=.o) $(CCSOURCES:.C=.o)
RUN := bin/synthetic_transits
COMMON := -Wno-write-strings
CFLAGS := -Iinclude -Imodelgen/include -Irgwtimer/include
LDFLAGS := -lcfitsio

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

.PHONY: clean
