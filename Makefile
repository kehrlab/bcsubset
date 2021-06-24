# Use version 4.9 of g++

# Set dependencies, s.t. make detects changes in header files.
%.o: %.cc %.h
	g++ -c $^

CXX=g++
CC=$(CXX)

SEQAN_LIB=.
SPOA_LIB=.

CXXFLAGS+=-I$(SEQAN_LIB) -I$(SPOA_LIB) -DSEQAN_HAS_ZLIB=1 -std=c++14 -DSEQAN_DISABLE_VERSION_CHECK
LDLIBS=-lz -lpthread /usr/local/lib/libspoa.a

DATE=on 2021-06-24
VERSION=0.0.1
CXXFLAGS+=-DDATE=\""$(DATE)"\" -DVERSION=\""$(VERSION)"\"

# Enable warnings, disable some
CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result -Wno-deprecated-copy -Wno-class-memaccess

HEADERS=argparse.h bamsubset.h workflow.h

.PHONY: all
all: CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
all: bcsubset

.PHONY: debug
debug: CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1
debug: bcsubset

PREFIX = /usr/local
.PHONY: install
install: CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
install: bcsubset
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp $< $(DESTDIR)$(PREFIX)/bin/bcsubset

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/bcsubset

ctprocess: bcsubset.o

main.o: bcsubset.cpp $(HEADERS)

.PHONY: clean
clean:
	rm -f *.o bcsubset
