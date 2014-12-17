ROOTCONFIG   := root-config

ARCH         := $(shell $(ROOTCONFIG) --arch)
PLATFORM     := $(shell $(ROOTCONFIG) --platform)

CXX          :=  g++

ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
HASTHREAD    := $(shell $(ROOTCONFIG) --has-thread)

LDFLAGS      += $(ROOTLDFLAGS)

LIBS = $(ROOTLIBS)  -L$(RKLIBS) -L/usr/lib64  -lgfortran 

MYINCLUDEFLAGS = -I$(RKLIBS)/../include  
CXXFLAGS = -Wall -g $(ROOTCFLAGS) $(MYINCLUDEFLAGS)
LINKFLAGS = -g $(LIBS)

%.o : %.cc
	$(CXX) $(CXXFLAGS) -MD -c $<
	@sed -i 's,\($*\.o\)[:]*\(.*\),$@: $$\(wildcard\2\)\n\1:\2,g' $*.d

all: main_makeClassFromRoot.o makeClassFromRoot 

main_makeClassFromRoot.o: main_makeClassFromRoot.cc
	$(CXX) $(CXXFLAGS) -c main_makeClassFromRoot.cc

makeClassFromRoot: main_makeClassFromRoot.o 
	$(CXX)  -o $@  main_makeClassFromRoot.o  $(LINKFLAGS) 


clean:
	rm -f *~ makeClassFromRoot main_makeClassFromRoot.o main_makeClassFromRoot.d 

-include main_makeClassFromRoot.d
