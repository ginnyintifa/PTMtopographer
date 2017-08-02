BINDIR=bin
SRCDIR=src
PROGS=program_feature_generation program_prediction_summary program_additional_annotation
PROG_FULLPATHS=$(addprefix $(BINDIR)/, $(PROGS))

all: $(PROG_FULLPATHS)

CC = $(CXX)
CXXFLAGS += -Wall -O2 -std=gnu++0x -fopenmp
LDFLAGS += -fopenmp
CPPFLAGS += -I .

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<

$(BINDIR)/%: $(SRCDIR)/%.o
	$(LINK.o) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	$(RM) $(PROG_FULLPATHS)
	$(RM) $(addsuffix .o, $(addprefix $(SRCDIR)/, $(PROGS)))


