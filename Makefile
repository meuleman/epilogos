BINDIR = bin
SRCDIR = src
CXX = g++
CXXFLAGS = -O3 -pedantic -Wall -ansi -static

TARGETS = computeEpilogosPart1_perChrom computeEpilogosPart2_perChrom computeEpilogosPart3_perChrom
EXE = $(addprefix $(BINDIR)/,$(TARGETS))

default: $(EXE)

$(BINDIR)/% : $(SRCDIR)/%.cpp
	mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(EXE)
