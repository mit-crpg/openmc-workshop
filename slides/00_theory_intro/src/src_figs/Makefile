templates = $(wildcard *.plot)
objects = $(templates:.plot=.pdf)

all: $(objects)
clean:
	rm -f *.pdf

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .plot .pdf
.PHONY: all clean

%.pdf: %.plot
	gnuplot $^
