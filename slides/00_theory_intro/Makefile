GNUPLOTSRC=$(wildcard *.plot)
slides.pdf: gnuplot src/slides.tex
	mkdir -p build
	cd build; pdflatex --shell-escape ../src/slides.tex; pdflatex --shell-escape ../src/slides.tex; cp slides.pdf ..

gnuplot: $(GNUPLOTSRC)
	cd src/src_figs; make

clean:
	rm -rf build; rm -f slides.pdf

distclean: clean
	cd src/src_figs; make clean

neat:
	rm -rf build

.PHONY: gnuplot clean distclean neat
