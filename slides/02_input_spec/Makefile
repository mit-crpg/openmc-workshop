slides.pdf: src/slides.tex beamercolorthememit.sty
	mkdir -p build
	pdflatex --output-directory=build src/slides.tex
	pdflatex --output-directory=build src/slides.tex
	cp build/slides.pdf .
	rm -f beamercolorthememit.sty

beamercolorthememit.sty: src/beamercolorthememit.sty
	cp src/$@ $@

clean:
	rm -rf build
