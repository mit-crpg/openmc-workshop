slides.pdf: src/slides.tex
	mkdir -p build
	pdflatex --output-directory=build src/slides.tex
	pdflatex --output-directory=build src/slides.tex
	cp build/slides.pdf .

clean:
	rm -rf build
