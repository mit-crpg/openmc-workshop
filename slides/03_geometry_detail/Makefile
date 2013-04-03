
slides = src/slides.tex \
src/univs.tex \
src/lattices.tex \
src/core.tex \
src/plotting.tex

slides.pdf: ${slides}
	mkdir -p build
	pdflatex --output-directory=build src/slides.tex
	pdflatex --output-directory=build src/slides.tex
	cp build/slides.pdf .

clean:
	rm -rf build
