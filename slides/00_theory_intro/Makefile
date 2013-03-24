slides.pdf: src/slides.tex
	mkdir -p build
	cd build; pdflatex --shell-escape ../src/slides.tex; pdflatex --shell-escape ../src/slides.tex; cp slides.pdf ..
	cd ..

clean:
	rm -rf build
