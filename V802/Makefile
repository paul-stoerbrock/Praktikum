all: build/V802.pdf

build/V803.pdf: fourier.tex build/plot.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	fourier.tex

build/plot.pdf: fourier.py data.txt  |build
	
	python fourier.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean