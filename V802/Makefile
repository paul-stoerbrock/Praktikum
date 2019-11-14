all: build/V802.pdf

build/V802pdf: V802.tex build/plot.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V802.tex

build/plot.pdf: V802.py data.txt  |build
	
	python fourier.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean