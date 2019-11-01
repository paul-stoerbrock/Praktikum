all: build/V803.pdf

build/V803.pdf: V803.tex build/plot.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V803.tex

build/plot.pdf: V803-1.py data.txt  |build
	
	python V803-1.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
