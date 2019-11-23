all: build/V353.pdf

build/V353.pdf: FORCE V353.tex build/plot.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V353.tex

build/plot.pdf: V353.py data.txt entladung.txt |build
	
	TEXINPUTS=$$(pwd): python V353.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
