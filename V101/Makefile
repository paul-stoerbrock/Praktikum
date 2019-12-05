all: build/V101.pdf

build/V101.pdf: FORCE V101.tex build/plotD.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V101.tex

build/plotD.pdf: V101.py data.txt |build
	
	TEXINPUTS=$$(pwd): python V101.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean