all: build/V354.pdf

build/V354.pdf: FORCE V354.tex build/plotc.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V354.tex

build/plotc.pdf: V354.py data.txt |build
	
	TEXINPUTS=$$(pwd): python V354.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
