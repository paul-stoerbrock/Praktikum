all: build/V353.pdf

build/V353.pdf: FORCE V353.tex build/plot.pdf build/mean_aRC.tex build/mean_bRC.tex build/mean_cRC.tex |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V353.tex

build/plot.pdf: V353.py data.txt entladung.txt true.txt|build
	
	TEXINPUTS=$$(pwd): python V353.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
