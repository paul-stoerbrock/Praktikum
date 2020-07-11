all: build/V703_P.pdf

build/V703_P.pdf: FORCE V703_P.tex V703_P.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V703_P.tex

build/plot.pdf:  V703_P.py .dat |build

	TEXINPUTS=$$(pwd): python V703_P.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean