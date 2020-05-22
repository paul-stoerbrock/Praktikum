all: build/V703.pdf

build/V703.pdf: FORCE V703_T.tex V703_T.py build/plot_kenn.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V703_T.tex

build/plot_kenn.pdf: V703_T.py Kennlinie.dat |build

	TEXINPUTS=$$(pwd): python V703_T.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean