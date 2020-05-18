all: build/V602.pdf

build/V602.pdf: FORCE V602_T.tex V602_T.py build/plot_Cu.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V602_T.tex

build/plot_Cu.pdf: V602_T.py Emissionsspektrum.dat |build

	TEXINPUTS=$$(pwd): python V602_T.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean