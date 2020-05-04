all: build/V603.pdf

build/V603.pdf: FORCE V603_T.tex V603_T.py build/plot_Cu.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V603_T.tex

build/plot_Cu.pdf: V603_T.py EmissionCu.dat |build

	TEXINPUTS=$$(pwd): python V603_T.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean