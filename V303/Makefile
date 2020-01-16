all: build/V303.pdf

build/V303.pdf: FORCE V303.tex V303.py build/plotphi_oS.pdf  |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V303.tex

build/plotphi_oS.pdf: V303.py aufg2.txt |build

	TEXINPUTS=$$(pwd): python V303.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
