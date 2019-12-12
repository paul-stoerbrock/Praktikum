all: build/V103.pdf

build/V103.pdf: FORCE V103.tex build/plotCuein.pdf build/Cu_ein.tex |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V103.tex

build/plotCuein.pdf: V103.py dataCuein.txt |build

	TEXINPUTS=$$(pwd): python V103.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean