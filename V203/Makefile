all: build/V203.pdf

build/V203.pdf: FORCE V203.tex |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V203.tex

build/plotCuein.pdf: V203.py datak1b.txt datag1b.txt |build

	TEXINPUTS=$$(pwd): python V203.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean