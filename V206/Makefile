all: build/V206.pdf

build/V206.pdf: FORCE V206.tex V206.py build/plotT1.pdf build/plotT2.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V206.tex

build/plotT1.pdf: V206.py data.txt |build

	TEXINPUTS=$$(pwd): python V206.py

build/plotT2.pdf: V206.py data.txt |build

	TEXINPUTS=$$(pwd): python V206.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
