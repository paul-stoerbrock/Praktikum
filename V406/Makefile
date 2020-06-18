all: build/V406.pdf

build/V406.pdf: FORCE V406.tex V406.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V406.tex

#build/plot.pdf:  V406.py .dat |build
#
#	TEXINPUTS=$$(pwd): python V406.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean