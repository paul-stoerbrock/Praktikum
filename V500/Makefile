all: build/V500.pdf

build/V500.pdf: FORCE V500.tex V500.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V500.tex

#build/plot.pdf:  V500.py .dat |build
#
#	TEXINPUTS=$$(pwd): python V500.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean