all: build/V602_P.pdf

build/V602_P.pdf: FORCE V602_P.tex V602_P.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V602_P.tex

#build/plot-example.pdf:  V602_P.py .dat |build
#
#	TEXINPUTS=$$(pwd): python V602_P.py
#
#build/plot-example.pdf:  V602_P.py .dat |build
#
#	TEXINPUTS=$$(pwd): python V602_P.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean