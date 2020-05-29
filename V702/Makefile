all: build/V702.pdf

build/V702.pdf: FORCE V702.tex V702.py build/plot_V.pdf  |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V702.tex

build/plot_V.pdf: V702.py Vanadium.dat |build

	TEXINPUTS=$$(pwd): python V702.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean