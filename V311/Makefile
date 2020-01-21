all: build/V311.pdf

build/V311.pdf: FORCE V311.tex V311.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V311.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
