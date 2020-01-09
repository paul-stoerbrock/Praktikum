all: build/V206.pdf

build/V203.pdf: FORCE V206.tex V206.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V206.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean