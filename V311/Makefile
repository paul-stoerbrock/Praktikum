all: build/V311.pdf

build/V311.pdf: FORCE V311.tex V311.py build/plothy.pdf|build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V311.tex

build/plothy.pdf: V311.py|build
	python V311.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
