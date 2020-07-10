all: build/V601.pdf

build/V601.pdf: FORCE V601.tex V601.py build/plot_UA.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V601.tex

build/plot_UA.pdf:  V601.py dataUA_T1.txt |build

	TEXINPUTS=$$(pwd): python V601.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean