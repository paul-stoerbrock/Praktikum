all: build/V903.pdf

build/V903.pdf: FORCE V903.tex V903.py build/plot15.pdf build/plot30.pdf build/plot45.pdf build/plotdepth.pdf build/plotsignal.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V903.tex

build/plot15.pdf:  V903.py 15degree.txt |build

	TEXINPUTS=$$(pwd): python V903.py

build/plot30.pdf:  V903.py 30degree.txt |build

	TEXINPUTS=$$(pwd): python V903.py

build/plot45.pdf:  V903.py 45degree.txt |build

	TEXINPUTS=$$(pwd): python V903.py

build/plotdepth.pdf:  V903.py data.txt |build

	TEXINPUTS=$$(pwd): python V903.py

build/plotsignal.pdf:  V903.py data.txt |build

	TEXINPUTS=$$(pwd): python V903.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean