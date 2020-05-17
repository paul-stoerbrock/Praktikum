all: build/V602_P.pdf

build/V602_P.pdf: FORCE V602_P.tex V602_P.py build/plotBr.pdf build/plotBragg.pdf build/plotCu.pdf build/plotE_K.pdf build/plotGa.pdf build/plotRb.pdf build/plotSr.pdf build/plotZn.pdf build/plotZr.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V602_P.tex

build/plotBr.pdf:  V602_P.py Brom.dat |build

	TEXINPUTS=$$(pwd): python V602_P.py

build/plotBragg.pdf:  V602_P.py Bragg.dat |build

	TEXINPUTS=$$(pwd): python V602_P.py

build/plotCu.pdf:  V602_P.py Emissionsspektrum.dat |build

	TEXINPUTS=$$(pwd): python V602_P.py

build/plotE_K.pdf:  V602_P.py |build

	TEXINPUTS=$$(pwd): python V602_P.py

build/plotGa.pdf:  V602_P.py Gallium.dat |build

	TEXINPUTS=$$(pwd): python V602_P.py

build/plotRb.pdf:  V602_P.py Rubidium.dat |build

	TEXINPUTS=$$(pwd): python V602_P.py

build/plotSr.pdf:  V602_P.py Strontium.dat |build

	TEXINPUTS=$$(pwd): python V602_P.py

build/plotZn.pdf:  V602_P.py Zink.dat |build

	TEXINPUTS=$$(pwd): python V602_P.py

build/plotZr.pdf:  V602_P.py Zirkonium.dat |build

	TEXINPUTS=$$(pwd): python V602_P.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean