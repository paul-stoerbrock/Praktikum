all: build/V500.pdf

build/V500.pdf: FORCE V500.tex V500.py build/plotRot.pdf build/plotGelb.pdf build/plotGruen.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V500.tex

build/plotRot.pdf:  V500.py Rot.txt |build

	TEXINPUTS=$$(pwd): python V500.py

build/plotGelb.pdf:  V500.py Gelb.txt |build

	TEXINPUTS=$$(pwd): python V500.py

build/plotGruen.pdf:  V500.py Gruen.txt |build

	TEXINPUTS=$$(pwd): python V500.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean