all: build/Übungsblatt_3.pdf

build/Übungsblatt_3.pdf: FORCE Übungsblatt_3.tex build/plot.pdf build/plot1.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	Übungsblatt_3.tex

build/plot.pdf: Übungsblatt_3.py |build

	TEXINPUTS=$$(pwd): python Übungsblatt_3.py

build/plot1.pdf: Übungsblatt_3.py |build

	TEXINPUTS=$$(pwd): python Übungsblatt_3.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean