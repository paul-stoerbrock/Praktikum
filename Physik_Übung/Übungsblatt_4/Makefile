all: build/Übungsblatt_4.pdf

build/Übungsblatt_4.pdf: FORCE Übungsblatt_4.tex |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	Übungsblatt_4.tex

#build/plot1.pdf: Übungsblatt_4.py |build
#
#	TEXINPUTS=$$(pwd): python Übungsblatt_4.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean