all: build/Übungsblatt_9.pdf

build/Übungsblatt_9.pdf: FORCE Übungsblatt_9.tex Übungsblatt_9.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	Übungsblatt_9.tex

#build/plot_d.pdf: Übungsblatt_9.py |build
#
#	TEXINPUTS=$$(pwd): python Übungsblatt_9.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
