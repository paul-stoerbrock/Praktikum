all: build/Übungsblatt_8.pdf

build/Übungsblatt_8.pdf: FORCE Übungsblatt_8.tex Übungsblatt_8.py build/plot_d.pdf  |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	Übungsblatt_8.tex

build/plot_d.pdf: Übungsblatt_8.py |build

	TEXINPUTS=$$(pwd): python Übungsblatt_8.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean