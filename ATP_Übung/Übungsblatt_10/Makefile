all: build/Übungsblatt_10.pdf

build/Übungsblatt_10.pdf: FORCE Übungsblatt_10.tex Übungsblatt_10.py   |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	python Übungsblatt_10.py
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	Übungsblatt_10.tex




build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean