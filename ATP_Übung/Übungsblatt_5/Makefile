all: build/Übungsblatt_5.pdf

build/Übungsblatt_5.pdf: FORCE Übungsblatt_5.tex Übungsblatt_5.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	python Übungsblatt_5.py
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	Übungsblatt_5.tex
	

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
