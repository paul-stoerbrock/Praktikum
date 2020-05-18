all: build/Übungsblatt_3.pdf

build/Übungsblatt_3.pdf: FORCE Übungsblatt_3.tex Übungsblatt_3.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	python Übungsblatt_3.py
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	Übungsblatt_3.tex
	

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean