all: build/Übungsblatt_4.pdf

build/Übungsblatt_4.pdf: FORCE Übungsblatt_4.tex Übungsblatt_4.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	python Übungsblatt_4.py
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	Übungsblatt_4.tex
	

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
