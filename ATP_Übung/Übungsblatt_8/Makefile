all: build/Übungsblatt_8.pdf

build/Übungsblatt_8.pdf: FORCE Übungsblatt_8.tex Übungsblatt_8.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	python Übungsblatt_8.py
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	Übungsblatt_8.tex
	

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
