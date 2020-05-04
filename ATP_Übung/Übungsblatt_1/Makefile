all: build/Übungsblatt_1.pdf

build/Übungsblatt_1.pdf: FORCE Übungsblatt_1.tex Übungsblatt_1.py |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	python Übungsblatt_1.py
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	Übungsblatt_1.tex
	

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean