all: build/Übungsblatt_2.pdf

build/Übungsblatt_2.pdf: FORCE Übungsblatt_2.tex |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	Übungsblatt_2.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean