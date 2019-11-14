all: build/V802.pdf

build/V802.pdf: V802.tex |build
	TEXINPUTS="$(call translate,build:)" \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	V802.tex

#build/plot.pdf: V802.py data.txt  |build
	
#	python 802.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean