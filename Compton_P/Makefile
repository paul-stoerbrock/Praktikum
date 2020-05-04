all: build/Compton_P.pdf

build/Compton_P.pdf: FORCE compton_p.tex compton_p.py build/plotAl.pdf build/plotCu.pdf |build
	TEXINPUTS="$(call translate,build:)" \
	BIBINPUTS=build:\
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=batchmode \
	  --halt-on-error \
	compton_p.tex

#build/plotCu.pdf: compton_p.py EmissionCu.dat |build
#
#	TEXINPUTS=$$(pwd): python compton_p.py

build/plotAl.pdf:  compton_p.py ComptonAl.txt |build

	TEXINPUTS=$$(pwd): python compton_p.py

build/plotCu.pdf:  compton_p.py EmissionCu.dat |build

	TEXINPUTS=$$(pwd): python compton_p.py


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean