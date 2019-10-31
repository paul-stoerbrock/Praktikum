all: build/V803.pdf

build/plot.pdf: V803.py data.txt  |build
	python V803.py

build/V803.pdf: V803.tex build/plot.pdf |build
	lualatex --outputdirectory=build V803.tex

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
