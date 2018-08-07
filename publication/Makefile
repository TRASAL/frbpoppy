
# * We use the modern lualatex engine which allows choosing opentype system fonts
# * All LaTeX output will go to the `buildplot` directory which is gitignored
# * We stop at the first error
# * We don't allow LaTeX to prompt us, Makefiles should be noninteractive
texopts = --lualatex \
	  --output-directory=buildplot \
	  --halt-on-error \
	  --interaction=nonstopmode


buildplot/rates.pdf: matplotlibrc matplotlib_header.tex rates.py | buildplot
	TEXINPUTS=$(shell pwd): python3 rates.py

buildplot:
	mkdir -p buildplot

clean:
	rm -rf buildplot/

FORCE:

.PHONY: clean all FORCE
