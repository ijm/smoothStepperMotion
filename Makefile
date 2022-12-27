PYTHON = python3
PDFLATEX = pdflatex
LYX=lyx
MKINXJSON=../scripts/mkInxJSON.py

PROFILES=profile1.in profileWobble.in
SECTS=$(patsubst %.in,%.section.tex,$(PROFILES))
MOTORS=-m 17Y408S.dat 

NAME = smoothStepper
all: $(NAME).pdf $(NAME).inx.json

distclean:
	rm ./*.tex ./*.pdf ./*.pgf ./*.c ./*.aux $(NAME).log $(NAME).out
	rm -rf __pycache__

polynomials.tex polynomials.%.pgf: polynomials.dat
	$(PYTHON) plotPolys.py -p polynomials.dat

%.c %.section.tex: %.in
	$(PYTHON) generate.py -t profile -p polynomials.dat $(MOTORS) -i $< 

plots.tex: $(SECTS)
	-rm $@
	for f in $(SECTS); do echo "\\input{$$f}" >> $@; done
        
$(NAME).pdf: $(NAME).tex plots.tex $(SECTS) polynomials.tex
	$(PDFLATEX) $(NAME).tex && $(PDFLATEX) $(NAME).tex && $(PDFLATEX) $(NAME).tex 

%.tex: %.lyx
	$(LYX) --export pdflatex $<

%.inx.json: %.tex
	if [ -f "$(MKINXJSON)" ]; then $(PYTHON) $(MKINXJSON) -i $< -o $@ ; fi

