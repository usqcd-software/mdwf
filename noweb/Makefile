NOWEAVE  = noweave -filter btdefn -index -latex -delay
NOTANGLE = notangle -filter btdefn
NOINDEX  = noindex
CPIF     = cpif
MZSCHEME = mzscheme
MZC      = mzc
INDENT   = indent
LATEX    = latex
DVIPDFM  = dvipdfm
#TFLAGS   = -L'\#line %L "%F"%N'

noweb = mdwf.nw \
        qa0.nw

sources = ../src/qop-mdwf3.h \
          ../utils/basis.ss

f.tex = $(noweb:%.nw=%.tex)
f.dvi = $(noweb:%.nw=%.dvi)
f.pdf = $(noweb:%.nw=../doc/%.pdf)

.PHONY: all clean realclean 

all: pdf sources

clean:
	$(RM) $(noweb:%.nw=%.aux) $(noweb:%.nw=%.dvi) $(noweb:%.nw=%.log) \
              $(noweb:%.nw=%.nwi) $(noweb:%.nw=%.out) $(noweb:%.nw=%.tex) \
              $(noweb:%.nw=%.toc)

realclean: clean
	$(RM) $(f.pdf) $(sources)

sources: $(sources)

pdf: $(f.pdf)

###
$(f.pdf): ../doc/%.pdf: %.dvi
	$(DVIPDFM) -o $@ $<

$(f.dvi): %.dvi: %.tex
	$(LATEX) $<
	$(NOINDEX) $(<:%.tex=%)
	$(LATEX) $<
	$(NOINDEX) $(<:%.tex=%)
	$(LATEX) $<

$(f.tex): %.tex: %.nw
	$(NOWEAVE) $< > $@

$(sources): mdwf.nw
	$(NOTANGLE) $(TFLAGS) -R'File [[$@]]' -c $< | $(CPIF) $@
