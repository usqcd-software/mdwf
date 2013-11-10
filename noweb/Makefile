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
        qcd-defs.nw

tex = eigen-cg.tex

incs = defs.tex

sources = ../port/qop-mdwf3.h

f.tex = $(noweb:%.nw=%.tex)
f.dvi = $(noweb:%.nw=%.dvi) $(tex:%.tex=%.dvi)
f.pdf = $(noweb:%.nw=../doc/%.pdf) $(tex:%.tex=../doc/%.pdf)

.PHONY: all clean realclean 

all: pdf sources

clean:
	$(RM) $(noweb:%.nw=%.aux) $(noweb:%.nw=%.dvi) $(noweb:%.nw=%.log) \
              $(noweb:%.nw=%.nwi) $(noweb:%.nw=%.out) $(noweb:%.nw=%.tex) \
              $(noweb:%.nw=%.toc)
	$(RM) $(tex:%.tex=%.aux) $(tex:%.tex=%.dvi) $(tex:%.tex=%.log) \
              $(tex:%.tex=%.nwi) $(tex:%.tex=%.out) $(tex:%.tex=%.toc) \
              $(tex:%.tex=%.loa)

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

$(tex:%.tex=%.dvi): $(incs)

$(sources): mdwf.nw
	$(NOTANGLE) $(TFLAGS) -R'File [[$@]]' -c $< | $(CPIF) $@
