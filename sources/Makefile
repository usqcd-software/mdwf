INCDIRS  = $$HOME/QDP/qmp-2.1.7/include
NOWEAVE  = noweave -filter btdefn -index -latex -delay
NOTANGLE = notangle -filter btdefn $(TFLAGS)
NOINDEX  = noindex
INDENT   = indent
LATEX    = latex
DVIPDFM  = dvipdfm
TFLAGS   = -L\#'line %L "%F"%N'

noweb = mdwf.nw

h.sources = qop-mdwf3.h \
            mdwf.h

c.sources = error.c \
            finit.c \
            performance.c \
            memory.c \
            state-dump.c \
            test-state.c

sources = $(h.sources) $(c.sources)

test.state.sources = test-state.c state-dump.c finit.c memory.c error.c

f.tex = $(noweb:%.nw=%.tex)
f.dvi = $(noweb:%.nw=%.dvi)
f.pdf = $(noweb:%.nw=%.pdf)


.PHONY: all clean realclean sources pdf

all: pdf sources

clean:
	$(RM) $(noweb:%.nw=%.aux) $(noweb:%.nw=%.dvi) $(noweb:%.nw=%.log) \
              $(noweb:%.nw=%.nwi) $(noweb:%.nw=%.out) $(noweb:%.nw=%.tex) \
              $(noweb:%.nw=%.toc) *.o

realclean: clean
	$(RM) $(f.pdf) $(sources) test-state

sources: $(sources)

pdf: $(f.pdf)

###
$(f.pdf): %.pdf: %.dvi
	$(DVIPDFM) $<

$(f.dvi): %.dvi: %.tex
	$(LATEX) $<
	$(NOINDEX) $(<:%.tex=%)
	$(LATEX) $<
	$(NOINDEX) $(<:%.tex=%)
	$(LATEX) $<

$(f.tex): %.tex: %.nw
	$(NOWEAVE) $< > $@

$(sources): $(noweb)
	($(NOTANGLE) -R'File [[$@]]' -c $< > $@) || ($(RM) $@; exit 1)


test-state: $(test.state.sources)
	$(CC) -g -Wall -DDEBUG_INIT -I. $(INCDIRS:%=-I%) \
              -o $@ $(test.state.sources)

test-state: $(h.sources)