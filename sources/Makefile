NOWEAVE  = noweave -filter btdefn -index -latex -delay
NOTANGLE = notangle -filter btdefn
NOINDEX  = noindex
CPIF     = cpif
MZSCHEME = mzscheme
MZC      = mzc
INDENT   = indent
LATEX    = latex
DVIPDFM  = dvipdfm
TFLAGS   = -L'\#line %L "%F"%N'
CC       = gcc
AR       = ar
QA0c     = ./qa0 c99
QA0h     = ./qa0 header
INCDIRS  = $$HOME/QDP/qmp-2.1.7/include ./c99
CFLAGS   = -Wall -I. $(INCDIRS:%=-I%)

noweb = mdwf.nw qa0.nw

h.sources = qop-mdwf3.h \
            mdwf.h

q.sources = error.c \
            finit.c \
            performance.c \
            memory.c \
            state-dump.c \
            comm.c \
            params.c

qx.sources = memory-g.c \
             memory-f.c \
             memory-h.c \
             import-g.c \
             import-f.c \
             import-h.c \
             export-f.c \
             export-h.c \
             operator-ddw.c \
             operator-ddwc.c \
             op-ApF.c \
             op-A.c \
             op-B.c

qi.sources = spin.c \
             ABF.c \
             simport-f.c

ss.sources = basis.ss

qs.sources = common.ss \
             parser.ss \
             cfolding.ss \
             q2complex.ss \
             qa0print.ss \
             mk-cheader.ss \
             backends.ss \
             qa0.xx

qa.sources = c99/spin.qa0 \
             c99/ABF.qa0 \
             c99/simport-f.qa0

qz.objects = $(q.sources:%.c=%.o) \
             $(qx.sources:%.c=%f.o) \
             $(qx.sources:%.c=%d.o) \
             $(qi.sources:%.c=%d.o) \
             $(qi.sources:%.c=%f.o) \

c.sources = $(q.sources) \
            $(qx.sources) \
            test-state.c

x.sources = $(h.sources) $(c.sources) $(qi.sources)

sources = $(x.sources) $(ss.sources) $(qs.sources) $(qa.sources)

f.lib = libqtest.a

test.state.sources = test-state.c state-dump.c finit.c memory.c error.c

f.tex = $(noweb:%.nw=%.tex)
f.dvi = $(noweb:%.nw=%.dvi)
f.pdf = $(noweb:%.nw=%.pdf)


.PHONY: all clean realclean sources pdf test-lib test-symbols

all: pdf sources

clean:
	$(RM) $(noweb:%.nw=%.aux) $(noweb:%.nw=%.dvi) $(noweb:%.nw=%.log) \
              $(noweb:%.nw=%.nwi) $(noweb:%.nw=%.out) $(noweb:%.nw=%.tex) \
              $(noweb:%.nw=%.toc) $(qz.objects) $(f.lib) qa0.xx qa0.ss

realclean: clean
	$(RM) $(f.pdf) $(sources) test-state qa0.ss qa0 c99/*

sources: $(sources)

pdf: $(f.pdf)

test-lib: $(f.lib)

test-symbols: test-lib
	@nm -fP $(f.lib) | \
        awk '/[^:]$$/{if (sym[$$1]!="T") sym[$$1]=$$2} \
            END{for (s in sym) if (sym[s]=="U") printf "%s %s\n", sym[s], s;}'|\
        sort

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

$(x.sources): mdwf.nw
	$(NOTANGLE) $(TFLAGS) -R'File [[$@]]' -c $< | $(CPIF) $@

$(qa.sources) $(ss.sources): mdwf.nw
	$(NOTANGLE) -R'File [[$@]]' -c $< | $(CPIF) $@

$(qs.sources): qa0.nw
	$(NOTANGLE) -R'File [[$@]]' -c $< | $(CPIF) $@

test-state: $(test.state.sources)
	$(CC) -g -Wall -DDEBUG_INIT -I. $(INCDIRS:%=-I%) \
              -o $@ $(test.state.sources)

test-state: $(h.sources)

$(f.lib): $(qz.objects)
	$(AR) cr $(f.lib) $(qz.objects)

$(q.sources:%.c=%.o): %.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(qi.sources:%.c=%d.o) $(qx.sources:%.c=%d.o): %d.o: %.c
	$(CC) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' $(CFLAGS) -c -o $@ $<

$(qi.sources:%.c=%f.o) $(qx.sources:%.c=%f.o): %f.o: %.c
	$(CC) -DQOP_MDWF_DEFAULT_PRECISION=\'F\' $(CFLAGS) -c -o $@ $<

$(qi.sources:%.c=%f.o): $(qa.sources:%.qa0=%f.c)
$(qi.sources:%.c=%d.o): $(qa.sources:%.qa0=%d.c)

$(qz.objects): $(h.sources)
$(qz.objects): $(qa.sources:%.qa0=%f.h)
$(qz.objects): $(qa.sources:%.qa0=%d.h)

qa0.ss: $(qs.sources) basis.ss mk-qa0.ss
	echo '(load "mk-qa0.ss")' | $(MZSCHEME)

qa0: qa0.ss
	$(MZC) --exe $@ $<

$(qa.sources:%.qa0=%d.c): %d.c: %.qa0
	$(QA0c) double $@ $<

$(qa.sources:%.qa0=%f.c): %f.c: %.qa0
	$(QA0c) float $@ $<

$(qa.sources:%.qa0=%d.h): %d.h: %.qa0
	$(QA0h) double $@ $<

$(qa.sources:%.qa0=%f.h): %f.h: %.qa0
	$(QA0h) float $@ $<


$(qa.sources:%.qa0=%d.c): qa0
$(qa.sources:%.qa0=%f.c): qa0
$(qa.sources:%.qa0=%d.h): qa0
$(qa.sources:%.qa0=%f.h): qa0

