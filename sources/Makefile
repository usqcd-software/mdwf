NOWEAVE  = noweave -filter btdefn -index -latex -delay
NOTANGLE = notangle -filter btdefn
NOINDEX  = noindex
MZSCHEME = mzscheme
MZC      = mzc
INDENT   = indent
LATEX    = latex
DVIPDFM  = dvipdfm
TFLAGS   = -L'\#line %L "%F"%N'
CC       = gcc
INCDIRS  = $$HOME/QDP/qmp-2.1.7/include
CFLAGS   = -Wall -I. $(INCDIRS:%=-I%)

noweb = mdwf.nw qa0.nw

h.sources = qop-mdwf3.h \
            mdwf.h

q.sources = error.c \
            finit.c \
            performance.c \
            memory.c \
            state-dump.c

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
             op-ApF.c

ss.sources = basis.ss

qs.sources = common.ss \
             parser.ss \
             cfolding.ss \
             q2complex.ss \
             qa0print.ss \
             mk-cheader.ss \
             backends.ss \
             qa0.xx

qz.objects = $(q.sources:%.c=%.o) \
             $(qx.sources:%.c=%f.o) \
             $(qx.sources:%.c=%d.o)

c.sources = $(q.sources) \
            $(qx.sources) \
            test-state.c

x.sources = $(h.sources) $(c.sources)

sources = $(x.sources) $(ss.sources) $(qs.sources)

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
	$(RM) $(f.pdf) $(sources) test-state qa0.ss qa0

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
	($(NOTANGLE) $(TFLAGS) -R'File [[$@]]' -c $< >$@) || ($(RM) $@; exit 1)

$(ss.sources): mdwf.nw
	($(NOTANGLE) -R'File [[$@]]' -c $< >$@) || ($(RM) $@; exit 1)

$(qs.sources): qa0.nw
	($(NOTANGLE) -R'File [[$@]]' -c $< >$@) || ($(RM) $@; exit 1)

test-state: $(test.state.sources)
	$(CC) -g -Wall -DDEBUG_INIT -I. $(INCDIRS:%=-I%) \
              -o $@ $(test.state.sources)

test-state: $(h.sources)

$(f.lib): $(qz.objects)
	$(AR) cr $(f.lib) $(qz.objects)

$(q.sources:%.c=%.o): %.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(qx.sources:%.c=%d.o): %d.o: %.c
	$(CC) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' $(CFLAGS) -c -o $@ $<

$(qx.sources:%.c=%f.o): %f.o: %.c
	$(CC) -DQOP_MDWF_DEFAULT_PRECISION=\'F\' $(CFLAGS) -c -o $@ $<

$(qz.objects): $(h.sources)

qa0.ss: $(qs.sources) mk-qa0.ss
	echo '(load "mk-qa0.ss")' | $(MZSCHEME)

qa0: qa0.ss
	$(MZC) --exe $@ $<
