NOWEAVE  = noweave -filter btdefn -index -latex -delay
NOTANGLE = notangle -filter btdefn $(TFLAGS)
NOINDEX  = noindex
INDENT   = indent
LATEX    = latex
DVIPDFM  = dvipdfm
TFLAGS   = -L\#'line %L "%F"%N'
CC       = gcc
INCDIRS  = $$HOME/QDP/qmp-2.1.7/include
CFLAGS   = -Wall -I. $(INCDIRS:%=-I%)

noweb = mdwf.nw

h.sources = qop-mdwf3.h \
            mdwf.h

q.sources = error.c \
            finit.c \
            performance.c \
            memory.c \
            state-dump.c

qx.sources = memory-g.c \
             memory-h.c \
             memory-f.c

qz.objects = $(q.sources:%.c=%.o) \
             $(qx.sources:%.c=%f.o) \
             $(qx.sources:%.c=%d.o)

c.sources = $(q.sources) \
            $(qx.sources) \
            test-state.c

sources = $(h.sources) $(c.sources)

f.lib = libqtest.a

test.state.sources = test-state.c state-dump.o finit.o memory.o error.o

f.tex = $(noweb:%.nw=%.tex)
f.dvi = $(noweb:%.nw=%.dvi)
f.pdf = $(noweb:%.nw=%.pdf)


.PHONY: all clean realclean sources pdf test-lib

all: pdf sources

clean:
	$(RM) $(noweb:%.nw=%.aux) $(noweb:%.nw=%.dvi) $(noweb:%.nw=%.log) \
              $(noweb:%.nw=%.nwi) $(noweb:%.nw=%.out) $(noweb:%.nw=%.tex) \
              $(noweb:%.nw=%.toc) $(qz.objects) $(f.lib)

realclean: clean
	$(RM) $(f.pdf) $(sources) test-state

sources: $(sources)

pdf: $(f.pdf)

test-lib: $(f.lib)

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

$(f.lib): $(qz.objects)
	$(AR) cr $(f.lib) $(qz.objects)

$(q.sources:%.c=%.o): %.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(qx.sources:%.c=%d.o): %d.o: %.c
	$(CC) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' $(CFLAGS) -c -o $@ $<

$(qx.sources:%.c=%f.o): %f.o: %.c
	$(CC) -DQOP_MDWF_DEFAULT_PRECISION=\'F\' $(CFLAGS) -c -o $@ $<

$(qz.objects): $(h.sources)