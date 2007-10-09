# This is a working copy (MacOS)
QMP_TOP    = $$HOME/QDP/qmp-2.1.7
QMP_CFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)
CC         = gcc
CFLAGS     = -Wall -g -I. $(QMP_CFLAGS)
AR         = ar
RANLIB     = ranlib

x.sources = error.c \
            alloc-aligned.c \
            alloc-eo.c \
            alloc.c \
            dealloc.c \
            perf.c \
            set-borichi.c \
            set-chiu.c \
            set-moebius.c \
            set-shamir.c \
            step-even.c \
            step-odd.c \
            version.c

p.sources =

x.objects = $(x.sources:%.c=%.o)

p.objects = $(p.sources:%.c=%f.o) \
            $(p.sources:%.c=%d.o)

objects = $(x.objects) $(p.objects)

headers = mdwf.h \
          qop-mdwf3.h

library = libqop-mdwf3.a

.PHONY: all clean realclean

all: $(library)

clean:
	$(RM) $(objects)

realclean: clean
	$(RM) $(library)


$(library): $(objects)
	$(AR) cr $@ $<
	$(RANLIB) $@

$(x.objects): %.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(p.objects): %f.o: %.c
	$(CC) $(CFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'F\' -c -o $@ $<

$(p.objects): %d.o: %.c
	$(CC) $(CFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(objects): $(headers)
