# This is a working copy (MacOS)
QMP_TOP    = $$HOME/QDP/qmp-2.1.7
QMP_CFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)
LIMP       = faux
CC         = gcc
CFLAGS     = -Wall -g -I. $(QMP_CFLAGS) -I$(LIMP)
AR         = ar
RANLIB     = ranlib

V          = 0

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
            version.c \
            fini.c

p.sources = f-alloc.c \
            f-import.c \
            f-export.c \
            f-free.c \
            h-alloc.c \
            h-import.c \
            h-export.c \
            h-free.c \
            g-free.c \
            g-import.c

x.objects = $(x.sources:%.c=%.o)

p.objects = $(p.sources:%.c=%f.o) \
            $(p.sources:%.c=%d.o)

objects = $(x.objects) $(p.objects)

headers = mdwf.h \
          qop-mdwf3.h

limp.headers = 

library = libqop-mdwf3.a

.PHONY: all clean realclean

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

all: $(library)

clean:
	$E RM objects
	$C$(RM) $(objects)

realclean: clean
	$E RM library
	$C$(RM) $(library)


$(library): $(objects)
	$E AR $@
	$C$(AR) cr $@ $^
	$C$(RANLIB) $@

$(x.objects): %.o: %.c
	$E CC $<
	$C$(CC) $(CFLAGS) -c -o $@ $<

$(p.sources:%.c=%f.o): %f.o: %.c
	$E CC-f $<
	$C$(CC) $(CFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'F\' -c -o $@ $<

$(p.sources:%.c=%d.o): %d.o: %.c
	$E CC-d $<
	$C$(CC) $(CFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(objects): $(headers) $(limp.headers:%=$(LIMP)/%)