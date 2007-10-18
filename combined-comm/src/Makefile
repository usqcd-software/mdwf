# This is a working copy (MacOS)

CONFIG     = macosx
TOP        = .
QMP_CFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)

include config/$(CONFIG)

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
            set-generic.c \
            step-even.c \
            step-odd.c \
            version.c \
            init.c \
            fini.c \
            comm-setup.c \
            comm-free.c

p.sources = f-alloc.c \
            f-import.c \
            f-export.c \
            f-free.c \
            h-alloc.c \
            h-import.c \
            h-export.c \
            h-free.c \
            g-free.c \
            g-import.c \
            x-import.c \
            x-export.c \
            f-madd.c \
            h-madd.c \
            f-dot.c \
            h-dot.c \
            f-norm2.c \
            h-norm2.c  \
            op-dirac.c

x.objects = $(x.sources:%.c=%.o)

p.objects = $(p.sources:%.c=%f.o) \
            $(p.sources:%.c=%d.o)

objects = $(x.objects) $(p.objects)

headers = mdwf.h \
          qop-mdwf3.h

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
	$C$(MAKE) -C $(LIMP) 'V=$V' 'CONFIG=$(CONFIG)' 'library=$(library)' $@
	$C$(MAKE) -C tests 'V=$V' 'CONFIG=$(CONFIG)' $@

clean:
	$E RM objects
	$C$(RM) $(objects)
	$C$(MAKE) -C $(LIMP) 'V=$V' 'CONFIG=$(CONFIG)' 'library=$(library)' $@
	$C$(MAKE) -C tests 'V=$V' 'CONFIG=$(CONFIG)' $@

realclean: clean
	$E RM library
	$C$(RM) $(library)
	$C$(MAKE) -C $(LIMP) 'V=$V' 'CONFIG=$(CONFIG)' 'library=$(library)' $@
	$C$(MAKE) -C tests 'V=$V' 'CONFIG=$(CONFIG)' $@

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