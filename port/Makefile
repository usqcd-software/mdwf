#
CONFIG     = macosx
LIBRARY    = libqop-mdwf3.a
LIMPDIR    = .

QMP_CFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)
lib        = $(LIMPDIR)/$(LIBRARY)
HDR        = qop-mdwf3.h

include ../config/$(CONFIG)

V          = 0

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

.PHONY: all clean realclean library dist

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
            p-free.c \
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
            cg-dirac.c \
            op-dirac.c \
#            op-dirac-conj.c

x.objects = $(x.sources:%.c=$(LIMPDIR)/%.o)

p.objects = $(p.sources:%.c=$(LIMPDIR)/%f.o) \
            $(p.sources:%.c=$(LIMPDIR)/%d.o)

objects = $(x.objects) $(p.objects)

headers = mdwf.h \
          $(HDR)


all: $(lib)
ifeq ("$(LIMPDIR)", ".")
	$E $(HDR) in place
else
	$E CP $(HDR)
	$C cp $(HDR) $(LIMPDIR)/
endif

clean:
	$E RM objects
	$C$(RM) $(objects)

realclean: clean
ifneq ("$(LIMPDIR)", ".")
	$C$(RM) $(LIMPDIR)/$(HDR)
endif


dist: clean

$(lib): $(objects)
	$E AR $@
	$C$(AR) cr $@ $^
	$C$(RANLIB) $@

$(x.objects): $(LIMPDIR)/%.o: %.c
	$E CC $<
	$C$(CC) $(CFLAGS) -I. -c -o $@ $<

$(p.sources:%.c=$(LIMPDIR)/%f.o): $(LIMPDIR)/%f.o: %.c
	$E CC-f $<
	$C$(CC) $(CFLAGS) -I. -DQOP_MDWF_DEFAULT_PRECISION=\'F\' -c -o $@ $<

$(p.sources:%.c=$(LIMPDIR)/%d.o): $(LIMPDIR)/%d.o: %.c
	$E CC-d $<
	$C$(CC) $(CFLAGS) -I. -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(objects): $(headers)