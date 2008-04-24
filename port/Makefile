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

p.sources = fermion-alloc.c \
            fermion-import.c \
            fermion-export.c \
            fermion-free.c \
            fermion-madd.c \
            fermion-dot.c \
            fermion-norm2.c \
            gauge-free.c \
            gauge-import.c \
            hfermion-alloc.c \
            hfermion-import.c \
            hfermion-export.c \
            hfermion-free.c \
            hfermion-madd.c \
            hfermion-dot.c \
            hfermion-norm2.c  \
            x-import.c \
            x-export.c \
            boundary-n.c \
            boundary-x.c \
            op-boundary.c \
            op-a.c \
            op-ax.c \
	    op-a1.c \
            op-a1x.c \
            op-b.c \
            op-bx.c \
	    op-b1.c \
            op-b1x.c \
            op-f.c \
            op-fx.c \
            op-a1xbx.c \
            op-apf.c \
            op-1mf.c \
            op-1mfx.c \
            op-1mfxn.c \
            op-1mba1f.c \
            op-1mba1fn.c \
            op-a1xbxfx.c \
            op-axpbxfx.c \
            op-ba1f.c \
            op-ba1.c \
            op-d.c \
            op-m.c \
            op-mn.c \
            op-mx.c \
            op-mxn.c \
            op-dirac.c \
            op-dirac-conj.c \
            op-precond.c \
            op-precond-conj.c \
#            cg-dirac.c \


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

$(lib): library

library: $(objects)
	$E AR $(lib)
	$C$(AR) cr $(lib) $^
	$C$(RANLIB) $(lib)

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
