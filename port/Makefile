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
            alloc.c \
            dealloc.c \
            perf.c \
            gauge-d2f.c \
            set-borichi.c \
            set-chiu.c \
            set-moebius.c \
            set-shamir.c \
            set-generic.c \
            p-free.c \
            version.c \
            init.c \
            fini.c \
            comm-setup.c \
            comm-free.c \
            deflator-la.c \
            deflator-alloc.c \
            deflator-free.c \
            deflator-reset.c \
            deflator-signature.c \
            deflator-is-stopped.c \
            deflator-resume.c \
            deflator-stop.c \
            deflator-eigen.c \
            deflator-preamble.c \
            deflator-postamble.c \
            deflator-update0.c \
            deflator-update1.c \
            deflator-dim.c \
            deflator-extract-vector.c \
            deflator-add-vector.c \
            deflator-start-load.c \
            deflator-stop-load.c \
            deflator-end-load.c \
	    poly-normalize.c \
            hfermion-size.c \
            cg-mixed.c \
            solver-dirac-mixed.c \
            solver-dirac-deflated-mixed.c

p.sources = alloc-eo.c \
            step-even.c \
            step-odd.c \
            fermion-alloc.c \
            fermion-import.c \
            fermion-import4d.c \
            fermion-midpoint.c \
            fermion-export.c \
            fermion-export4d.c \
            fermion-free.c \
            fermion-madd.c \
            fermion-dot.c \
            fermion-norm2.c \
            gauge-free.c \
            gauge-import.c \
            hfermion-alloc.c \
            hfermion-import.c \
            hfermion-import-blas.c \
            hfermion-export.c \
            hfermion-export-blas.c \
            hfermion-free.c \
            hfermion-madd.c \
            hfermion-dot.c \
            hfermion-norm2.c  \
            vfermion-alloc.c \
            vfermion-free.c \
            vfermion-get.c \
            vfermion-put.c \
            x-import.c \
            x4-import.c \
            x-midpoint.c \
            x-export.c \
            x4-export.c \
            boundary-n.c \
            boundary-x.c \
            op-boundary.c \
            op-down-boundary.c \
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
            op-apfn.c \
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
            op-dn.c \
            op-m.c \
            op-mn.c \
            op-mx.c \
            op-mxn.c \
            axial-current.c \
            op-axial-current.c \
            op-dirac.c \
            op-dirac-conj.c \
            op-precond.c \
            op-precond-conj.c \
            op-mxm-poly.c \
            op-print.c \
	    cg-precondition.c \
            cg-solver.c \
            cg-inflate.c \
            cg-log.c \
            cg-dirac-error.c \
            scg-solver.c \
            solver-dirac.c \
            solver-mxm.c \
            solver-dxd.c \
            solver-shifted.c


x.objects = $(x.sources:%.c=$(LIMPDIR)/%.o)

p.objects = $(p.sources:%.c=$(LIMPDIR)/%f.o) \
            $(p.sources:%.c=$(LIMPDIR)/%d.o)

objects = $(x.objects) $(p.objects)

headers = mdwf.h \
          deflator-la.h \
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
	$C$(CC) $(CFLAGS) $(COPTS) -I. -c -o $@ $<

$(p.sources:%.c=$(LIMPDIR)/%f.o): $(LIMPDIR)/%f.o: %.c
	$E CC-f $<
	$C$(CC) $(CFLAGS) $(COPTS) -I. -DQOP_MDWF_DEFAULT_PRECISION=\'F\' -c -o $@ $<

$(p.sources:%.c=$(LIMPDIR)/%d.o): $(LIMPDIR)/%d.o: %.c
	$E CC-d $<
	$C$(CC) $(CFLAGS) $(COPTS) -I. -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(objects): $(headers)
