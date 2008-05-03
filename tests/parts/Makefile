TOP         = ../..
V           = 0
CONFIG      = ../target
library     = qop-mdwf3

QMP_CFLAGS  = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)
QMP_LDFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --ldflags)
QMP_LIBS    = $(shell $(QMP_TOP:%=%/bin/qmp-config) --libs)

TCFLAGS = $(CFLAGS) $(QMP_CFLAGS) -I$(TOP)/$(LIMPDIR)

include $(TOP)/config/$(CONFIG)

tests = state.c \
        statex.c \
        conjgrad.c \
        conjgrad-2.c \
        conjgrad-f.c \
        conjgrad-fx.c \
        cg-op.c

x.tests =

#         heater.c \
#         localheat.c

opdriver = optest.c

optests = test-op-dummy.c \
          test-op-a.c \
          test-op-a1.c \
          test-op-ax.c \
          test-op-ax1.c \
          test-op-b.c \
          test-op-b1.c \
          test-op-bx.c \
          test-op-bx1.c \
          test-op-f0.c \
          test-op-f.c

opytests = test-op-a-2.c \
           test-op-b-2.c \
           test-op-f-2.c \
           test-op-ddw-2.c \

#           test-op-ddwq-2.c

opxdriver = opxtest.c

opxtests = test-op-axa.c \
           test-op-a1a.c \
           test-op-ax1ax.c \
           test-op-bxb.c \
           test-op-b1b.c \
           test-op-bx1bx.c \
           test-op-fxf.c

opztests = test-op-axa-2.c \
           test-op-ax1ax-2.c \
           test-op-a1a-2.c \
           test-op-bxb-2.c \
           test-op-bx1bx-2.c \
           test-op-b1b-2.c \
           test-op-fxf-2.c \
           test-op-af-2.c \
           test-op-ba-2.c \
           test-op-baf-2.c \
           test-op-abf-2.c \
           test-op-1f-2.c \
           test-op-1fx-2.c \
           test-op-1baf-2.c \
           test-op-1bafn-2.c \
           test-op-ab-2.c \
           test-op-md-2.c \
           test-op-mdx-2.c \
           test-op-mdxmd-2.c \
           test-op-axa-3.c \
           test-op-ax1ax-3.c \
           test-op-a1a-3.c \
           test-op-bxb-3.c \
           test-op-bx1bx-3.c \
           test-op-b1b-3.c \
           test-op-fxf-3.c \
           test-op-af-3.c \
           test-op-ba-3.c \
           test-op-baf-3.c \
           test-op-abf-3.c \
           test-op-1f-3.c \
           test-op-1fx-3.c \
           test-op-1baf-3.c \
           test-op-1bafn-3.c \
           test-op-m-3.c \
           test-op-mn-3.c \
           test-op-mx-3.c \
           test-op-mxn-3.c \

oproutines = op-a.c \
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
             op-project-n.c \
             op-project-x.c \
             op-boundary.c \
             op-apf.c \
             op-1mf.c \
             op-1mfx.c \
             op-1mba1f.c \
             op-1mba1fn.c \
             op-a1xbxfx.c \
             op-ba1f.c \
	     op-ba1.c \
             dot-fermion.c

oplib = op

libname     = $(oplib:%=lib%.a)
objects     = $(tests:%.c=%.o)
x.objects   = $(x.tests:%.c=%-d.o) $(x.tests:%.c=%-f.o)
op.objects  = $(optests:%.c=%.o) $(opdriver:%.c=%.o)
opx.objects = $(opxtests:%.c=%.o) $(opxdriver:%.c=%.o)
opy.objects = $(opytests:%.c=%.o)
opz.objects = $(opztests:%.c=%.o)
opl.objects = $(oproutines:%.c=%.o)
all.objects = $(objects) \
              $(x.objects) \
              $(op.objects) \
              $(opx.objects) \
              $(opy.objects) \
              $(opz.objects) \
              $(opl.objects)

programs = $(objects:%.o=%) \
           $(x.objects:%.o=%) \
           $(optests:%.c=%) \
           $(opytests:%.c=%) \
           $(opxtests:%.c=%) \
           $(opztests:%.c=%) \

.PHONY: all clean realclean dist optests opxtests opytests opztests

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

all: $(programs) optests opxtests opytests opztests

optests: $(optests:%.c=%)

opxtests: $(opxtests:%.c=%)

opytests: $(opxtests:%.c=%)

opztests: $(opztests:%.c=%)

clean:
	$E RM tests/objects
	$C$(RM) $(all.objects)

dist realclean: clean
	$E RM tests/programs
	$C$(RM) $(programs)
	$E RM tests/libraries
	$C$(RM) $(libname)

$(objects:%.o=%) $(x.objects:%.o=%): %: %.o
	$E LD $@
	$C$(CC) $(TCFLAGS) $(QMP_LDFLAGS) -L$(TOP)/$(LIMPDIR) -o $@ $< \
              -l$(library) $(QMP_LIBS)

$(objects): %.o: %.c
	$E CC $@
	$C$(CC) $(TCFLAGS) -I$(TOP)/port -c -o $@ $<

$(x.tests:%.c=%-f.o): %-f.o: %.c
	$E CC-f $@
	$C$(CC) $(TCFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'F\' -c -o $@ $<

$(x.tests:%.c=%-d.o): %-d.o: %.c
	$E CC-d $@
	$C$(CC) $(TCFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

###
$(opdriver:%.c=%.o): %.o: %.c
	$E CC-d $@
	$C$(CC) $(TCFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

###
$(optests:%.c=%.o): %.o: %.c
	$E CC-d $@
	$C$(CC) $(TCFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(optests:%.c=%): %: %.o
	$E CC $@
	$C$(CC) $(TCFLAGS) $(QMP_LDFLAGS) -L$(TOP)/$(LIMPDIR) -o $@ $^ \
              -l$(library) $(QMP_LIBS)

$(optests:%.c=%): $(opdriver:%.c=%.o)

$(optests:%.c=%.o) $(opdriver:%.c=%.o): optest.h

###
$(opytests:%.c=%.o): %.o: %.c
	$E CC-d $@
	$C$(CC) $(TCFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(opytests:%.c=%): %: %.o
	$E CC $@
	$C$(CC) $(TCFLAGS) $(QMP_LDFLAGS) -L$(TOP)/$(LIMPDIR) -o $@ $^ \
              -L. -lop -l$(library) $(QMP_LIBS)

$(opytests:%.c=%): $(opdriver:%.c=%.o)

$(opytests:%.c=%.o): op-routines.h optest.h \
                     $(TOP)/port/mdwf.h $(TOP)/$(LIMPDIR)/qop-mdwf3.h

$(opytests:%.c=%): $(libname)

###
$(opxtests:%.c=%.o) $(opxdriver:%.c=%.o): %.o: %.c
	$E CC-d $@
	$C$(CC) $(TCFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(opxtests:%.c=%): %: %.o
	$E CC $@
	$C$(CC) $(TCFLAGS) $(QMP_LDFLAGS) -L$(TOP)/$(LIMPDIR) -o $@ $^ \
              -l$(library) $(QMP_LIBS)

$(opxtests:%.c=%): $(opxdriver:%.c=%.o)

$(opxtests:%.c=%) $(optests:%.c=%): $(TOP)/$(LIMPDIR)/lib$(library).a

$(opxtests:%.c=%.o) $(opxdriver:%.c=%.o): opxtest.h
###

$(opztests:%.c=%.o): %.o: %.c
	$E CC-d $@
	$C$(CC) $(TCFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(opztests:%.c=%): %: %.o
	$E CC $@
	$C$(CC) $(TCFLAGS) $(QMP_LDFLAGS) -L$(TOP)/$(LIMPDIR) -o $@ $^ \
              -l$(library) $(QMP_LIBS)

$(opztests:%.c=%): $(opxdriver:%.c=%.o)

$(opztests:%.c=%): $(TOP)/$(LIMPDIR)/lib$(library).a $(libname)

$(opztests:%.c=%.o): opxtest.h op-routines.h
###

$(libname): $(opl.objects)
	$E AR $(libname)
	$C$(AR) cr $(libname) $(opl.objects)
	$C$(RANLIB) $(libname)

$(opl.objects): %.o: %.c
	$E CC $@
	$C$(CC) $(TCFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(objl.objects): $(TOP)/port/mdwf.h $(TOP)/$(LIMPDIR)/qop-mdwf3.h op-routines.h

###
$(all.objects): $(TOP)/port/mdwf.h $(TOP)/$(LIMPDIR)/qop-mdwf3.h

$(programs): $(TOP)/$(LIMPDIR)/lib$(library).a

conjgrad-2.o: clear-mdwf.h

conjgrad.o: clear-mdwf.h do-cg.c

conjgrad-f.o: clear-mdwf.h

conjgrad-fx.o: clear-mdwf.h
