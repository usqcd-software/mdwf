TOP         = ../..
LIMPDIR     = c99
V           = 0
CONFIG      = macosx
library     = qop-mdwf3

QMP_CFLAGS  = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)
QMP_LDFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --ldflags)
QMP_LIBS    = $(shell $(QMP_TOP:%=%/bin/qmp-config) --libs)

TCFLAGS = $(CFLAGS) $(QMP_CFLAGS) -I$(TOP)/$(LIMPDIR)

include $(TOP)/config/$(CONFIG)

tests = operator.c \
        conjgrad.c \
        state.c \
        statex.c

x.tests = heater.c \
          localheat.c

optests = test-op-dummy.c \
          test-op-ddw.c \
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

opxtests = test-op-axa.c \
           test-op-a1a.c \
           test-op-ax1ax.c \
           test-op-bxb.c \
           test-op-b1b.c \
           test-op-bx1bx.c \
           test-op-fxf.c

opdriver = optest.c

opxdriver = opxtest.c

objects = $(tests:%.c=%.o)
x.objects = $(x.tests:%.c=%-d.o) $(x.tests:%.c=%-f.o)
op.objects = $(optests:%.c=%.o) $(opdriver:%.c=%.o)
opx.objects = $(opxtests:%.c=%.o) $(opxdriver:%.c=%.o)
all.objects = $(objects) $(x.objects) $(op.objects) $(opx.objects)

programs = $(objects:%.o=%) \
           $(x.objects:%.o=%) \
           $(optests:%.c=%) \
           $(opxtests:%.c=%)

.PHONY: all clean realclean dist optests

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

all: $(programs) optests

optests: $(optests:%.c=%)

opxtests: $(opxtests:%.c=%)

clean:
	$E RM tests/objects
	$C$(RM) $(all.objects)

dist realclean: clean
	$E RM tests/programs
	$C$(RM) $(programs)

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

$(optests:%.c=%.o) $(opdriver:%.c=%.o): %.o: %.c
	$E CC-d $@
	$C$(CC) $(TCFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(optests:%.c=%): %: %.o
	$E CC $@
	$C$(CC) $(TCFLAGS) $(QMP_LDFLAGS) -L$(TOP)/$(LIMPDIR) -o $@ $^ \
              -l$(library) $(QMP_LIBS)

$(optests:%.c=%): $(opdriver:%.c=%.o)

$(opxtests:%.c=%.o) $(opxdriver:%.c=%.o): %.o: %.c
	$E CC-d $@
	$C$(CC) $(TCFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(opxtests:%.c=%): %: %.o
	$E CC $@
	$C$(CC) $(TCFLAGS) $(QMP_LDFLAGS) -L$(TOP)/$(LIMPDIR) -o $@ $^ \
              -l$(library) $(QMP_LIBS)

$(opxtests:%.c=%): $(opxdriver:%.c=%.o)

$(all.objects): $(TOP)/port/mdwf.h $(TOP)/$(LIMPDIR)/qop-mdwf3.h

conjgrad.o: clear-mdwf.h do-cg.c

$(opxtests:%.c=%) $(optests:%.c=%) $(programs): $(TOP)/$(LIMPDIR)/lib$(library).a

$(optests:%.c=%.o) $(opdriver:%.c=%.o): optest.h

$(opxtests:%.c=%.o) $(opxdriver:%.c=%.o): opxtest.h

