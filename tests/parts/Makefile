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

optests = op-dummy.c \
          op-ddw.c \
          op-a.c \
          op-b.c \
          op-f.c \
          op-ax.c \
          op-bx.c \
          op-f0.c

opdriver = optest.c

objects = $(tests:%.c=%.o)
x.objects = $(x.tests:%.c=%-d.o) $(x.tests:%.c=%-f.o)
op.objects = $(optests:%.c=%.o) $(opdriver:%.c = %.o)

programs = $(objects:%.o=%) $(x.objects:%.o=%)

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

clean:
	$E RM tests/objects
	$C$(RM) $(objects) $(x.objects)

dist realclean: clean
	$E RM tests/programs
	$C$(RM) $(programs)

$(programs): %: %.o
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

$(x.objects) $(objects): $(TOP)/port/mdwf.h $(TOP)/$(LIMPDIR)/qop-mdwf3.h

conjgrad.o: clear-mdwf.h do-cg.c

$(optests:%.c=%) $(programs): $(TOP)/$(LIMPDIR)/lib$(library).a

$(optests:%.c=%.o) $(opdriver:%.c=%.o): optest.h
