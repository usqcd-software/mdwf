TOP         = ..
V           = 0
CONFIG      = macosx
library     = qop-mdwf3

QMP_CFLAGS  = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)
QMP_LDFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --ldflags)
QMP_LIBS    = $(shell $(QMP_TOP:%=%/bin/qmp-config) --libs)

include $(TOP)/config/$(CONFIG)

tests = operator.c \
        state.c \
        statex.c

x.tests = heater.c

objects = $(tests:%.c=%.o)
x.objects = $(x.tests:%.c=%-d.o) $(x.tests:%.c=%-f.o)

programs = $(objects:%.o=%) $(x.objects:%.o=%)

.PHONY: all clean realclean dist

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

all: $(programs)

clean:
	$E RM tests/objects
	$C$(RM) $(objects) $(x.objects)

dist realclean: clean
	$E RM tests/programs
	$C$(RM) $(programs)

$(programs): %: %.o
	$E LD $@
	$C$(CC) $(CFLAGS) $(QMP_LDFLAGS) -L$(TOP) -o $@ $< \
              -l$(library) $(QMP_LIBS)

$(objects): %.o: %.c
	$E CC $<
	$C$(CC) $(CFLAGS) $(QMP_CFLAGS) -c -o $@ $<

$(x.tests:%.c=%-f.o): %-f.o: %.c
	$E CC-f $<
	$C$(CC) $(CFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'F\' -c -o $@ $<

$(x.tests:%.c=%-d.o): %-d.o: %.c
	$E CC-f $<
	$C$(CC) $(CFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<


$(x.objects) $(objects): ../mdwf.h ../qop-mdwf3.h

$(programs): ../lib$(library).a
