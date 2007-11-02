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

objects = $(tests:%.c=%.o)

programs = $(tests:%.c=%)

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
	$C$(RM) $(objects)

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

$(objects): ../mdwf.h ../qop-mdwf3.h

$(programs): ../lib$(library).a