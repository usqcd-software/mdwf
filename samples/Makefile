TOP     = ..
CONFIG  = ../target
V       = 0

QMP_CFLAGS  = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)
QMP_LDFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --ldflags)
QMP_LIBS    = $(shell $(QMP_TOP:%=%/bin/qmp-config) --libs)

include ../config/$(CONFIG)

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

.PHONY: all clean realclean

sources = solver-ddw.c \
          solver-mxm.c \
          solver-shifted.c

xsources = solver-mixed.c

common = common.c

headers = common.h

f.objects = $(sources:%.c=%-f.o)
d.objects = $(sources:%.c=%-d.o)
x.objects = $(xsources:%.c=%.o) $(common:%.c=%.o)

f.programs = $(sources:%.c=%-f)
d.programs = $(sources:%.c=%-d)
x.programs = $(xsources:%.c=%)

programs = $(f.programs) $(d.programs) $(x.programs)

all: $(programs)

clean:
	$(RM) $(f.objects) $(d.objects) $(x.objects)

realclean: clean
	$(RM) $(programs)


$(programs): %: %.o
	$E LD $@
	$C$(CC) $(CFLAGS) $(QMP_LDFLAGS) -L$(TOP)/$(LIMPDIR) -o $@ $^ \
	  -lqop-mdwf3 $(QMP_LIBS)

$(programs): $(common:%.c=%.o)

$(f.objects): %-f.o: %.c
	$E CC $@
	$C$(CC) $(CFLAGS) -I$(TOP)/$(LIMPDIR) \
            -DQOP_MDWF_DEFAULT_PRECISION=\'F\' \
            -c -o $@ $<

$(d.objects): %-d.o: %.c
	$E CC $@
	$C$(CC) $(CFLAGS) -I$(TOP)/$(LIMPDIR) $(QMP_CFLAGS) \
            -DQOP_MDWF_DEFAULT_PRECISION=\'D\' \
            -c -o $@ $<

$(x.objects): %.o: %.c
	$E CC $@
	$C$(CC) $(CFLAGS) -I$(TOP)/$(LIMPDIR) \
            -c -o $@ $<

$(x.objects) $(d.objects) $(f.objects): $(headers)

