QMP_TOP = $$HOME/QDP/qmp-2.1.7

QMP_CFLAGS  = $(shell $(QMP_TOP)/bin/qmp-config --cflags)
QMP_LDFLAGS = $(shell $(QMP_TOP)/bin/qmp-config --ldflags)
QMP_LIBS    = $(shell $(QMP_TOP)/bin/qmp-config --libs)
CC          = gcc
CFLAGS      = -I.. -L.. -g
CLIBS       = -lqtest -lm

sources = operator.c

extras = 

programs = $(sources:%.c=%)

objects = $(sources:%.c=%.o) $(extras:%.c=%.o)

.PHONY: all clean realclean

all: $(programs)

clean:
	$(RM) $(objects)

realclean: clean
	$(RM) $(programs)

$(programs): %: %.o
	$(CC) $(CFLAGS) $(QMP_LDFLAGS) -o $@ $^ $(extras:%.c=%.o) $(CLIBS) $(QMP_LIBS)

$(programs): $(objects)

$(objects): %.o: %.c
	$(CC) $(CFLAGS) $(QMP_CFLAGS) -c -o $@ $^
