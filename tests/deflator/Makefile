TOP = ../..
V = 0
CONFIG = $(TOP)/target
library = qop-mdwf3
CDEBUG = -g
# PREC = D

TCFLAGS  = $(CFLAGS) -I$(TOP)/port -g
TLDFLAGS = $(LDFLAGS)
TCLIBS   = $(LIBS)

include $(CONFIG)

.PHONY: all clean realclean thelib dist

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null " "
   C=
endif

programs = check-vfsize \
           check-vgp \
           check-vf-dot-f \
           check-vf-dot-vf \
           check-vf-dot-vz \
           check-vf-dot-mz

programs.f = $(programs:%=%F)
programs.d = $(programs:%=%D)

programs.all = $(programs.f) $(programs.d)

utils = vfermion.c

headers = vfermion-test.h \
          $(TOP)/$(LIMPDIR)/qop-mdwf3.h

sources = $(programs:%=%.c) $(utils)

objects.f = $(sources:%.c=%F.o)
objects.d = $(sources:%.c=%D.o)

###############################################################################
all: $(programs.all)

clean:
	$E RM objects
	$C$(RM) $(objects.d) $(objects.f)

realclean: clean
	$E RM programs
	$C$(RM) $(programs.all)

dist: realclean

thelib:
	$(MAKE) -C $(TOP)


$(programs.f):%: %.o
	$E LD $@
	$C$(CC) $(CDEBUG) $(TLDFLAGS)  \
		-DQOP_MDWF_DEFAULT_PRECISION=\'F\' \
		-L$(TOP)/$(LIMPDIR)  -o $@ $< $(utils:%.c=%F.o) \
		-l$(library) $(TCLIBS)

$(programs.d):%: %.o
	$E LD $@
	$C$(CC) $(CDEBUG) $(TLDFLAGS)  \
		-DQOP_MDWF_DEFAULT_PRECISION=\'D\' \
		-L$(TOP)/$(LIMPDIR)  -o $@ $< $(utils:%.c=%D.o) \
		-l$(library) $(TCLIBS)

$(programs.f): $(utils:%.c=%F.o)
$(programs.d): $(utils:%.c=%D.o)

$(objects.f): %F.o: %.c
	$E CC $@
	$C$(CC) $(CDEBUG) -DQOP_MDWF_DEFAULT_PRECISION=\'F\' \
		$(TCFLAGS) -c -o $@ $<

$(objects.d): %D.o: %.c
	$E CC $@
	$C$(CC) $(CDEBUG) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' \
		$(TCFLAGS) -c -o $@ $<

$(objects.f): $(headers)
$(objects.d): $(headers)

$(programs.all): $(TOP)/$(LIMPDIR)/lib$(library).a
