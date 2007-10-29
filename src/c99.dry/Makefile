TOP    = ..
V      = 0
CONFIG = macosx
DEST   = $(TOP)/c99.dry
QMP_CFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)

include $(TOP)/config/$(CONFIG)

.PHONY: all clean realclean

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

headers = ../mdwf.h

sources =  do-abd.c \
           do-abf.c \
           do-apfd.c \
           do-apff.c \
           dot-fermiond.c \
           dot-fermionf.c \
           fix-down.c \
           fix-up.c \
           get-down-f.c \
           get-fermiond.c \
           get-fermionf.c \
           get-up-f.c \
           madd-fermiond.c \
           madd-fermionf.c \
           norm2-fermiond.c \
           norm2-fermionf.c \
           proj-minusd.c \
           proj-minusf.c \
           proj-plusd.c \
           proj-plusf.c \
           proj-u-minusd.c \
           proj-u-minusf.c \
           proj-u-plusd.c \
           proj-u-plusf.c \
           put-ab.c \
           put-down.c \
           put-fermiond.c \
           put-fermionf.c \
           put-gauged.c \
           put-gaugef.c \
           put-neighbor.c \
           put-up.c \
           sizeof-ab-table.c \
           sizeof-down-pack.c \
           sizeof-fermiond.c \
           sizeof-fermionf.c \
           sizeof-gauged.c \
           sizeof-gaugef.c \
           sizeof-neighbor.c \
           sizeof-pfermiond.c \
           sizeof-pfermionf.c \
           sizeof-up-pack.c

objects = $(sources:%.c=%.o)

all: $(objects)
	$E AR $@ / c99.dry
	$C$(AR) cr $(TOP)/$(library) $^
	$C$(RANLIB) $(TOP)/$(library)

clean:
	$E RM c99.dry/objects
	$C$(RM) $(objects)

realclean: clean
	$E RM c99.dry/sources
	$C$(RM) $(sources)

$(sources:%.c=%.o): %.o: %.c
	$E CC $@
	$C$(CC) $(CFLAGS) -c -o $@ $<

