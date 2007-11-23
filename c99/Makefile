TOP    = ..
V      = 0
CONFIG = macosx
QMP_CFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)

include $(TOP)/config/$(CONFIG)

.PHONY: all clean realclean dist

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

headers = ../mdwf.h

i.sources = put-ab.c \
            put-down.c \
            get-down.c \
            get-down-f.c \
            put-up.c \
            get-up.c \
            get-up-f.c \
            put-neighbor.c \
            get-neighbor.c \
            fix-down.c \
            fix-up.c \
            sizeof-ab-table.c \
            sizeof-down-pack.c \
            sizeof-neighbor.c \
            sizeof-up-pack.c

x.sources = do-A \
            do-A-conj \
            do-A-inv \
            do-A-conj-inv \
            do-ApF \
            do-1pA \
            do-ApB \
            do-AxpBxFx \
            do-F \
            dot-fermion \
            madd-fermion \
            norm2-fermion \
            proj-minus \
            proj-plus \
            proj-u-minus \
            proj-u-plus \
            put-fermion \
            get-fermion \
            put-gauge \
            sizeof-fermion \
            sizeof-gauge \
            sizeof-pfermion

sources = $(i.sources) \
          $(x.sources:%=%f.c) \
          $(x.sources:%=%d.c)

objects = $(sources:%.c=%.o)

all: $(objects)
	$E AR $@/c99
	$C$(AR) cr $(TOP)/$(library) $^
	$C$(RANLIB) $(TOP)/$(library)

dist clean:
	$E RM c99/objects
	$C$(RM) $(objects)

realclean: clean
	$E RM c99.dry/sources
	$C$(RM) $(sources)

$(sources:%.c=%.o): %.o: %.c
	$E CC $<
	$C$(CC) $(CFLAGS) -c -o $@ $<

