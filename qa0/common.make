# common makefile for standard targets
TOP        = ..
V          = 0
LIBRARY    = libqop-mdwf3.a
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

headers = ../port/mdwf.h

i.sources = put-ab.c \
            put-abi.c \
            put-abi-z.c \
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
            sizeof-abi-table.c \
            sizeof-down-pack.c \
            sizeof-neighbor.c \
            sizeof-up-pack.c

x.sources = get-fermion \
            put-fermion \
            put-gauge \
            sizeof-fermion \
            sizeof-gauge \
            sizeof-pfermion \
            f-copy \
            f-norm \
            f-add3 \
            f-add2 \
            f-add2-norm \
            f-add2x \
	    f-dot \
            cg-xp \
            proj-minus \
            proj-plus \
            proj-u-minus \
            proj-u-plus \
            do-A \
            do-A-conj \
            do-A-inv \
            do-A-inv-conj \
            do-F \
            do-F-conj \
            do-A-plus-F \
            do-B-A-inv \
            do-B-A-inv-F \
            do-1-sub-B-A-inv-F \
            do-1-sub-F \
            do-1-sub-F-conj \
            do-1-sub-F-conj-norm \
            do-1-sub-B-A-inv-F-norm \
            do-A-conj-plus-B-conj-F-conj \
            do-A-inv-conj-B-conj \
            do-A-inv-conj-B-conj-F-conj \

sources = $(i.sources) \
          $(x.sources:%=%f.c) \
          $(x.sources:%=%d.c)

objects = $(sources:%.c=%.o)

all: $(objects)
	$E AR $@/$(LIMP)
	$C$(AR) cr $(LIBRARY) $^
	$C$(RANLIB) $(LIBRARY)
	$C$(MAKE) CONFIG='$(CONFIG)' V='$V' LIBRARY='$(LIBRARY)' \
                LIMPDIR=../$(LIMP) -C ../port $@

dist clean:
	$E RM $(LIMP)/objects
	$C$(RM) $(objects)
	$C$(MAKE) CONFIG='$(CONFIG)' V='$V' LIBRARY='$(LIBRARY)' \
                LIMPDIR=../$(LIMP) -C ../port $@

realclean: clean
	$E RM $(LIMP)/sources
	$C$(RM) $(sources)
	$C$(MAKE) CONFIG='$(CONFIG)' V='$V' LIBRARY='$(LIBRARY)' \
                LIMPDIR=../$(LIMP) -C ../port $@
	$C$(RM) $(LIBRARY)

$(sources:%.c=%.o): %.o: %.c
	$E CC $<
	$C$(CC) $(CFLAGS) -I../port -c -o $@ $<

