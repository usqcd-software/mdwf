TOP    = ..
V      = 0
CONFIG = macosx
QMP_CFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)

include $(TOP)/config/$(CONFIG)

x.sources = do_A.c \
            do_ApF.c \
            dot_fermion.c \
            get_fermion.c \
            madd_fermion.c \
            norm2_fermion.c \
            proj_Ug0plus.c \
            proj_Ug1plus.c \
            proj_Ug2plus.c \
            proj_Ug3plus.c \
            proj_g0minus.c \
            proj_g1minus.c \
            proj_g2minus.c \
            proj_g3minus.c \
            put_fermion.c \
            put_gauge.c \
            sizeof_gauge.c \
            sizeof_fermion.c \
            sizeof_pfermion.c

sources = fix_neighbor_f_down.c \
          fix_neighbor_f_up.c \
          get_neighbor.c \
          get_down_pack.c \
          get_down_pack_f.c \
          get_up_pack.c \
          get_up_pack_f.c \
          put_neighbor.c \
          put_down_pack.c \
          put_up_pack.c \
          put_abtable.c \
          sizeof_neighbor.c \
          sizeof_up_pack.c \
          sizeof_down_pack.c \
          sizeof_abtable.c 

headers = $(TOP)/mdwf.h \
          $(TOP)/qop-mdwf3.h

objects = $(x.sources:%.c=%f.o) \
          $(x.sources:%.c=%d.o) \
          $(sources:%.c=%.o)

.PHONY: all clean realclean

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

all: $(objects)
	$E AR $@ / faux
	$C$(AR) cr $(TOP)/$(library) $^
	$C$(RANLIB) $(TOP)/$(library)

clean:
	$E RM faux/objects
	$C$(RM) $(objects)

realclean: clean

$(sources:%.c=%.o): %.o: %.c
	$E CC $@
	$C$(CC) $(CFLAGS) -c -o $@ $<

$(x.sources:%.c=%f.o): %f.o: %.c
	$E CC-f $<
	$C$(CC) $(CFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'F\' -c -o $@ $<

$(x.sources:%.c=%d.o): %d.o: %.c
	$E CC-d $<
	$C$(CC) $(CFLAGS) -DQOP_MDWF_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(objects): $(headers)