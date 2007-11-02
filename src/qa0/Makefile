TOP    = ..
V      = 0
CONFIG = macosx
DEST   = $(TOP)/c99
QA0    = $(TOP)/../utils/qa0

include $(TOP)/config/$(CONFIG)

.PHONY: all clean realclean dist

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

headers = types.qa0

sources = fix-down.qa0 \
          fix-up.qa0 \
          put-ab.qa0 \
          put-down.qa0 \
          get-down.qa0 \
          get-down-f.qa0 \
          put-up.qa0 \
          get-up.qa0 \
          get-up-f.qa0 \
          put-neighbor.qa0 \
          get-neighbor.qa0 \
          sizeof-ab-table.qa0 \
          sizeof-down-pack.qa0 \
          sizeof-up-pack.qa0 \
          sizeof-neighbor.qa0

x.sources = do-ab.qa0 \
            do-apf.qa0 \
            dot-fermion.qa0 \
            get-fermion.qa0 \
            madd-fermion.qa0 \
            norm2-fermion.qa0 \
            proj-minus.qa0 \
            proj-plus.qa0 \
            proj-u-minus.qa0 \
            proj-u-plus.qa0 \
            put-fermion.qa0 \
            put-gauge.qa0 \
            sizeof-fermion.qa0 \
            sizeof-gauge.qa0 \
            sizeof-pfermion.qa0


c99.sources = $(sources:%.qa0=$(DEST)/%.c) \
              $(x.sources:%.qa0=$(DEST)/%d.c) \
              $(x.sources:%.qa0=$(DEST)/%f.c)

dist all: $(c99.sources)

clean:

realclean: clean
	$E RM c99/sources
	$C $(RM) $(c99.sources)

$(sources:%.qa0=$(DEST)/%.c): $(DEST)/%.c: %.qa0
	$E QA0 $@
	$C $(QA0) c99 double $@ $<

$(x.sources:%.qa0=$(DEST)/%d.c): $(DEST)/%d.c: %.qa0
	$E QA0-d $@
	$C $(QA0) c99 double $@ $<

$(x.sources:%.qa0=$(DEST)/%f.c): $(DEST)/%f.c: %.qa0
	$E QA0-f $@
	$C $(QA0) c99 float $@ $<


$(c99.sources): $(headers)
