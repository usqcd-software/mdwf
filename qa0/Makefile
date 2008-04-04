TOP     = ..
V       = 0
CONFIG  = macosx
TARGETS = c99 c99-64 cee cee-64 bgl/xlc dry
TARGET  =
DEST    =
QA0     = qa0

include $(TOP)/config/$(CONFIG)

.PHONY: all clean realclean dist

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

DEST=$(TOP)/$(shell echo $(TARGET) | sed -e sx/x-xg)

headers = types.qa0 \
          defs-A.qa0 \
          defs-F.qa0 \
          defs-spin.qa0 \

sources = fix-down.qa0 \
          fix-up.qa0 \
          put-ab.qa0 \
          put-abi.qa0 \
          put-abi-z.qa0 \
          put-down.qa0 \
          get-down.qa0 \
          get-down-f.qa0 \
          put-up.qa0 \
          get-up.qa0 \
          get-up-f.qa0 \
          put-neighbor.qa0 \
          get-neighbor.qa0 \
          sizeof-ab-table.qa0 \
          sizeof-abi-table.qa0 \
          sizeof-down-pack.qa0 \
          sizeof-up-pack.qa0 \
          sizeof-neighbor.qa0

x.sources = get-fermion.qa0 \
            put-fermion.qa0 \
            put-gauge.qa0 \
            sizeof-fermion.qa0 \
            sizeof-gauge.qa0 \
            sizeof-pfermion.qa0 \
            f-copy.qa0 \
            f-norm.qa0 \
            f-add3.qa0 \
            f-add2.qa0 \
            f-add2-norm.qa0 \
            f-add2x.qa0 \
            f-dot.qa0 \
            proj-minus.qa0 \
            proj-plus.qa0 \
            proj-u-minus.qa0 \
            proj-u-plus.qa0 \
            do-A.qa0 \
            do-A-conj.qa0 \
            do-A-inv.qa0 \
            do-A-inv-conj.qa0 \
            do-F.qa0 \
            do-F-conj.qa0 \
            do-A-plus-F.qa0 \
            do-B-A-inv-F.qa0 \
            do-1-sub-B-A-inv-F.qa0 \
            do-1-sub-B-A-inv-F-norm.qa0 \
            do-A-conj-plus-B-conj-F-conj.qa0 \
            do-A-inv-conj-B-conj.qa0 \
            do-A-inv-conj-B-conj-F-conj.qa0 \

# Who include defs-A.qa0
A.sources = do-A.qa0 \
            do-A-conj.qa0 \
            do-A-inv.qa0 \
            do-A-conj-inv.qa0 \
            do-A-plus-F.qa0 \
            do-A-conj-plus-B-conj-F-conj.qa0 \
            do-B-A-inv-F.qa0 \
            do-1-sub-B-A-inv-F.qa0 \
            do-1-sub-B-A-inv-F-norm.qa0 \
            do-A-inv-conj-B-conj.qa0 \
            do-A-inv-conj-B-conj-F-conj.qa0 \

# Who include defs-F.qa0 (they also need defs-spin.qa0)
F.sources = do-F.qa0 \
            do-F-conj.qa0 \
            do-A-plus-F.qa0 \
            do-B-A-inv-F.qa0 \
            do-A-conj-plus-B-conj-F-conj.qa0 \
            do-1-sub-B-A-inv-F.qa0 \
            do-1-sub-B-A-inv-F-norm.qa0 \
            do-A-inv-conj-B-conj-F-conj.qa0 \

# Who include defs-spin.qa0
spin.sources = $(F.sources) \
               proj-plus.qa0 \
               proj-u-plus.qa0 \
               proj-minus.qa0 \
               proj-u-minus.qa0

# Convert qa0 to C names
dst.sources = $(sources:%.qa0=$(DEST)/%.c) \
              $(x.sources:%.qa0=$(DEST)/%d.c) \
              $(x.sources:%.qa0=$(DEST)/%f.c)

c.A.sources = $(A.sources:%.qa0=$(DEST)/%d.c) \
              $(A.sources:%.qa0=$(DEST)/%f.c)

c.F.sources = $(F.sources:%.qa0=$(DEST)/%d.c) \
              $(F.sources:%.qa0=$(DEST)/%f.c)

c.spin.sources = $(spin.sources:%.qa0=$(DEST)/%d.c) \
                 $(spin.sources:%.qa0=$(DEST)/%f.c)

.PHONY: do.all

dist all:
	@for t in $(TARGETS) ; do \
	  echo making target $$t; \
	  $(MAKE) "V=$(V)" CONFIG="$(CONFIG)" TARGET="$$t" do.all ; \
	done

do.all: $(dst.sources)

clean:

.PHONY: do.realclean

realclean: clean
	@for t in $(TARGETS); do \
	  echo making realclean $$t; \
	  $(MAKE) "V=$(V)" CONFIG="$(CONFIG)" TARGET="$$t" do.realclean ; \
	done

do.realclean:
	$E RM $(DEST)/sources
	$C $(RM) $(dst.sources)

$(sources:%.qa0=$(DEST)/%.c): $(DEST)/%.c: %.qa0
	$E QA0 $@
	$C $(QA0) $(TARGET) double $@ $<

$(x.sources:%.qa0=$(DEST)/%d.c): $(DEST)/%d.c: %.qa0
	$E QA0-d $@
	$C $(QA0) $(TARGET) double $@ $<

$(x.sources:%.qa0=$(DEST)/%f.c): $(DEST)/%f.c: %.qa0
	$E QA0-f $@
	$C $(QA0) $(TARGET) float $@ $<

$(dst.sources): types.qa0
$(c.A.sources): defs-A.qa0
$(c.F.sources): defs-F.qa0
$(c.spin.sources): defs-spin.qa0

