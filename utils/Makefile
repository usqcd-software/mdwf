SC = mzc

sources = ast.ss \
          attr.ss \
          basis.ss \
          cenv.ss \
          cfold.ss \
          common.ss \
          parser.ss \
          q2complex.ss \
          cheader.ss \
          qa0print.ss \
          verbose.ss \
          backend.ss \
          be-c99.ss

main = qa0.ss

program = $(main:%.ss=%)
objects = $(sources:%.ss=%.o)
kobjects = $(sources:%.ss=%.kp)

.PHONY: all clean realclean

all: $(program)

clean:
	$(RM) $(objects) $(kobjects)

realclean: clean
	$(RM) $(program)

$(objects): %.o: %.ss
	$(SC) -o $<

$(program): $(main) $(sources)
	$(SC) --exe $@ $(main)

$(objects): $(sources)
