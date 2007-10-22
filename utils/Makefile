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

ast.o:        common.ss
attr.o:       common.ss ast.ss cenv.ss
cenv.o:       common.ss ast.ss basis.ss
cfold.o:      common.ss ast.ss cenv.ss basis.ss
parser.o:     common.ss ast.ss cenv.ss
q2complex.o:  common.ss ast.ss cenv.ss attr.ss
qa0print.o:   common.ss ast.ss
backend.o:    cenv.ss
be-c99.o:     common.ss ast.ss cenv.ss attr.ss cheader.ss
cheader.o:    common.ss ast.ss cenv.ss attr.ss
