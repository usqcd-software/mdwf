
xlc-1f.o:     file format elf32-powerpc-blrts

Disassembly of section .text:

00000000 <foo1>:
   0:	c0 04 00 00 	lfs     f0,0(r4)
   4:	c0 24 00 04 	lfs     f1,4(r4)
   8:	d0 03 00 00 	stfs    f0,0(r3)
   c:	d0 23 00 04 	stfs    f1,4(r3)
  10:	4e 80 00 20 	blr

00000014 <foo2>:
  14:	c0 04 00 00 	lfs     f0,0(r4)
  18:	c0 24 00 04 	lfs     f1,4(r4)
  1c:	c0 45 00 00 	lfs     f2,0(r5)
  20:	c0 65 00 04 	lfs     f3,4(r5)
  24:	ec 00 10 2a 	fadds   f0,f0,f2
  28:	ec 21 18 2a 	fadds   f1,f1,f3
  2c:	d0 03 00 00 	stfs    f0,0(r3)
  30:	d0 23 00 04 	stfs    f1,4(r3)
  34:	4e 80 00 20 	blr

00000038 <foo3>:
  38:	c0 04 00 04 	lfs     f0,4(r4)
  3c:	c0 25 00 00 	lfs     f1,0(r5)
  40:	c0 45 00 04 	lfs     f2,4(r5)
  44:	c0 64 00 00 	lfs     f3,0(r4)
  48:	ec 80 00 72 	fmuls   f4,f0,f1
  4c:	ec 00 00 b2 	fmuls   f0,f0,f2
  50:	ec 43 20 ba 	fmadds  f2,f3,f2,f4
  54:	ec 03 00 78 	fmsubs  f0,f3,f1,f0
  58:	d0 43 00 04 	stfs    f2,4(r3)
  5c:	d0 03 00 00 	stfs    f0,0(r3)
  60:	4e 80 00 20 	blr

00000064 <foo4>:
  64:	c0 04 00 00 	lfs     f0,0(r4)
  68:	c0 24 00 04 	lfs     f1,4(r4)
  6c:	c0 45 00 00 	lfs     f2,0(r5)
  70:	c0 65 00 04 	lfs     f3,4(r5)
  74:	ec 21 10 2a 	fadds   f1,f1,f2
  78:	ec 00 18 28 	fsubs   f0,f0,f3
  7c:	d0 23 00 04 	stfs    f1,4(r3)
  80:	d0 03 00 00 	stfs    f0,0(r3)
  84:	4e 80 00 20 	blr

00000088 <foo5>:
  88:	c0 44 00 00 	lfs     f2,0(r4)
  8c:	38 c0 00 04 	li      r6,4
  90:	c0 05 00 00 	lfs     f0,0(r5)
  94:	7c 44 31 1c 	lfssx   f2,r4,r6
  98:	7c 05 31 1c 	lfssx   f0,r5,r6
  9c:	00 01 10 24 	fxcpmadd f0,f1,f0,f2
  a0:	00 20 06 40 	fsmtp   f1,f0
  a4:	fc 00 00 18 	frsp    f0,f0
  a8:	fc 20 08 18 	frsp    f1,f1
  ac:	d0 03 00 00 	stfs    f0,0(r3)
  b0:	d0 23 00 04 	stfs    f1,4(r3)
  b4:	4e 80 00 20 	blr
