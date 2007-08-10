
xlc-1d.o:     file format elf32-powerpc-blrts

Disassembly of section .text:

00000000 <foo1>:
   0:	c8 04 00 00 	lfd     f0,0(r4)	# f0=Re(a)
   4:	c8 24 00 08 	lfd     f1,8(r4)	# f1=Im(a)	
   8:	d8 03 00 00 	stfd    f0,0(r3)	# Re(result) = Re(a)
   c:	d8 23 00 08 	stfd    f1,8(r3)	# Im(result) = Im(a)
  10:	4e 80 00 20 	blr

00000014 <foo2>:
  14:	c8 24 00 00 	lfd     f1,0(r4)	# f1 = Re(a)
  18:	38 c0 00 08 	li      r6,8
  1c:	c8 05 00 00 	lfd     f0,0(r5)	# f0 = Re(b)
  20:	7c 24 31 9c 	lfsdx   f1,r4,r6	# p1 = a
  24:	7c 05 31 9c 	lfsdx   f0,r5,r6	# p2 = b
  28:	00 01 00 18 	fpadd   f0,f1,f0	# p0 = a + b
  2c:	d8 03 00 00 	stfd    f0,0(r3)	# store Re(result)
  30:	7c 03 35 9c 	stfsdx  f0,r3,r6	# store Im(result)
  34:	4e 80 00 20 	blr

00000038 <foo3>:
  38:	c8 04 00 00 	lfd     f0,0(r4)	# f0 = Re(a)
  3c:	38 e0 00 08 	li      r7,8
  40:	3c c0 00 00 	lis     r6,0
     42: R_PPC_ADDR16_HA	.const_dr
  44:	7c 45 39 9c 	lfsdx   f2,r5,r7	# s2 = Im(b)
  48:	38 c6 00 00 	addi    r6,r6,0
     4a: R_PPC_ADDR16_LO	.const_dr
  4c:	7c 26 3b 1c 	lfpsx   f1,r6,r7	# p1 = 0 (?)
  50:	7c 04 39 9c 	lfsdx   f0,r4,r7	# p0 = a
  54:	c8 45 00 00 	lfd     f2,0(r5)	# p3 = b
  58:	10 22 08 3a 	fxcxnpma f1,f2,f0,f1	# p1 = (a*b)[0]
  5c:	00 02 08 24 	fxcpmadd f0,f2,f0,f1	# p1 = a*b
  60:	d8 03 00 00 	stfd    f0,0(r3)	# store Re(result)
  64:	7c 03 3d 9c 	stfsdx  f0,r3,r7	# store Im(result)
  68:	4e 80 00 20 	blr

0000006c <foo4>:
  6c:	c8 24 00 00 	lfd     f1,0(r4)	# f1 = Re(a)
  70:	3c c0 00 00 	lis     r6,0
     72: R_PPC_ADDR16_HA	.const_dr
  74:	38 e0 00 08 	li      r7,8
  78:	c8 05 00 00 	lfd     f0,0(r5)	# f0 = Re(b)
  7c:	38 c6 00 00 	addi    r6,r6,0
     7e: R_PPC_ADDR16_LO	.const_dr
  80:	7c 40 31 1c 	lfssx   f2,r0,r6	# s2 = 1
  84:	7c 24 39 9c 	lfsdx   f1,r4,r7	# p1 = a
  88:	7c 05 39 9c 	lfsdx   f0,r5,r7	# p0 = b
  8c:	10 02 08 3a 	fxcxnpma f0,f2,f0,f1	# p0 = a + b * I
  90:	d8 03 00 00 	stfd    f0,0(r3)	# store Re(result)
  94:	7c 03 3d 9c 	stfsdx  f0,r3,r7	# store Im(result)
  98:	4e 80 00 20 	blr

0000009c <foo5>:
  9c:	c8 44 00 00 	lfd     f2,0(r4)	# f2 = Re(a)
  a0:	38 c0 00 08 	li      r6,8
  a4:	c8 05 00 00 	lfd     f0,0(r5)	# f0 = Re(b)
  a8:	7c 44 31 9c 	lfsdx   f2,r4,r6	# p2 = a
  ac:	7c 05 31 9c 	lfsdx   f0,r5,r6	# p0 = b
  b0:	00 01 10 24 	fxcpmadd f0,f1,f0,f2	# f0 = a + scale * b
  b4:	d8 03 00 00 	stfd    f0,0(r3)	# store Re(result)
  b8:	7c 03 35 9c 	stfsdx  f0,r3,r6	# store Im(result)
  bc:	4e 80 00 20 	blr
