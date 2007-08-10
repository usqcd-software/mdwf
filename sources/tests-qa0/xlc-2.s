
xlc-2.o:     file format elf32-powerpc-blrts

Disassembly of section .text:

00000000 <foo1>:
   0:	7c 00 23 9c 	lfpdx   f0,r0,r4
   4:	7c 00 1f 9c 	stfpdx  f0,r0,r3
   8:	4e 80 00 20 	blr

0000000c <foo2>:
   c:	7c 00 23 9c 	lfpdx   f0,r0,r4
  10:	00 00 01 40 	fpneg   f0,f0
  14:	7c 00 1f 9c 	stfpdx  f0,r0,r3
  18:	4e 80 00 20 	blr

0000001c <foo3>:
  1c:	7c 20 23 9c 	lfpdx   f1,r0,r4
  20:	3c 80 00 00 	lis     r4,0
			22: R_PPC_ADDR16_HA	.const_dr
  24:	38 84 00 00 	addi    r4,r4,0
			26: R_PPC_ADDR16_LO	.const_dr
  28:	c0 04 00 04 	lfs     f0,4(r4)
  2c:	7c 00 21 1c 	lfssx   f0,r0,r4
  30:	00 01 00 12 	fxmul   f0,f1,f0
  34:	7c 00 1f 9c 	stfpdx  f0,r0,r3
  38:	4e 80 00 20 	blr

0000003c <foo4>:
  3c:	7c 20 23 9c 	lfpdx   f1,r0,r4
  40:	3c 80 00 00 	lis     r4,0
			42: R_PPC_ADDR16_HA	.const_dr
  44:	38 a0 00 04 	li      r5,4
  48:	38 84 00 00 	addi    r4,r4,0
			4a: R_PPC_ADDR16_LO	.const_dr
  4c:	c0 04 00 00 	lfs     f0,0(r4)
  50:	7c 04 29 1c 	lfssx   f0,r4,r5
  54:	00 01 00 12 	fxmul   f0,f1,f0
  58:	7c 00 1f 9c 	stfpdx  f0,r0,r3
  5c:	4e 80 00 20 	blr

00000060 <foo5>:
  60:	7c 20 23 9c 	lfpdx   f1,r0,r4
  64:	7c 00 2b 9c 	lfpdx   f0,r0,r5
  68:	00 01 00 18 	fpadd   f0,f1,f0
  6c:	7c 00 1f 9c 	stfpdx  f0,r0,r3
  70:	4e 80 00 20 	blr

00000074 <foo6>:
  74:	7c 20 23 9c 	lfpdx   f1,r0,r4
  78:	7c 00 2b 9c 	lfpdx   f0,r0,r5
  7c:	00 01 00 1a 	fpsub   f0,f1,f0
  80:	7c 00 1f 9c 	stfpdx  f0,r0,r3
  84:	4e 80 00 20 	blr

00000088 <foo7>:
  88:	7c 00 2b 9c 	lfpdx   f0,r0,r5
  8c:	7c 40 23 9c 	lfpdx   f2,r0,r4
  90:	00 22 00 14 	fxpmul  f1,f2,f0
  94:	10 02 08 32 	fxcsnpma f0,f2,f0,f1
  98:	7c 00 1f 9c 	stfpdx  f0,r0,r3
  9c:	4e 80 00 20 	blr

000000a0 <foo8>:
  a0:	7c 00 2b 9c 	lfpdx   f0,r0,r5
  a4:	7c 40 23 9c 	lfpdx   f2,r0,r4
  a8:	00 22 00 14 	fxpmul  f1,f2,f0
  ac:	10 02 08 36 	fxcsnsma f0,f2,f0,f1
  b0:	7c 00 1f 9c 	stfpdx  f0,r0,r3
  b4:	4e 80 00 20 	blr

000000b8 <foo9>:
  b8:	7c 40 23 9c 	lfpdx   f2,r0,r4
  bc:	7c 00 33 9c 	lfpdx   f0,r0,r6
  c0:	7c 20 2b 9c 	lfpdx   f1,r0,r5
  c4:	00 41 10 24 	fxcpmadd f2,f1,f0,f2
  c8:	10 01 10 32 	fxcsnpma f0,f1,f0,f2
  cc:	7c 00 1f 9c 	stfpdx  f0,r0,r3
  d0:	4e 80 00 20 	blr

000000d4 <foo10>:
  d4:	7c 40 23 9c 	lfpdx   f2,r0,r4
  d8:	7c 00 33 9c 	lfpdx   f0,r0,r6
  dc:	7c 20 2b 9c 	lfpdx   f1,r0,r5
  e0:	00 41 10 24 	fxcpmadd f2,f1,f0,f2
  e4:	10 01 10 36 	fxcsnsma f0,f1,f0,f2
  e8:	7c 00 1f 9c 	stfpdx  f0,r0,r3
  ec:	4e 80 00 20 	blr

000000f0 <foo11>:
  f0:	7c 20 23 9c 	lfpdx   f1,r0,r4
  f4:	3c 80 00 00 	lis     r4,0
			f6: R_PPC_ADDR16_HA	.const_dr
  f8:	7c 00 2b 9c 	lfpdx   f0,r0,r5
  fc:	38 84 00 00 	addi    r4,r4,0
			fe: R_PPC_ADDR16_LO	.const_dr
 100:	c0 44 00 00 	lfs     f2,0(r4)
 104:	10 02 08 30 	fxcpnpma f0,f2,f0,f1
 108:	7c 00 1f 9c 	stfpdx  f0,r0,r3
 10c:	4e 80 00 20 	blr

00000110 <foo12>:
 110:	7c 20 23 9c 	lfpdx   f1,r0,r4
 114:	3c 80 00 00 	lis     r4,0
			116: R_PPC_ADDR16_HA	.const_dr
 118:	7c 00 2b 9c 	lfpdx   f0,r0,r5
 11c:	38 84 00 00 	addi    r4,r4,0
			11e: R_PPC_ADDR16_LO	.const_dr
 120:	c0 44 00 00 	lfs     f2,0(r4)
 124:	10 02 08 34 	fxcpnsma f0,f2,f0,f1
 128:	7c 00 1f 9c 	stfpdx  f0,r0,r3
 12c:	4e 80 00 20 	blr

00000130 <foo13>:
 130:	7c 40 23 9c 	lfpdx   f2,r0,r4
 134:	7c 00 2b 9c 	lfpdx   f0,r0,r5
 138:	10 01 10 30 	fxcpnpma f0,f1,f0,f2
 13c:	7c 00 1f 9c 	stfpdx  f0,r0,r3
 140:	4e 80 00 20 	blr
 144:	00 00 00 00 	.long 0x0
