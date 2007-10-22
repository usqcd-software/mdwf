	.section __TEXT,__text,regular,pure_instructions
	.section __TEXT,__picsymbolstub1,symbol_stubs,pure_instructions,32
	.machine ppc
	.text
	.align 2
	.p2align 4,,15
	.globl _foo1
_foo1:
	lfd f0,8(r4)
	lfd f13,0(r4)
	stfd f0,8(r3)
	stfd f13,0(r3)
	blr
	.align 2
	.p2align 4,,15
	.globl _foo2
_foo2:
	lfd f0,8(r4)
	lfd f13,0(r4)
	fneg f0,f0
	fneg f13,f13
	stfd f0,8(r3)
	stfd f13,0(r3)
	blr
	.literal8
	.align 3
LC0:
	.long	0
	.long	0
	.text
	.align 2
	.p2align 4,,15
	.globl _foo3
_foo3:
	mflr r0
	lfd f12,0(r4)
	lfd f0,8(r4)
	bcl 20,31,"L00000000001$pb"
"L00000000001$pb":
	mflr r10
	mtlr r0
	addis r2,r10,ha16(LC0-"L00000000001$pb")
	lfd f11,lo16(LC0-"L00000000001$pb")(r2)
	fmsub f13,f12,f11,f0
	fmadd f0,f0,f11,f12
	stfd f13,0(r3)
	stfd f0,8(r3)
	blr
	.literal8
	.align 3
LC1:
	.long	-2147483648
	.long	0
	.text
	.align 2
	.p2align 4,,15
	.globl _foo4
_foo4:
	mflr r0
	lfd f12,0(r4)
	lfd f0,8(r4)
	bcl 20,31,"L00000000002$pb"
"L00000000002$pb":
	mflr r10
	mtlr r0
	addis r2,r10,ha16(LC1-"L00000000002$pb")
	lfd f11,lo16(LC1-"L00000000002$pb")(r2)
	fmadd f13,f12,f11,f0
	fmsub f0,f0,f11,f12
	stfd f13,0(r3)
	stfd f0,8(r3)
	blr
	.align 2
	.p2align 4,,15
	.globl _foo5
_foo5:
	lfd f13,8(r4)
	lfd f11,8(r5)
	lfd f0,0(r4)
	lfd f12,0(r5)
	fadd f13,f13,f11
	fadd f0,f0,f12
	stfd f13,8(r3)
	stfd f0,0(r3)
	blr
	.align 2
	.p2align 4,,15
	.globl _foo6
_foo6:
	lfd f13,8(r4)
	lfd f11,8(r5)
	lfd f0,0(r4)
	lfd f12,0(r5)
	fsub f13,f13,f11
	fsub f0,f0,f12
	stfd f13,8(r3)
	stfd f0,0(r3)
	blr
	.align 2
	.p2align 4,,15
	.globl _foo7
_foo7:
	lfd f12,8(r4)
	lfd f10,0(r5)
	lfd f0,8(r5)
	lfd f13,0(r4)
	fmul f11,f12,f10
	fmul f12,f12,f0
	fmadd f0,f13,f0,f11
	fmsub f13,f13,f10,f12
	stfd f0,8(r3)
	stfd f13,0(r3)
	blr
	.align 2
	.p2align 4,,15
	.globl _foo8
_foo8:
	lfd f0,8(r4)
	lfd f12,0(r5)
	lfd f10,8(r5)
	lfd f13,0(r4)
	fneg f0,f0
	fmul f11,f0,f10
	fmul f0,f0,f12
	fmsub f12,f13,f12,f11
	fmadd f13,f13,f10,f0
	stfd f12,0(r3)
	stfd f13,8(r3)
	blr
	.align 2
	.p2align 4,,15
	.globl _foo9
_foo9:
	lfd f12,8(r5)
	lfd f8,0(r6)
	lfd f0,8(r6)
	lfd f13,0(r5)
	lfd f10,8(r4)
	lfd f9,0(r4)
	fmul f11,f12,f8
	fmul f12,f12,f0
	fmadd f0,f13,f0,f11
	fmsub f13,f13,f8,f12
	fadd f0,f0,f10
	fadd f13,f13,f9
	stfd f0,8(r3)
	stfd f13,0(r3)
	blr
	.align 2
	.p2align 4,,15
	.globl _foo10
_foo10:
	lfd f0,8(r5)
	lfd f8,0(r6)
	lfd f13,8(r6)
	lfd f12,0(r5)
	lfd f10,8(r4)
	lfd f9,0(r4)
	fneg f0,f0
	fmul f11,f0,f8
	fmul f0,f0,f13
	fmadd f13,f12,f13,f11
	fmsub f12,f12,f8,f0
	fadd f13,f13,f10
	fadd f12,f12,f9
	stfd f13,8(r3)
	stfd f12,0(r3)
	blr
	.literal8
	.align 3
LC2:
	.long	-2147483648
	.long	0
	.text
	.align 2
	.p2align 4,,15
	.globl _foo11
_foo11:
	mflr r0
	lfd f13,0(r5)
	lfd f11,8(r5)
	bcl 20,31,"L00000000003$pb"
"L00000000003$pb":
	lfd f10,8(r4)
	lfd f12,0(r4)
	mflr r10
	mtlr r0
	addis r2,r10,ha16(LC2-"L00000000003$pb")
	lfd f9,lo16(LC2-"L00000000003$pb")(r2)
	fmsub f0,f11,f9,f13
	fmadd f13,f13,f9,f11
	fadd f0,f0,f10
	fadd f13,f13,f12
	stfd f0,8(r3)
	stfd f13,0(r3)
	blr
	.literal8
	.align 3
LC3:
	.long	-2147483648
	.long	0
	.text
	.align 2
	.p2align 4,,15
	.globl _foo12
_foo12:
	mflr r0
	lfd f13,0(r5)
	lfd f11,8(r5)
	bcl 20,31,"L00000000004$pb"
"L00000000004$pb":
	lfd f10,8(r4)
	lfd f12,0(r4)
	mflr r10
	mtlr r0
	addis r2,r10,ha16(LC3-"L00000000004$pb")
	lfd f9,lo16(LC3-"L00000000004$pb")(r2)
	fmsub f0,f11,f9,f13
	fmadd f13,f13,f9,f11
	fadd f0,f0,f10
	fadd f13,f13,f12
	stfd f0,8(r3)
	stfd f13,0(r3)
	blr
	.literal8
	.align 3
LC4:
	.long	0
	.long	0
	.text
	.align 2
	.p2align 4,,15
	.globl _foo13
_foo13:
	mflr r0
	lfd f9,0(r7)
	lfd f0,8(r7)
	bcl 20,31,"L00000000005$pb"
"L00000000005$pb":
	lfd f11,8(r4)
	lfd f10,0(r4)
	mflr r10
	mtlr r0
	addis r2,r10,ha16(LC4-"L00000000005$pb")
	lfd f13,lo16(LC4-"L00000000005$pb")(r2)
	fmul f12,f9,f13
	fmul f13,f0,f13
	fmadd f0,f1,f0,f12
	fmsub f1,f1,f9,f13
	fadd f0,f0,f11
	fadd f1,f1,f10
	stfd f0,8(r3)
	stfd f1,0(r3)
	blr
	.subsections_via_symbols
