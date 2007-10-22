	.file	"gcc-1.c"
	.section	".text"
	.align 2
	.globl foo1
	.type	foo1, @function
foo1:
	lfd 13,0(4)                   # load  Re(a)
	stfd 13,0(3)                  # store Re(a)
	lfd 0,8(4)                    # load  Im(a)
	stfd 0,8(3)                   # store Im(a)
	blr
	.size	foo1, .-foo1
	.align 2
	.globl foo2
	.type	foo2, @function
foo2:
	lfd 13,8(4)                   # f13 = Im(a)
	lfd 11,8(5)                   # f11 = Im(b)
	lfd 0,0(4)                    # f0 = Re(a)
	lfd 12,0(5)                   # f12 = Re(b)
	fadd 13,13,11                 # f13 = Im(a) + Im(b)
	fadd 0,0,12                   # f0 = Re(a) + Re(b)
	stfd 13,8(3)                  # store Im(result)
	stfd 0,0(3)                   # store Re(result)
	blr
	.size	foo2, .-foo2
	.align 2
	.globl foo3
	.type	foo3, @function
foo3:
	lfd 12,8(4)                   # f12 = Im(a)
	lfd 10,8(5)                   # f10 = Im(b)
	lfd 13,0(5)                   # f13 = Re(b)
	fmul 11,12,10                 # f11 = Im(a) * Im(b)
	lfd 0,0(4)                    # f0 = Re(a)
	fmul 12,13,12                 # f12 = Im(a) * Re(b)
	fmsub 13,0,13,11              # f13 = Re(a) * Re(b) - Im(a) * Im(b)
	fmadd 0,0,10,12               # f0 = Re(a) * Im(b) + Im(a) * Re(b)
	stfd 13,0(3)                  # store Re(result)
	stfd 0,8(3)                   # store Im(result)
	blr
	.size	foo3, .-foo3
	.section	.rodata
	.align 3
.LC0:
	.long	0
	.long	0
	.long	1072693248
	.long	0
	.section	".text"
	.align 2
	.globl foo4
	.type	foo4, @function
foo4:
	lis 9,.LC0@ha
	lfd 12,8(5)                  # f12 = Im(b)
	la 11,.LC0@l(9)
	lfd 9,.LC0@l(9)              # f9 = 0
	lfd 11,8(11)                 # f11 = 1
	fmul 0,9,12                  # .....
	lfd 10,0(5)
	lfd 13,8(4)
	fmul 12,12,11
	fmadd 11,10,11,0
	lfd 0,0(4)
	fmsub 10,10,9,12
	fadd 13,13,11
	fadd 0,0,10
	stfd 13,8(3)
	stfd 0,0(3)
	blr
	.size	foo4, .-foo4
	.section	.rodata.cst8,"aM",@progbits,8
	.align 3
.LC1:
	.long	0
	.long	0
	.section	".text"
	.align 2
	.globl foo5
	.type	foo5, @function
foo5:
	lis 9,.LC1@ha
	lfd 11,8(5)
	la 9,.LC1@l(9)
	lfd 10,0(5)
	lfd 13,0(9)
	lfd 12,8(4)
	fmul 9,10,13
	lfd 0,0(4)
	fmul 13,11,13
	fmadd 11,1,11,9
	fmsub 1,1,10,13
	fadd 12,12,11
	fadd 0,0,1
	stfd 12,8(3)
	stfd 0,0(3)
	blr
	.size	foo5, .-foo5
	.ident	"GCC: (GNU) 3.4.3"
