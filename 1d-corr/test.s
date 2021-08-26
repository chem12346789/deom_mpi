	.file	"test.cpp"
	.text
	.p2align 4,,15
	.globl	_Z3addRSt6vectorIiSaIiEES2_
	.type	_Z3addRSt6vectorIiSaIiEES2_, @function
_Z3addRSt6vectorIiSaIiEES2_:
.LFB1930:
	.cfi_startproc
	endbr64
	movq	(%rdi), %rax
	movq	(%rsi), %rdx
	movl	$10, (%rax)
	movl	$12, (%rdx)
	movl	(%rax), %eax
	addl	$12, %eax
	ret
	.cfi_endproc
.LFE1930:
	.size	_Z3addRSt6vectorIiSaIiEES2_, .-_Z3addRSt6vectorIiSaIiEES2_
	.section	.text.startup,"ax",@progbits
	.p2align 4,,15
	.globl	main
	.type	main, @function
main:
.LFB1931:
	.cfi_startproc
	endbr64
	movl	$0, 0
	ud2
	.cfi_endproc
.LFE1931:
	.size	main, .-main
	.p2align 4,,15
	.type	_GLOBAL__sub_I__Z3addRSt6vectorIiSaIiEES2_, @function
_GLOBAL__sub_I__Z3addRSt6vectorIiSaIiEES2_:
.LFB2450:
	.cfi_startproc
	endbr64
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	leaq	_ZStL8__ioinit(%rip), %rdi
	call	_ZNSt8ios_base4InitC1Ev@PLT
	movq	_ZNSt8ios_base4InitD1Ev@GOTPCREL(%rip), %rdi
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	leaq	__dso_handle(%rip), %rdx
	leaq	_ZStL8__ioinit(%rip), %rsi
	jmp	__cxa_atexit@PLT
	.cfi_endproc
.LFE2450:
	.size	_GLOBAL__sub_I__Z3addRSt6vectorIiSaIiEES2_, .-_GLOBAL__sub_I__Z3addRSt6vectorIiSaIiEES2_
	.section	.init_array,"aw"
	.align 8
	.quad	_GLOBAL__sub_I__Z3addRSt6vectorIiSaIiEES2_
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.hidden	__dso_handle
	.ident	"GCC: (Ubuntu 8.4.0-3ubuntu2) 8.4.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	 1f - 0f
	.long	 4f - 1f
	.long	 5
0:
	.string	 "GNU"
1:
	.align 8
	.long	 0xc0000002
	.long	 3f - 2f
2:
	.long	 0x3
3:
	.align 8
4:
