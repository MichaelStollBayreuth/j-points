	.file	"j-points-sift0-0.6-a.c"
	.version	"01.01"
gcc2_compiled.:
.text
	.align 16
.globl sift0
	.type	 sift0,@function
sift0:
	subl $24,%esp
	pushl %ebp
	pushl %edi
	pushl %esi
	pushl %ebx
	movl 56(%esp),%eax
	subl 52(%esp),%eax
	movl %eax,32(%esp)
	movl $0,36(%esp)
	movl 36(%esp),%edx
	cmpl sieve_primes1,%edx
	jge .Lx74
	movl 56(%esp),%eax
	notl %eax
	movl %eax,28(%esp)
	.align 4
.Lx76:
	movl 36(%esp),%edx
	movl pnn(,%edx,4),%eax
	sall $2,%eax
	movl prime(%eax),%esi
	movl sieves1+4(,%edx,8),%ecx
	cmpl $0,52(%esp)
	jg .Lx77
	movl 52(%esp),%eax
	negl %eax
	cltd
	idivl %esi
	movl %eax,%edi
	negl %edi
	jmp .Lx78
	.align 16
.Lx77:
	movl 52(%esp),%eax
	decl %eax
	cltd
	idivl %esi
	leal 1(%eax),%edi
.Lx78:
	cmpl $0,56(%esp)
	jge .Lx79
	movl 28(%esp),%eax
	cltd
	idivl %esi
	movl %eax,%ebp
	notl %ebp
	jmp .Lx80
	.align 16
.Lx79:
	movl 56(%esp),%eax
	cltd
	idivl %esi
	movl %eax,%ebp
.Lx80:
	movl survivors,%ebx
	cmpl %edi,%ebp
	jge .Lx81
	imull %ebp,%esi
	movl 52(%esp),%edx
	subl %esi,%edx
	leal (%ecx,%edx,4),%esi
	movl 32(%esp),%ecx
	testl %ecx,%ecx
	je .Lx75
	.align 4
.Lx85:
	movl (%esi),%eax
	andl %eax,(%ebx)
	addl $4,%esi
	addl $4,%ebx
	decl %ecx
	jnz .Lx85
	jmp .Lx75
	.align 16
.Lx81:
	movl %esi,%edx
	imull %edi,%edx
	movl 52(%esp),%eax
	subl %eax,%edx
	movl %esi,%eax
	subl %edx,%eax
	leal (%ecx,%eax,4),%ecx
	testl %edx,%edx
	je .Lx89
	.align 4
.Lx91:
	movl (%ecx),%eax
	andl %eax,(%ebx)
	addl $4,%ecx
	addl $4,%ebx
	decl %edx
	jnz .Lx91
.Lx89:
	leal 0(,%esi,4),%eax
	imull %ebp,%esi
	subl %edi,%ebp
	je .Lx94
	.align 4
.Lx96:
	movl %ecx,%edx
	subl %eax,%ecx
	.align 4
.Lx99:
	movl (%ecx),%edi
	andl %edi,(%ebx)
	addl $4,%ecx
	addl $4,%ebx
	cmpl %edx,%ecx
	jne .Lx99
	decl %ebp
	jnz .Lx96
.Lx94:
	subl %eax,%ecx
	movl 56(%esp),%edx
	subl %esi,%edx
	je .Lx75
	.align 4
.Lx105:
	movl (%ecx),%eax
	andl %eax,(%ebx)
	addl $4,%ecx
	addl $4,%ebx
	decl %edx
	jnz .Lx105
.Lx75:
	movl 36(%esp),%edx
	incl %edx
        movl %edx,36(%esp)
	cmpl sieve_primes1,%edx
	jl .Lx76
.Lx74:
	movl survivors,%eax
	movl %eax,24(%esp)
	movl 52(%esp),%ebp
	cmpl 56(%esp),%ebp
	jge .Lx109
	.align 4
.Lx111:
	movl 24(%esp),%edx
	movl (%edx),%edi
	addl $4,%edx
	movl %edx,24(%esp)
	testl %edi,%edi
	je .Lx110
	movl $sieves2p,%ecx
	testl %ebp,%ebp
	jge .Lx115
	movl $sieves2n,%ecx
.Lx115:
	incl num_surv1
	movl %ebp,%esi
	movl sieve_primes1,%eax
	movl sieve_primes2,%ebx
	subl %eax,%ebx
	je .Lx120
	.align 4
.Lx117:
	testl %edi,%edi
	je .Lx110
	movl %esi,%eax
	cltd
	idivl (%ecx)
	movl 4(%ecx),%eax
	andl (%eax,%edx,4),%edi
	addl $8,%ecx
	decl %ebx
	jnz .Lx117
.Lx120:
	testl %edi,%edi
	je .Lx110
	movl %ebp,%ebx
	sall $5,%ebx
	movl $1,%esi
	.align 4
.Lx125:
	testl %edi,%esi
	je .Lx126
	incl num_surv2
	movl 44(%esp),%eax
	movl 48(%esp),%ecx
	testl %ecx,%ecx
        je .Lgcd1end
	jge .Lx130
	negl %ecx
	.align 4
.Lx130:
	cltd
	idivl %ecx                 # division m/n
        movl %ecx,%eax             # m = n
	movl %edx,%ecx             # n = m_old%n 
	testl %ecx,%ecx            # n = 0 ?
	jne .Lx130
.Lgcd1end:                         # gcd(gcd,c)
        movl %ebx,%ecx             # %ecx = n = c
	testl %ecx,%ecx            # n = / > 0 ?
        je .Lgcd2end
	jge .L132a
	negl %ecx                  # n = -n
        .align 4
.L132a:
	cltd
	idivl %ecx                 # division m/n
        movl %ecx,%eax             # m = n
	movl %edx,%ecx             # n = m_old%n 
	testl %ecx,%ecx            # n = 0 ?
	jne .L132a
.Lgcd2end:
	cmpl $1,%eax
	jne .Lx126
	pushl %ebx
	movl 52(%esp),%eax
	pushl %eax
	movl 52(%esp),%edx
	pushl %edx
	call check_one_point
	addl $12,%esp
	testl %eax,%eax
	je .Lx126
	cmpl $0,one_point
	jne .Lx150
.Lx126:
	incl %ebx
	addl %esi,%esi
	jne .Lx125
	jmp .Lx110
	.align 16
.Lx150:
	movl $1,%eax
	jmp .Lx149
	.align 16
.Lx110:
	incl %ebp
	cmpl 56(%esp),%ebp
	jl .Lx111
.Lx109:
	xorl %eax,%eax
.Lx149:
	popl %ebx
	popl %esi
	popl %edi
	popl %ebp
	addl $24,%esp
	ret
.Lfe1:
	.size	 sift0,.Lfe1-sift0
	.align 16
.globl sift
	.type	 sift,@function
sift:
	subl $44,%esp
	pushl %ebp
	pushl %edi
	pushl %esi
	pushl %ebx
	movl height,%ebx
	addl $32,%ebx              # %ebx = c0 = height + 32
        movl 64(%esp),%ecx         # %ecx = a
	testl %ecx,%ecx            # a = 0 ?
	je .L0151
	movl 68(%esp),%eax         # %eax = b
	imull %eax,%eax            # %eax = b*b
	sall $2,%ecx               # %ecx = 4*a
	cltd
	idivl %ecx                 # %eax = (b*b)/(4*a), %edx = (b*b)%(4*a)
	testl %edx,%edx
	jne .L0151
	movl %eax,%ebx             # c0 = (b*b)/(4*a)
.L0151:
	movl %ebx,%eax
        sarl $5,%ebx
	movl %ebx,32(%esp)         # 32(%esp) = c0>>5
	andl $31,%eax
	movl %eax,40(%esp)         # 40(%esp) = c0 & 0x1f

	xorl %ebx,%ebx             # %ebx = n = 0
	movl sieve_primes1,%esi    # %esi = sieve_primes1
	movl %esi,24(%esp)         # 24(%esp) = sieve_primes1
	cmpl %esi,%ebx             # n >= sive_primes1 ?
	jge .L0153
	movl $sieves1+4,%ebp       # %ebp = &sieves1[0].ptr
	movl $0,16(%esp)           # 16(%esp) = 0
	.align 4
.L0155:
	movl pnn(,%ebx,4),%eax     # %eax = pn = pnn[n]
	sall $2,%eax
	movl prime(%eax),%edi      # %edi = p = prime[pn]
	movl 68(%esp),%eax         # %eax = b
	cltd
	idivl %edi
	movl %edx,%ecx             # %ecx = bp = b%p
	testl %ecx,%ecx            # bp >= 0 ?
	jge .L0156
	addl %edi,%ecx             # bp += p
.L0156:
	movl 64(%esp),%eax         # %eax = a
	cltd
	idivl %edi                 # %eax = a/p, %edx = a%p
	sall $14,%edx
	movl 16(%esp),%eax
	leal sieve(%eax,%edx),%eax # %eax = &sieve[a%p + 16%(esp)]
	sall $8,%ecx
	addl %eax,%ecx             # ... + 8*bp
	movl %ecx,(%ebp)           # sieves1[n].ptr = ...
        addl $8,%ebp
	addl $1048576,16(%esp)     # increment base address in 16(%esp)
	incl %ebx                  # n++
	cmpl %esi,%ebx             # n < sieve_primes1 ?
	jl .L0155
.L0153:
        movl sieve_primes2,%esi
	cmpl %esi,%ebx             # n >= sieve_primes2 ?
	jge .L0159
	movl $sieves2p+4,%ebp      # %ebp = &sieves2p[0].ptr
	.align 4
.L0161:
	movl pnn(,%ebx,4),%eax     # %eax = pn = pnn[n]
	sall $2,%eax
	movl prime(%eax),%edi      # %edi = p = prime[pn]
	movl 68(%esp),%eax         # %eax = b
	cltd
	idivl %edi
	movl %edx,%ecx             # %ecx = bp = b%p
	testl %ecx,%ecx            # bp >= 0 ?
	jge .L0162
	addl %edi,%ecx             # bp += p
.L0162:
	movl 64(%esp),%eax         # %eax = a
	cltd
	idivl %edi                 # %eax = a/p, %edx = a%p
	sall $14,%edx
	movl 16(%esp),%eax
	leal sieve(%eax,%edx),%eax
	sall $8,%ecx
	addl %eax,%ecx
	movl %ecx,(%ebp)           # sieves2p[n].ptr = ...
	leal (%ecx,%edi,4),%edi    # address + p
        movl %edi,160(%ebp)        # NB: 160 = sieves2n-sieves2p !!!
	addl $1048576,16(%esp)
	addl $8,%ebp
	incl %ebx
	cmpl %esi,%ebx
	jl .L0161
.L0159:
	movl height,%eax           # %eax = height
        movl %eax,%edx
	negl %edx
        sarl $5,%edx
	movl %edx,56(%esp)         # 56(%esp) = w_low = (-height)>>5
	sarl $5,%eax
	incl %eax
	movl %eax,52(%esp)         # 52(%esp) = w_high = (height>>5)+1
	movl %edx,16(%esp)         # 16(%esp) = w_low0 = w_low
	cmpl %eax,%edx
	jge .L0165
        movl $-1,%edx              # %edx = mask = ~zero
        testb $1,64(%esp)          # a&1
        jne .Lweiter
        testb $1,68(%esp)          # b&1
        jne .Lweiter
        movl $-1431655766,%edx     # mask = HALF_MASK
.Lweiter:
        movl %edx,36(%esp)         # 36(%esp) = mask
	.align 4
.L0167:
	movl 16(%esp),%ebx
        movl %ebx,%eax             # %eax = w_low0
	addl array_size,%ebx       # %ebx = w_high0 = w_low0 + array_size
        movl 52(%esp),%ecx         # %ecx = w_high
	cmpl %ecx,%ebx             # w_high0 > w_high ?
	jle .L0168
	movl %ecx,%ebx             # w_high0 = w_high
.L0168:
	movl survivors,%edx
	movl %ebx,%edi
	subl %eax,%edi             # %edi = i = range = w_high0 - w_low0
	# je .L0175                # have always range positive
        movl %edx,%ecx             # %ecx = surv = survivors
        movl 36(%esp),%esi         # %esi = mask
	.align 4
.L0173:
	movl %esi,(%ecx)           # *surv = mask
	addl $4,%ecx               # surv++
	decl %edi                  # i--
	jnz .L0173
.L0175:                            # now %ecx = &survivors[range]
	cmpl 56(%esp),%eax         # w_low0 = wlow ?
	jne .L0181
	movl begmask,%edi
	andl %edi,(%edx)           # survivors[0] &= begmask
.L0181:
	cmpl 52(%esp),%ebx         # w_high0 = w_high ?
	jne .L0182
	movl endmask,%edi
	andl %edi,-4(%ecx)         # survivors[range-1] &= endmask
.L0182:
	movl 32(%esp),%esi         # %esi = c0>>5
	cmpl %ebx,%esi             # %esi >= w_high0 ?
	jge .L0183
        subl %eax,%esi             # %esi = c0>>5 - w_low0
	jl .L0183
	movl 40(%esp),%ecx         # %ecx = c0 & 0x1f
        movl $1,%eax               # %eax = 1UL
        sall %cl,%eax              # (1UL)<<(c0 & 0x1f)
        notl %eax
	andl %eax,(%edx,%esi,4)    # survivors[...] &= ...
.L0183:
	pushl %ebx
	movl 20(%esp),%edx
	pushl %edx
	movl 76(%esp),%esi
	pushl %esi
	movl 76(%esp),%eax
	pushl %eax
	call sift0
	addl $16,%esp
	testl %eax,%eax
	je .L0166
	cmpl $0,one_point
	je .L0166
	movl $1,%eax
	jmp .L0186
	.align 16
.L0166:
	movl array_size,%edx
	addl %edx,16(%esp)
	movl 52(%esp),%esi
	cmpl %esi,16(%esp)
	jl .L0167
.L0165:
	xorl %eax,%eax
.L0186:
	popl %ebx
	popl %esi
	popl %edi
	popl %ebp
	addl $44,%esp
	ret
.Lfe2:
	.size	 sift,.Lfe2-sift
	.comm	sieves1,136,32
	.comm	sieves2p,136,32
	.comm	sieves2n,136,32
	.ident	"GCC: (GNU) egcs-2.91.66 19990314/Linux (egcs-1.1.2 release)"
