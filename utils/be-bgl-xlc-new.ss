(module be-bgl-xlc-new
	mzscheme
  (require "common.ss")
  (require "ast.ss")
  (require "parser.ss")
  (require "attr.ss")
  (require "backend.ss")
  (require "cx2dh.ss")
  (require "cenv.ss")
  (require "cheader.ss")
  (require "verbose.ss")
  (require "be-ckind.ss")

  (provide machine-bgl/xlc)

  (define op-type-table
    '((dh-make                   dh-double)
      (dh-real                   double)
      (dh-imag                   double)
      (dh-move                   dh-double)
      (dh-neg                    dh-double)
      (dh-times-plus-i           dh-double)
      (dh-times-minus-i          dh-double)
      (dh-add                    dh-double)
      (dh-sub                    dh-double)
      (dh-rmul                   dh-double)
      (dh-rmadd                  dh-double)
      (dh-add-i                  dh-double)
      (dh-sub-i                  dh-double)
      (dh-norm-init              dh-double)
      (dh-norm-add               dh-double)
      (dh-norm-fini              double)
      (dh-dot-init               dh-double)
      (dh-dot-a                  dh-double)
      (dh-dot-b                  dh-double)
      (dh-dot-fini               dh-double)
      (dh-mul-a                  dh-double)
      (dh-mul-b                  dh-double)
      (dh-madd-a                 dh-double)
      (dh-madd-b                 dh-double)
      (dh-cmul-a                 dh-double)
      (dh-cmul-b                 dh-double)
      (dh-cmadd-a                dh-double)
      (dh-cmadd-b                dh-double)
      (dh->float                 dh-double)
      (double-add                double)
      (double-div                double)
      (double-madd               double)
      (double-move               double)
      (double-mul                double)
      (double-neg                double)
      (double-sub                double)
      (int-add                   int)
      (int-and                   int)
      (int-div                   int)
      (int-mod                   int)
      (int-move                  int)
      (int-mul                   int)
      (int-or                    int)
      (int-sub                   int)
      (int-xor                   int)
      (pointer-add               pointer)
      (pointer-move              pointer)
      (nop                       )))
  (define op-emit-table
    `((pointer-move            1 "$0 = %0"                                0)
      (pointer-add             2 "$0 = %0 + (%1)"                         0)
      (int-move                1 "$0 = %0"                                0)
      (int-mul                 2 "$0 = %0 * (%1)"                         0)
      (int-div                 2 "$0 = %0 / (%1)"                         0)
      (int-mod                 2 "$0 = %0 % (%1)"                         0)
      (int-add                 2 "$0 = %0 + (%1)"                         0)
      (int-sub                 2 "$0 = %0 - (%1)"                         0)
      (int-neg                 1 "$0 = -(%0)"                             0)
      (int-and                 2 "$0 = %0 & (%1)"                         0)
      (int-or                  2 "$0 = %0 | (%1)"                         0)
      (int-xor                 2 "$0 = %0 ^ (%1)"                         0)
      (int-not                 1 "$0 = !(%0)"                             0)
      (double-move             1 "$0 = %0"                                0)
      (double-neg              1 "$0 = - (%0)"                            1)
      (double-add              2 "$0 = %0 + (%1)"                         1)
      (double-sub              2 "$0 = %0 - (%1)"                         1)
      (double-div              2 "$0 = %0 / (%1)"                         1)
      (double-mul              2 "$0 = %0 * (%1)"                         1)
      (double-madd             3 "$0 = %0 + (%1) * (%2)"                  2)
      (dh-make                 2 "$0 = __cmplx(%0, %1)"                   0)
      (dh-real                 1 "$0 = __creal(%0)"                       0)
      (dh-imag                 1 "$0 = __cimag(%0)"                       0)
      (dh-move                 1 "$0 = (%0)"                              0)
      (dh-neg                  1 "$0 = __fpneg(%0)"                       2)
      (dh-times-plus-i         1 "$0 = __fxcxnpma(gZERO, %0, gONE)"       1)
      (dh-times-minus-i        1 "$0 = __fxcxnsma(gZERO, %0, gONE)"       1)
      (dh-add                  2 "$0 = __fpadd(%0, %1)"                   2)
      (dh-sub                  2 "$0 = __fpsub(%0, %1)"                   2)
      (dh-rmul                 2 "$0 = __fxpmul(%1, %0)"                  2)
      (dh-rmadd                3 "$0 = __fxcpmadd(%0, %2, %1)"            4)
      (dh-add-i                2 "$0 = __fxcxnpma(%0, %1, gONE)"          2)
      (dh-sub-i                2 "$0 = __fxcxnsma(%0, %1, gONE)"          2)
      (dh-norm-init            0 "$0 = gZERO"                             0)
      (dh-norm-add             2 "$0 = __fpmadd(%1, %1, %0)"              4)
      (dh-norm-fini            1 "$0 = __creal(%0) + __cimag(%0)"         1)
      (dh-dot-init             0 "$0 = gZERO"                             0)
      (dh-dot-a                3 "$0 = __fxcpmadd(%0, %2, __creal(%1))"   4)
      (dh-dot-b                3 "$0 = __fxcxnsma(%0, %2, __cimag(%1))"   4)
      (dh-dot-fini             1 "$0 = (%0)"                              0)
      (dh-mul-a                2 "$0 = __fxpmul(%1, __creal(%0))"         2)
      (dh-mul-b                3 "$0 = __fxcxnpma(%0, %2, __cimag(%1))"   4)
      (dh-madd-a               3 "$0 = __fxcpmadd(%0, %2, __creal(%1))"   4)
      (dh-madd-b               3 "$0 = __fxcxnpma(%0, %2, __cimag(%1))"   4)
      (dh-cmul-a               2 "$0 = __fxpmul(%1, __creal(%0))"         2)
      (dh-cmul-b               3 "$0 = __fxcxnsma(%0, %2, __cimag(%1))"   4)
      (dh-cmadd-a              3 "$0 = __fxcpmadd(%0, %2, __creal(%1))"   4)
      (dh-cmadd-b              3 "$0 = __fxcxnsma(%0, %2, __cimag(%1))"   4)
      (dh->float               1 "$0 = __fprsp(%0)"                       0)))
  (define load-table
    '((int               int)
      (pointer           pointer)
      (float             double)
      (double            double)
      (dh-float          dh-double)
      (dh-double         dh-double)))
  (define op-needing-zero*
    '(dh-times-plus-i
      dh-times-minus-i
      dh-norm-init
      dh-dot-init))
  (define op-needing-one*
    '(dh-times-plus-i
      dh-times-minus-i
      dh-add-i
      dh-sub-i))
  (define (emit-dh-double-store level addr* value env)
    (do-emit level (format "__stfpd((double *)(~a), ~a);"
			   (preemit-addr* addr* env)
			   (preemit-input value env))))
  (define (emit-dh-float-store level addr* value env)
    (do-emit level (format "__stfps((float *)(~a), ~a);"
			   (preemit-addr* addr* env)
			   (preemit-input value env))))
  (define (emit-dh-double-load level output addr* env)
    (do-emit level (format "~a = __lfpd((double *)(~a));"
			   (preemit-output output env)
			   (preemit-addr* addr* env))))
  (define (emit-dh-float-load level output addr* env)
    (do-emit level (format "~a = __lfps((float *)(~a));"
			   (preemit-output output env)
			   (preemit-addr* addr* env))))
  (define (emit-load level type output addr* env f)
    (case type
      [(dh-float) (emit-dh-float-load level output addr* env)]
      [(dh-double) (emit-dh-double-load level output addr* env)]
      [else (do-emit level "~a = *(~a *)(~a);"
		     (preemit-output output env)
		     (ce-lookup-x env 'name-of type "BGL/XLC name of type ~a"
				  type)
		     (preemit-addr* addr* env))])
    f)
  (define (emit-store level type addr* value env f)
    (case type
      [(dh-float) (emit-dh-float-store level addr* value env)]
      [(dh-double) (emit-dh-double-store level addr* value env)]
      [else (do-emit level "*(~a *)(~a) = ~a;"
		     (ce-lookup-x env 'name-of type "BGL/XLC name of type ~a"
				  type)
		     (preemit-addr* addr* env)
		     (preemit-input value env))])
    f)
  (define (collect-output* code* env)
    (define (add-output* output* type* env)
      (cond
       [(null? output*) env]
       [else (add-output* (cdr output*) (cdr type*)
			  (add-output (car output*) (car type*) env))]))
    (define (add-output output type env)
      (variant-case output
        [reg (name)
          (ce-search-x env 'back-end name
		       (lambda (v) env)
		       (lambda () (ce-bind-x env 'back-end name
					     (list (ck-new-var) type))))]))
    (define (add-zero name env)
      (if (not (memq name op-needing-zero*)) env
	  (ce-search-x env 'bgl/xlc 'zero
		       (lambda (x) env)
		       (lambda ()
			 (ce-bind-x env 'bgl/xlc 'zero (ck-new-var))))))
    (define (add-one name env)
      (if (not (memq name op-needing-one*)) env
	  (ce-search-x env 'bgl/xlc 'one
		       (lambda (x) env)
		       (lambda () (ce-bind-x env 'bgl/xlc 'one (ck-new-var))))))
    (walk-code* code*
		;; operation
		(lambda (env name attr* output* input*)
		  (let* ([env (add-zero name env)]
			 [env (add-one name env)])
		    (add-output* output* (be-out-type* env name) env)))
		;; load
		(lambda (env type attr* output addr*)
		  (add-output output (be-load-type env type) env))
		;; store
		(lambda (env type attr* addr* value) env)
		;; loop
		(lambda (env attr* var low high)
		  (add-output var 'int env))
		;; if
		(lambda (env var) env)
		env))
  (define (bgl/xlc-back-end ast env) (values (complex->double-hummer ast) env))
  (define (extra-env env)
    (let* ([env (ce-add-type env 'dh-double   "double _Complex" 16 16)]
	   [env (ce-add-type env 'dh-float    "float _Complex"   8  8)])
      env))
  (define (extra-decl* arg-name* arg-type* arg-c-name* arg-c-type* env)
    (define (zo-var*)
      (define (do-constant id init)
	(ce-search-x env 'bgl/xlc id
		     (lambda (v) (do-emit 1 "const double _Complex ~a = ~a;"
					  v init))
		     (lambda () #t)))
      (do-constant 'zero "__cmplx(0.0, 0.0)")
      (do-constant 'one  "__cmplx(1.0, 1.0)"))
    (define (emit-disjoint*)
      (do-emit 0 "")
      (let loop-a ([a* arg-name*] [t* arg-type*] [c* arg-c-name*])
	(cond
	 [(null? a*)]
	 [(not (eq? 'pointer (car t*))) (loop-a (cdr a*) (cdr t*) (cdr c*))]
	 [else
	  (let loop-b ([b* (cdr a*)] [s* (cdr t*)] [d* (cdr c*)])
	    (cond
	     [(null? b*) (loop-a (cdr a*) (cdr t*) (cdr c*))]
	     [(not (eq? 'pointer (car s*))) (loop-b (cdr b*) (cdr s*) (cdr d*))]
	     [else (do-emit 0 "#pragma disjoint(*~a,*~a)"
			    (preemit-param (car c*))
			    (preemit-param (car d*)))
		   (loop-b (cdr b*) (cdr s*) (cdr d*))]))])))
    (define (emit-align*)
      (do-emit 0 "")
      (let loop-a ([a* arg-name*] [t* arg-type*] [c* arg-c-name*])
	(cond
	 [(null? a*)]
	 [(not (eq? 'pointer (car t*))) (loop-a (cdr a*) (cdr t*) (cdr c*))]
	 [else (do-emit 1 "__alignx(16, ~a);" (preemit-param (car c*)))
	       (loop-a (cdr a*) (cdr t*) (cdr c*))])))
    ;;;
    (zo-var*)
    (emit-disjoint*)
    (emit-align*))

  (define (zo-define* env)
    (define (do-define id name def)
      (ce-search-x env 'bgl/xlc id
 		   (lambda (v) (do-emit 0 "#define ~a ~a" name (format def v)))
 		   (lambda () #t)))
     (do-define 'zero "gZERO" "(~a)")
     (do-define 'one "gONE" "(__cimag(~a))"))
  (define (zo-undefine* env)
    (define (do-undefine id name)
      (ce-search-x env 'bgl/xlc id
 		   (lambda (v) (do-emit 0 "#undef ~a" name))
 		   (lambda () #t)))
    (do-undefine 'zero "gZERO")
    (do-undefine 'one "gONE"))
  (define machine-bgl/xlc
    (build-ckind-back-end 'bgl/xlc          ; target-name
			  op-emit-table     ; op-emit-table
			  op-type-table     ; op-type-table
			  load-table        ; ld-type-table
			  emit-load         ; emit-load
			  emit-store        ; emit-store
			  collect-output*   ; collect-outputs
			  bgl/xlc-back-end  ; the-back-end
			  extra-env         ; extra-env
			  extra-decl*       ; extra-decl*
			  zo-define*        ; extra-def*
			  zo-undefine*)))   ; extra-undef*
