(module be-dry
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

  (provide machine-dry)

  (define op-type-table
    '((complex-add               complex-double)
      (complex-sub               complex-double)
      (complex-add-i             complex-double)
      (complex-sub-i             complex-double)
      (complex-times-plus-i      complex-double)
      (complex-times-minus-i     complex-double)
      (complex-rmadd             complex-double)
      (complex-rmsub             complex-double)
      (complex-cmadd             complex-double)
      (complex-cmul              complex-double)
      (complex-madd              complex-double)
      (complex-move              complex-double)
      (complex-rmul              complex-double)
      (complex-mul               complex-double)
      (complex-neg               complex-double)
      (complex                   complex-double)
      (complex-real              double)
      (complex-imag              double)
      (complex-norm-init         double)
      (complex-norm-add          double)
      (complex-norm-fini         double)
      (complex-dot-init          complex-double)
      (complex-dot-add           complex-double)
      (complex-dot-fini          complex-double)
      (double-add                double)
      (double-div                double)
      (double-madd               double)
      (double-msub               double)
      (double-zero               double)
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
    ;; op i-count fmt flop-count
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
      (double-zero             0 "$0 = 0.0"                               0)
      (double-move             1 "$0 = %0"                                0)
      (double-neg              1 "$0 = - (%0)"                            1)
      (double-add              2 "$0 = %0 + (%1)"                         1)
      (double-sub              2 "$0 = %0 - (%1)"                         1)
      (double-div              2 "$0 = %0 / (%1)"                         1)
      (double-mul              2 "$0 = %0 * (%1)"                         1)
      (double-madd             3 "$0 = %0 + (%1) * (%2)"                  2)
      (double-msub             3 "$0 = %0 - (%1) * (%2)"                  2)
      (complex-move            1 "$0 = %0"                                0)
      (complex                 2 "$0 = %0 + I * (%1)"                     0)
      (complex-real            1 "$0 = creal(%0)"                         0)
      (complex-imag            1 "$0 = cimag(%0)"                         0)
      (complex-neg             1 "$0 = -(%0)"                             2)
      (complex-times-plus-i    1 "$0 = I * (%0)"                          1)
      (complex-times-minus-i   1 "$0 = -I * (%0)"                         1)
      (complex-add             2 "$0 = %0 + (%1)"                         2)
      (complex-sub             2 "$0 = %0 - (%1)"                         2)
      (complex-rmul            2 "$0 = %0 * (%1)"                         2)
      (complex-mul             2 "$0 = %0 * (%1)"                         6)
      (complex-madd            3 "$0 = %0 + (%1) * (%2)"                  8)
      (complex-rmadd           3 "$0 = %0 + (%1) * (%2)"                  4)
      (complex-rmsub           3 "$0 = %0 - (%1) * (%2)"                  4)
      (complex-cmul            2 "$0 = conj(%0) * (%1)"                   6)
      (complex-cmadd           3 "$0 = %0 + conj(%1) * (%2)"              8)
      (complex-add-i           2 "$0 = %0 + I * (%1)"                     2)
      (complex-sub-i           2 "$0 = %0 - I * (%1)"                     2)
      (complex-norm-init       0 "$0 = 0.0"                               0)
      (complex-norm-add        2 "$0 = %0 + creal(conj(%1) * (%1))"       4)
      (complex-norm-fini       1 "$0 = %0"                                0)
      (complex-dot-init        0 "$0 = 0.0"                               0)
      (complex-dot-add         3 "$0 = %0 + conj(%1) * (%2)"              8)
      (complex-dot-fini        1 "$0 = (%0)"                              0)))
  (define load-table
    '((int               int)
      (pointer           pointer)
      (float             double)
      (double            double)
      (complex-float     complex-double)
      (complex-double    complex-double)))
  (define (emit-load level type output addr* env f)
    (let ([value (case type
		   [(int) "0"]
		   [(pointer) "(void *)0"]
		   [(float double)   "0.0"]
		   [(complex-float complex-double) "0.0"]
		   [else (error 'be-dry "emit-load type=~a, addr*=~a"
				type addr*)])])
      (do-emit level "~a = ~a;"
	       (preemit-output output env)
	       value))
    f)
  (define (emit-store level type addr* value env f)
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
				   (lambda ()
				     (ce-bind-x env 'back-end name
						(list (ck-new-var) type))))]))
    (walk-code* code*
		; operation
		(lambda (env name attr* output* input*)
		  (add-output* output* (be-out-type* env name) env))
		; load
		(lambda (env type attr* output addr*)
		  (add-output output (be-load-type env type) env))
		; store
		(lambda (env type attr* addr* value) env)
		; loop
		(lambda (env attr* var low high)
		  (add-output var 'int env))
		; if
		(lambda (env var) env)
		env))
  (define (c99-back-end ast env) (values ast env))
  (define (extra-env env) env)
  (define (extra-decl* arg-name* arg-type* arg-c-name* arg-c-type* env) #t)
  (define (extra-def* env) #t)
  (define (extra-undef* env) #t)
  (define machine-dry
    (build-ckind-back-end 'dry              ; target-name
                          4                 ; pointer-size
                          4                 ; pointer-align
			  op-emit-table     ; op-emit-table
			  op-type-table     ; op-type-table
			  load-table        ; ld-type-table
			  emit-load         ; emit-load
			  emit-store        ; emit-store
			  collect-output*   ; collect-outputs
			  c99-back-end      ; the-back-end
			  extra-env         ; extra-env
			  extra-decl*       ; extra-decl*
			  extra-def*        ; extra-def*
			  empty-postparam*  ; extra-postparam*
			  extra-undef*)))   ; extra-undef*
