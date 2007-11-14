(module be-c99
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

  (provide machine-c99-32)

  (define *var-count* 0)
  (define (new-var)
    (let ([x (gen-reg 'g *var-count*)])
      (set! *var-count* (+ *var-count* 1))
      x))
  (define (preemit-param p)
    (format "a_~a" p))

  (define c99-op*
    '((complex-add               complex-double)
      (complex-sub               complex-double)
      (complex-add-i             complex-double)
      (complex-sub-i             complex-double)
      (complex-times-plus-i      complex-double)
      (complex-times-minus-i     complex-double)
      (complex-rmadd             complex-double)
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
  (define c99-load*
    '((int               int)
      (pointer           pointer)
      (float             double)
      (double            double)
      (complex-float     complex-double)
      (complex-double    complex-double)))
  (define (c99-back-end ast env) (values ast env))
  (define (c99-emit qa0 env)
    (variant-case qa0
      [qa0-top (decl*) (let loop ([decl* decl*])
			 (cond
			  [(null? decl*)]
			  [else (emit-decl (car decl*) env)
				(loop (cdr decl*))]))]))
  (define (emit-decl decl env)
    (variant-case decl
      [qa0-proc (attr* name arg-name* arg-type* arg-c-name* arg-c-type* code*)
	(emit-proc attr* name arg-name* arg-type* arg-c-name* arg-c-type*
		   code* env)]
      [qa0-verbose (target* data*)
	(emit-verbose 'c99 target* data*)]
      [else #t]))
  (define (collect-outputs code* env)
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
						(list (new-var) type))))]))
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
  (define (check-inputs! code* env)
    (define (chk-input*! input* env)
      (cond
       [(null? input*) env]
       [else (chk-input! (car input*) env)
	     (chk-input*! (cdr input*) env)]))
    (define (chk-input! input env)
      (variant-case input
        [reg (name)
	  (ce-lookup-x env 'back-end name
		       "C99: register ~a is not defined" name)
	  env]
	[else env]))
    (walk-code* code*
		(lambda (env name attr* output* input*)
		  (chk-input*! input* env))
		; load
		(lambda (env type attr* output addr*)
		  (chk-input*! addr* env))
		; store
		(lambda (env type attr* addr* value)
		  (chk-input*! addr* env)
		  (chk-input! value env))
		; loop
		(lambda (env attr* var low high)
		  (chk-input! var env)
		  (chk-input! low env)
		  (chk-input! high env))
		; if
		(lambda (env var)
		  (chk-input! var env))
		env))
  (define (do-emit level fmt . arg*)
    (let loop ([n level])
      (cond
       [(zero? n)]
       [else (printf "  ") (loop (- n 1))]))
    (apply printf fmt arg*)
    (newline))
  (define (emit-proc-decl cl? rt attr* name arg-c-name* arg-c-type* env)
    (printf "~a~%" (build-proc-type name attr* env))
    (let* ([v (build-proc-name name attr* env)]
	   [x (make-string (+ 1 (string-length v)) #\space)])
      (printf "~a" v)
      (if (null? arg-c-name*) (printf "(void)~%")
	  (begin
	    (printf "(")
	    (let loop ([name* arg-c-name*] [type* arg-c-type*] [p ""])
	      (cond
	       [(null? name*)]
	       [else (printf "~a~a ~a~a" p (car type*)
			     (preemit-param (car name*))
			     (if (null? (cdr name*)) "" ",\n"))
		     (loop (cdr name*) (cdr type*) x)]))
	    (printf ")~%")))))
  (define (emit-variables env)
    (ce-for-each env
		 (lambda (k v) (and (list? k) (eq? (car k) 'back-end)))
		 (lambda (k v)
		   (do-emit 1 "~a ~a;"
			    (ce-lookup-x env 'name-of (cadr v)
					 "C99: type for ~a" v)
			    (car v)))))
  (define (emit-param* arg-name* arg-c-name* env)
    (do-emit 0 "")
    (let loop ([n* arg-name*] [c* arg-c-name*])
      (cond
       [(null? n*)]
       [else (let* ([rn (ce-lookup-x env 'back-end (car n*)
				     "C99: name for ~a" (car n*))]
		    [rt (ce-lookup-x env 'name-of (cadr rn)
				     "C99: type for ~a" (car n*))])
	       (do-emit 1 "~a = (~a)~a;" (car rn) rt (preemit-param (car c*)))
	       (loop (cdr n*) (cdr c*)))])))
  (define (preemit-output output env)
    (variant-case output
      [reg (name) (car (ce-lookup-x env 'back-end name
				    "C99: name of ~a not found" name))]))
  (define (preemit-input input env)
    (variant-case input
      [reg (name) (car (ce-lookup-x env 'back-end name
				    "C99: name of ~a not found" name))]
      [c-expr-number (number) number]))
  (define (preemit-addr* addr* env)
    (let loop ([r (preemit-input (car addr*) env)] [addr* (cdr addr*)])
      (cond
       [(null? addr*) r]
       [else (loop (format "~a + (~a)" r (preemit-input (car addr*) env))
		   (cdr addr*))])))
  (define (emit-count level cf? counter f)
    (if (and cf? (not (zero? f)))
	(do-emit level "~a += ~a; /* count flops */" counter f)))
  (define op-table
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
      (complex-move            1 "$0 = %0;"                               0)
      (complex                 2 "$0 = %0 + I * (%1);"                    0)
      (complex-real            1 "$0 = creal(%0);"                        0)
      (complex-imag            1 "$0 = cimag(%0);"                        0)
      (complex-neg             1 "$0 = -(%0);"                            2)
      (complex-times-plus-i    1 "$0 = I * (%0);"                         1)
      (complex-times-minus-i   1 "$0 = -I * (%0);"                        1)
      (complex-add             2 "$0 = %0 + (%1);"                        2)
      (complex-sub             2 "$0 = %0 - (%1);"                        2)
      (complex-rmul            2 "$0 = %0 * (%1);"                        2)
      (complex-mul             2 "$0 = %0 * (%1);"                        6)
      (complex-madd            3 "$0 = %0 + (%1) * (%2);"                 8)
      (complex-rmadd           3 "$0 = %0 + (%1) * (%2);"                 4)
      (complex-cmul            2 "$0 = conj(%0) * (%1);"                  6)
      (complex-cmadd           3 "$0 = %0 + conj(%1) * (%2);"             8)
      (complex-add-i           2 "$0 = %0 + I * (%1);"                    2)
      (complex-sub-i           2 "$0 = %0 - I * (%1);"                    2)
      (complex-norm-init       0 "$0 = 0.0;"                              0)
      (complex-norm-add        2 "$0 = %0 + creal(conj(%1) * (%1));"      4)
      (complex-norm-fini       1 "$0 = %0;"                               0)
      (complex-dot-init        0 "$0 = 0.0;"                              0)
      (complex-dot-add         3 "$0 = %0 + conj(%1) * (%2);"             8)
      (complex-dot-fini        1 "$0 = (%0);"                             0)))
  (define (get-element r* f)
    (if (char-numeric? f) (format "~a" (list-ref r* (- (char->integer f)
						       (char->integer #\0))))
	(string f)))
  (define (do-emit-op level name out* in* outx* inx* in-count fmt flops more)
    (let loop ([r ""] [f* (string->list fmt)])
      (cond
       [(null? f*) (do-emit level "~a;" r) (+ flops more)]
       [(and (eq? (car f*) #\$) (not (null? (cdr f*))))
	(loop (string-append r (get-element outx* (cadr f*))) (cddr f*))]
       [(and (eq? (car f*) #\%) (not (null? (cdr f*))))
	(loop (string-append r (get-element inx* (cadr f*))) (cddr f*))]
       [else (loop (string-append r (string (car f*))) (cdr f*))])))
  (define (emit-op level name attr* output* input* env f)
    (let ([in* (map (lambda (in) (preemit-input in env)) input*)]
	  [out* (map (lambda (out) (preemit-output out env)) output*)])
      (cond
       [(eq? name 'nop)
	(do-emit level "/* NOP: ~a */" (map qa0-attr->name attr*)) f]
       [(assq name op-table)
	=> (lambda (op)
	     (let ([i-count (cadr op)] [fmt (caddr op)] [f-count (cadddr op)])
	       (do-emit-op level name output* input* out* in* i-count fmt
			   f f-count)))]
       [else (error 'qa0-c99 "UNKNOWN op ~a, out* ~a, in* ~a, attr* ~a"
		    name output* input* attr*)])))
  (define (emit-load level type output addr* env f)
    (do-emit level "~a = *(~a *)(~a);"
	     (preemit-output output env)
	     (ce-lookup-x env 'name-of type "C99: name of type ~a" type)
	     (preemit-addr* addr* env))
    f)
  (define (emit-store level type addr* value env f)
    (do-emit level "*(~a *)(~a) = ~a;"
	     (ce-lookup-x env 'name-of type "C99: name of type ~a" type)
	     (preemit-addr* addr* env)
	     (preemit-input value env))
    f)
  (define (emit-if level var true-code* false-code* env cf? counter f)
    (emit-count level cf? counter f)
    (do-emit level "if (~a) {" (preemit-input var env))
    (emit-code* (+ level 1) true-code* env cf? counter 0)
    (if (not (null? false-code*))
	(begin (do-emit level "} else {")
	       (emit-code* (+ level 1) false-code* env cf? counter 0)))
    (do-emit level "}")
    0)
  (define (emit-loop level var low high code* env cf? counter f)
    (let ([v (preemit-output var env)]
	  [lo (preemit-input low env)]
	  [hi (preemit-input high env)])
      (cond
       [(null? code*) (do-emit level "~a = ~a; /* empty loop */" v hi) f]
       [else (emit-count level cf? counter f)
	     (do-emit level "for (~a = ~a; ~a < ~a; ~a++) {" v lo v hi v)
	     (emit-code* (+ level 1) code* env cf? counter 0)
	     (do-emit level "}")
	     0])))
  (define (emit-code level code env cf? counter f)
    (variant-case code
      [qa0-operation (name attr* output* input*)
        (emit-op level name attr* output* input* env f)]
      [qa0-load (type output addr*)
	(emit-load level type output addr* env f)]
      [qa0-store (type addr* value)
	(emit-store level type addr* value env f)]
      [qa0-if (var true-code* false-code*)
	(emit-if level var true-code* false-code* env cf? counter f)]
      [qa0-loop (var low high code*)
        (emit-loop level var low high code* env cf? counter f)]))
  (define (emit-code* level code* env cf? counter f)
    (cond
     [(null? code*) (emit-count level cf? counter f)]
     [else (let ([f (emit-code level (car code*) env cf? counter f)])
	     (emit-code* level (cdr code*) env cf? counter f))]))
  (define (emit-proc attr* name arg-name* arg-type* arg-c-name* arg-c-type*
		     code* env)
    (let* ([env (C-collect-args arg-name* arg-type* new-var env)]
	   [env (collect-outputs code* env)]
	   [cf? (attr-search attr* 'count-flops (lambda (v) #t) (lambda () #f))]
	   [rv (attr-search attr* 'return (lambda (v*) v*) (lambda () #f))]
	   [counter (new-var)])
      (check-inputs! code* env)
      (if (and cf? rv)
	  (error 'qa0-c99
		 "Both count-flops and return attributes in procedure ~a" name))
      (emit-proc-decl cf? rv attr* name arg-c-name* arg-c-type* env)
      (do-emit 0 "{")
      (if cf? (do-emit 1 "int ~a = 0;" counter))
      (emit-variables env)
      (emit-param* arg-name* arg-c-name* env)
      (emit-code* 1 code* env cf? counter 0)
      (if cf? (do-emit 1 "return ~a;" counter))
      (if rv (if (< (length rv) 1)
		 (error 'qa0-c99 "return without a name in ~a" name)
		 (do-emit 1 "return ~a;"
			  (preemit-output (make-reg (user-reg (car rv))) env))))
      (do-emit 0 "}~%")))
  (define (machine-c99-32 env)
    (let* ([env (machine-*-32 env c99-op* c99-load*)]
	   [env (ce-bind env 'back-end c99-back-end)]
	   [env (ce-bind env 'be-emit c99-emit)])
      env)))
	   