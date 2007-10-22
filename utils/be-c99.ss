(module be-c99
	mzscheme
   (require "common.ss")
   (require "ast.ss")
   (require "cenv.ss")
   (require "attr.ss")
   (require "backend.ss")
   (require "parser.ss")

   (provide machine-c99-32)
   
   (define new-var
     (let ([x 0]) (lambda () (let ([s (format "g~a" x)])
			       (set! x (+ x 1))
			       s))))
   (define (c99-back-end ast env) (values ast env))
   (define (c99-emit qa0 env)
     (define (emit-indent n)
       (cond
	[(zero? n)]
	[else (printf "  ") (emit-indent (- n 1))]))
     (define (do-emit n fmt . arg*)
       (emit-indent n)
       (apply printf fmt arg*)
       (newline))
     (define (emit-decl decl)
       (define (proc-outputs arg-name* arg-type* code* env)
	 (define (collect-args arg-name* arg-type* env)
	   (cond
	    [(null? arg-name*) env]
	    [else (let* ([env (ce-bind-x env 'back-end (car arg-name*)
					 (list (new-var) (car arg-type*)))])
		    (collect-args (cdr arg-name*) (cdr arg-type*) env))]))
	 (define (collect-output* code* env)
	   (cond
	    [(null? code*) env]
	    [else
	     (collect-output* (cdr code*) (c99-output-code (car code*) env))]))
	 (define (c99-output-code code env)
	   (variant-case code
             [qa0-operation (name output*)
               (add-output* output* (be-out-type* env name) env)]
	     [qa0-load (output type)
               (add-output output (be-load-type env type) env)]
	     [qa0-store () env]
	     [qa0-loop (var code*) (let ([env (collect-output* code* env)])
				     (add-output var 'int env))]
	     [qa0-if (true-code* false-code*)
	       (let* ([env (collect-output* true-code* env)]
		      [env (collect-output* false-code* env)])
		 env)]))
	 (define (add-output* output* type* env)
	   (cond
	    [(null? output*) env]
	    [else (add-output* (cdr output*) (cdr type*)
			       (add-output (car output*) (car type*) env))]))
	 (define (add-output output type env)
	   (variant-case output
             [reg (name)
               (ce-search env (list 'back-end name)
			  (lambda (v) env)
			  (lambda () (ce-bind-x env 'back-end name
						(list (new-var) type))))]))
	 (let* ([env (collect-args arg-name* arg-type* env)]
		[env (collect-output* code* env)])
	   env))
       (define (input-code* code* env)
	 (define (input-code code env)
	   (variant-case code
             [qa0-operation (input*) (c99-input* input* env)]
	     [qa0-load (addr*) (c99-input* addr* env)]
	     [qa0-store (addr* value)
	       (c99-input* addr* env) (c99-input value env)]
	     [qa0-loop (low high code*) (input-code* code* env)
               (c99-input low env) (c99-input high env)]
	     [qa0-if (var true-code* false-code*) (c99-input var env)
               (input-code* true-code* env) (input-code* false-code* env)]))
	 (define (c99-input* input* env)
	   (for-each (lambda (input) (c99-input input env)) input*))
	 (define (c99-input input env)
	   (variant-case input
             [reg (name)
	       (ce-lookup-x env 'back-end name
			    "C99 looking for definition of ~a" name)]
	     [else #t]))
	 (for-each (lambda (code) (input-code code env)) code*))
       (variant-case decl
         [qa0-proc (attr* name
		    arg-name* arg-type* arg-c-name* arg-c-type*
		    code*)
           (let ([env (proc-outputs arg-name* arg-type* code* env)])
	     (input-code* code* env)
	     (emit-proc attr* name
			arg-name* arg-type*
			arg-c-name* arg-c-type*
			code* env))]
	 [else #t]))
     (define (emit-proc attr* name
			arg-name* arg-type*
			arg-c-name* arg-c-type*
			code* env)
       (let ([cf? (attr-search attr* 'count-flops
			       (lambda (v) #t) (lambda () #f))]
	     [counter (new-var)])
	 (emit-proc-decl cf? attr* name arg-c-name* arg-c-type* env)
	 (printf "{~%")
	 (if cf? (do-emit 1 "int ~a = 0;" counter))
	 (emit-variables env)
	 (emit-param* arg-name* arg-c-name* env)
	 (emit-code* 1 code* env cf? counter 0)
	 (if cf? (do-emit 1 "return ~a;" counter))
	 (printf "}~%~%")
	 ))
     (define (emit-proc-decl count-flops? attr* name
			     arg-c-name* arg-c-type* env)
       (printf "~a~%~a" (if count-flops? "int" "void")
	       (build-proc-name
		(attr-lookup attr* 'stem "missing stem in proc ~a" name) env))
       (if (null? arg-c-name*) (printf "(void)~%")
	   (begin (printf "(")
		  (let loop ([name* arg-c-name*] [type* arg-c-type*])
		    (cond
		     [(null? name*)]
		     [else (printf "~a~a ~a" (if (eq? name* arg-c-name*)
						 "" ", ")
				   (car type*) (car name*))
			   (loop (cdr name*) (cdr type*))]))
		  (printf ")~%"))))
     (define (emit-variables env)
       (ce-for-each env
		    (lambda (k v) (and (list? k) (eq? (car k) 'back-end)))
		    (lambda (k v)
		      (do-emit 1  "~a ~a;"
			       (ce-lookup-x env 'name-of (cadr v)
					    "c99 type for ~a" v)
			       (car v)))))
     (define (emit-param* arg-name* arg-c-name* env)
       (do-emit 0 "")
       (let loop ([name* arg-name*] [c-name* arg-c-name*])
	 (cond
	  [(null? name*) (do-emit 0 "")]
	  [else (let ([rt (ce-lookup-x env 'back-end (car name*)
				       "C99 name for ~a" (car name*))])
		  (do-emit 1 "~a = (~a)~a;"
			   (car rt)
			   (ce-lookup-x env 'name-of (cadr rt)
					"C99 type for ~a" rt)
			   (car c-name*))
		  (loop (cdr name*) (cdr c-name*)))])))
     (define (emit-code* level code* env count-flops? counter f)
       (cond
	[(null? code*) (emit-count level count-flops? counter f)]
	[else (let ([f (emit-code level (car code*) env count-flops?
				  counter f)])
		(emit-code* level (cdr code*) env count-flops? counter f))]))
     
     (define (emit-count level count-flops? counter flops)
       (if (and count-flops? (not (zero? flops)))
	   (do-emit level "~a += ~a;" counter flops)))
     (define (emit-code l code env count-flops? counter flops) ; => flops'
       (variant-case code
         [qa0-operation (name output* input*)
           (emit-op l name output* input* env flops)]
	 [qa0-load (type output addr*)
	   (emit-load l type output addr* env flops)]
	 [qa0-store (type addr* value)
           (emit-store l type addr* value env flops)]
	 [qa0-if (var true-code* false-code*)
  	   (emit-if l var true-code*
		    false-code* env count-flops? counter flops)]
	 [qa0-loop (var low high code*)
           (emit-loop l var low high code* env count-flops? counter flops)]))
     (define (emit-loop l var low high code* env count-flops? counter flops)
       (cond
	[(null? code*) (do-emit l "~a = ~a; /* empty loop */" var high) flops]
	[else (let ([v (preemit-output var env)]
		    [lo (preemit-input low env)]
		    [hi (preemit-input high env)])
		(emit-count l count-flops? counter flops)
		(do-emit l "for (~a = ~a; ~a < ~a; ~a++) {" v lo v hi v)
		(emit-code* (+ l 1) code* env count-flops? counter 0)
		(do-emit l "}")
		0)]))
     (define (emit-if l var true-code*
		      false-code* env count-flops? counter flops)
       (emit-count l count-flops? counter flops)
       (do-emit l "if (~a) {" (preemit-input var env))
       (emit-code* (+ l 1) true-code* env count-flops? counter 0)
       (if (not (null? false-code*))
	   (begin (do-emit l "} else {")
		  (emit-code* (+ l 1) false-code* env count-flops? counter 0)))
       (do-emit l "}")
       0)
     (define (emit-store l type addr* value env flops)
       (do-emit l "*(~a *)(~a) = ~a;"
		(ce-lookup-x env 'name-of type "C99 name of type ~a" type)
		(preemit-addr* addr* env)
		(preemit-input value env))
       flops)
     (define (emit-load l type output addr* env flops)
       (do-emit l "~a = *(~a *)(~a);"
		(preemit-output output env)
		(ce-lookup-x env 'name-of type "C99 name of type ~a" type)
		(preemit-addr* addr* env))
       flops)
     (define (preemit-addr* addr* env)
       (let loop ([v (preemit-input (car addr*) env)] [addr* (cdr addr*)])
	 (cond
	  [(null? addr*) v]
	  [else (loop (format "~a + ~a" v (preemit-input (car addr*) env))
		      (cdr addr*))])))
     (define (preemit-input in env)
       (variant-case in
         [reg (name) (car (ce-lookup-x env 'back-end name
				       "C99 of reg ~a at input" name))]
	 [c-expr-number (number) number]))
     (define (preemit-output out env)
       (variant-case out
         [reg (name) (car (ce-lookup-x env 'back-end name
				       "C99 of reg ~a at output" name))]))
     (define op-table
       `((pointer-move            1 "~a = ~a;"                         0)
	 (pointer-add             2 "~a = ~a + ~a;"                    0)
	 (int-move                1 "~a = ~a;"                         0)
	 (int-mul                 2 "~a = ~a * ~a;"                    0)
	 (int-div                 2 "~a = ~a / ~a;"                    0)
	 (int-mod                 2 "~a = ~a % ~a;"                    0)
	 (int-add                 2 "~a = ~a + ~a;"                    0)
	 (int-sub                 2 "~a = ~a - ~a;"                    0)
	 (int-neg                 1 "~a = -~a;"                        0)
	 (int-and                 2 "~a = ~a & ~a;"                    0)
	 (int-or                  2 "~a = ~a | ~a;"                    0)
	 (int-xor                 2 "~a = ~a ^ ~a;"                    0)
	 (int-not                 1 "~a = !~a;"                        0)
	 (double-move             1 "~a = ~a;"                         0)
	 (double-neg              1 "~a = -~a;"                        1)
	 (double-add              2 "~a = ~a + ~a;"                    1)
	 (double-sub              2 "~a = ~a - ~a;"                    1)
	 (double-div              2 "~a = ~a / ~a;"                    1)
	 (double-mul              2 "~a = ~a * ~a;"                    1)
	 (double-madd             3 "~a = ~a + ~a * ~a;"               2)
	 (complex-move            1 "~a = ~a;"                         0)
	 (commplex                2 "~a = ~a + I * ~a;"                0)
	 (complex-real            1 "~a = creal(~a);"                  0)
	 (complex-imag            1 "~a = cimag(~a)"                   0)
	 (complex-neg             1 "~a = -~a;"                        2)
	 (complex-times-plus-i    1 "~a = I * ~a;"                     1)
	 (complex-times-minus-i   1 "~a = -I * ~a;"                    1)
	 (complex-add             2 "~a = ~a + ~a;"                    2)
	 (complex-sub             2 "~a = ~a - ~a;"                    2)
	 (complex-rmul            2 "~a = ~a * ~a;"                    2)
	 (complex-mul             2 "~a = ~a * ~a;"                    6)
	 (complex-madd            3 "~a = ~a + ~a * ~a;"               8)
	 (complex-rmadd           3 "~a = ~a + ~a * ~a;"               4)
	 (complex-cmul            2 "~a = conj(~a) * ~a;"              6)
	 (complex-cmadd           2 "~a = ~a + conj(~a) * ~a;"         8)
	 (complex-add-i           2 "~a = ~a + I * ~a;"                2)
	 (complex-sub-i           2 "~a = ~a - I * ~a;"                2)
	 (complex-norm-init       0 "~a = 0.0;"                        0)
	 (complex-norm-add        2 "~a = ~a + creal(conj(~a) * ~a);"  4) ;;
	 (complex-norm-fini       1 "~a = ~a;"                         0)
	 (complex-dot-init        0 "~a = 0.0;"                        0)
	 (complex-dot-add         3 "~a = ~a + conj(~a) * ~a;"         8)
	 (complex-dot-fini        1 "~a = ~a;"                         0)))
     (define (emit-op l name output* input* env flops)
       (let ([inx* (if (eq? name 'complex-norm-add)
		       (list (car input*) (cadr input*) (cadr input*))
		       input*)])
	 (cond
	  [(assq name op-table)
	   => (lambda (op-rec)
		(let ([i-count (cadr op-rec)]
		      [fmt (caddr op-rec)]
		      [f-count (cadddr op-rec)])
		  (do-emit-op l name output* input* inx* i-count fmt
			      flops f-count)))]
	  [else (error 'c99-emit "UNKNOWN op ~a, out* ~a, in* ~a"
		       name output* input*)])))
     (define (check-* name in* out* p* size)
       (if (not (= size (length p*)))
	   (error 'c99-backend "Illformed operation: ~a <- ~a ~a"
		  out* name in*)))
     (define (do-emit-op l name out* in* inx* in-count fmt flops more-flops)
       (check-* name in* out* out* 1)
       (check-* name in* out* in* in-count)
       (apply do-emit l fmt (car out*) inx*)
       (+ flops more-flops))
     (variant-case qa0
       [qa0-top (decl*) (let loop ([decl* decl*])
			  (cond
			   [(null? decl*)]
			   [else (emit-decl (car decl*))
				 (loop (cdr decl*))]))]))
   (define (machine-c99-32 env)
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
	 (pointer-move              pointer)))
     (define c99-load*
       '((int               int)
	 (pointer           pointer)
	 (float             double)
	 (double            double)
	 (complex-float     complex-double)
	 (complex-double    complex-double)))
     (let* ([env (machine-*-32 env c99-op* c99-load*)]
	    [env (ce-bind env 'back-end c99-back-end)]
	    [env (ce-bind env 'be-emit c99-emit)])
       env)))
  
