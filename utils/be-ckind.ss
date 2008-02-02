(module be-ckind
	mzscheme
  (require "common.ss")
  (require "ast.ss")
  (require "parser.ss")
  (require "attr.ss")
  (require "backend.ss")
  (require "cenv.ss")
  (require "cheader.ss")
  (require "verbose.ss")

  (provide build-ckind-back-end
	   empty-postparam*
	   do-emit
	   preemit-addr*
	   preemit-input
	   preemit-output
	   preemit-param
	   ck-new-var)

  (define (empty-postparam* arg-name* type* arg-c-name* c-type* env) #t)
  (define *var-count* 0)
  (define (ck-new-var)
    (let ([x (gen-reg 'g *var-count*)])
      (set! *var-count* (+ *var-count* 1))
      x))
  (define (do-emit level fmt . arg*)
    (let loop ([n level])
      (cond
       [(zero? n)]
       [else (printf "  ") (loop (- n 1))]))
    (apply printf fmt arg*)
    (newline))
  (define (preemit-output output env)
    (variant-case output
      [reg (name) (car (ce-lookup-x env 'back-end name
				    "be-ckind: name of ~a not found" name))]))
  (define (preemit-input input env)
    (variant-case input
      [reg (name) (car (ce-lookup-x env 'back-end name
				    "be-ckind: name of ~a not found" name))]
      [c-expr-number (number) number]))
  (define (preemit-addr* addr* env)
    (let loop ([r (preemit-input (car addr*) env)] [addr* (cdr addr*)])
      (cond
       [(null? addr*) r]
       [else (loop (format "~a + (~a)" r (preemit-input (car addr*) env))
		   (cdr addr*))])))
  (define (preemit-param p)
    (format "a_~a" p))
  (define (build-ckind-back-end target-name
				pointer-size
				pointer-align
				op-emit-table
				op-type-table
				ld-type-table
				emit-load
				emit-store
				collect-outputs
				the-back-end
				extra-env
				extra-decl*
				extra-def*
				extra-postparam*
				extra-undef*)
    (define (ckind-emit qa0 env)
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
  	  (emit-verbose target-name target* data*)]
	[else #t]))
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
			 "~a: register ~a is not defined" target-name name)
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
					   "~a type for ~a" target-name v)
			      (car v)))))
    (define (emit-param* arg-name* arg-c-name* env)
      (do-emit 0 "")
      (let loop ([n* arg-name*] [c* arg-c-name*])
	(cond
	 [(null? n*)]
	 [else (let* ([rn (ce-lookup-x env 'back-end (car n*)
				       "~a name for ~a" target-name (car n*))]
		      [rt (ce-lookup-x env 'name-of (cadr rn)
				       "~a type for ~a" target-name (car n*))])
		 (do-emit 1 "~a = (~a)~a;" (car rn) rt (preemit-param (car c*)))
		 (loop (cdr n*) (cdr c*)))])))
    (define (emit-count level cf? counter f)
      (if (and cf? (not (zero? f)))
	  (do-emit level "~a += ~a; /* count flops */" counter f)))
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
	 [(assq name op-emit-table)
	  => (lambda (op)
	       (let ([i-count (cadr op)] [fmt (caddr op)] [f-count (cadddr op)])
		 (do-emit-op level name output* input* out* in* i-count fmt
			     f f-count)))]
	 [else (error 'qa0-ckind "UNKNOWN op ~a, out* ~a, in* ~a, attr* ~a"
		      name output* input* attr*)])))
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
      (let* ([env (C-collect-args arg-name* arg-type* ck-new-var env)]
	     [env (collect-outputs code* env)]
	     [cf? (attr-search attr* 'count-flops (lambda (v) #t)
			       (lambda () #f))]
	     [rv (attr-search attr* 'return (lambda (v*) v*) (lambda () #f))]
	     [counter (ck-new-var)])
	(check-inputs! code* env)
	(if (and cf? rv)
	    (error 'qa0-ckind
		   "Both count-flops and return attributes in procedure ~a"
		   name))
	(emit-proc-decl cf? rv attr* name arg-c-name* arg-c-type* env)
	(do-emit 0 "{")
	(if cf? (do-emit 1 "int ~a = 0;" counter))
	(emit-variables env)
	(extra-decl* arg-name* arg-type* arg-c-name* arg-c-type* env)
	(emit-param* arg-name* arg-c-name* env)
	(extra-def* env)
	(extra-postparam* arg-name* arg-type* arg-c-name* arg-c-type* env)
	(emit-code* 1 code* env cf? counter 0)
	(extra-undef* env)
	(if cf? (do-emit 1 "return ~a;" counter))
	(if rv (if (< (length rv) 1)
		   (error 'qa0-ckind "return without a name in ~a" name)
		   (do-emit 1 "return ~a;"
			    (preemit-output (make-reg (user-reg (car rv)))
					    env))))
	(do-emit 0 "}~%")))
    (define (machine-ckind env)
      (let* ([env (machine-*-* pointer-size pointer-align env op-type-table ld-type-table)]
	     [env (extra-env env)]
	     [env (ce-bind env 'back-end the-back-end)]
	     [env (ce-bind env 'be-emit ckind-emit)])
	env))
    machine-ckind))
