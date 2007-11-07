(module q2complex
        mzscheme
   (require "common.ss")
   (require "ast.ss")
   (require "cenv.ss")
   (require "attr.ss")
   (provide qcd->complex)

   (define new-reg
     (let ([n 0])
       (lambda () (let ([r (string->symbol (format "~a~a" 'c n))])
		    (set! n (+ n 1))
		    r))))
   (define (q2c-decl decl env)
     (variant-case decl
       [qa0-proc (attr* name arg-name* arg-type* arg-c-name* arg-c-type* code*)
         (q2c-proc attr* name arg-name* arg-type* arg-c-name* arg-c-type*
		   code* env)]         
       [else (values decl env)]))
   (define (q2c-proc attr* name arg-name* arg-type* arg-c-name* arg-c-type*
		     code* env)
     (let-values* ([(c* env) (q2c-code* code* env)])
       (values (make-qa0-proc attr* name arg-name* arg-type*
			      arg-c-name* arg-c-type* c*)
	       env)))
   (define (q2c-code* code* env)
     (let loop ([r* '()] [code* code*] [env env])
       (cond
	[(null? code*) (values (reverse r*) env)]
	[else (let-values* ([(r* env) (q2c-code (car code*) r* env)])
		(loop r* (cdr code*) env))])))
   (define (q2c-code c r* env)
     (variant-case c
       [qa0-operation (attr* name output* input*)
	 (q2c-operation c attr* name output* input* r* env)]
       [qa0-load (attr* type output addr*)
         (q2c-load c attr* type output addr* r* env)]
       [qa0-store (attr* type addr* value)
         (q2c-store c attr* type addr* value r* env)]
       [qa0-loop (attr* var low high code*)
         (q2c-loop c attr* var low high code* r* env)]
       [qa0-if (var true-code* false-code*)
	 (q2c-if c var true-code* false-code* r* env)]))
   (define (q2c-if c var true-code* false-code* r* env)
     (let-values* ([(t* env) (q2c-code* true-code* env)]
		   [(f* env) (q2c-code* false-code* env)])
       (values (cons (make-qa0-if var t* f*) r*)
	       env)))
   (define (q2c-loop c attr* var low high code* r* env)
     (let-values* ([(c* env) (q2c-code* code* env)])
       (values (cons (make-qa0-loop attr* var low high c*) r*)
	       env)))
   (define (q2c-store-su-n attr* addr* value r* env)
     (q2c-store-xy attr* addr* value '*colors* '*colors*
		   'all 'gauge r* env))
   (define (q2c-store-fermion attr* addr* value r* env)
     (q2c-store-xy attr* addr* value '*colors* '*fermion-dim*
		   'all 'fermion r* env))
   (define (q2c-store-fermion-lo attr* addr* value r* env)
     (q2c-store-xy attr* addr* value '*colors* '*fermion-dim*
		   'low 'fermion r* env))
   (define (q2c-store-fermion-hi attr* addr* value r* env)
     (q2c-store-xy attr* addr* value '*colors* '*fermion-dim*
		   'high 'fermion r* env))
   (define (q2c-store-projected-fermion attr* addr* value r* env)
     (q2c-store-xy attr* addr* value '*colors* '*projected-fermion-dim*
		   'all 'projected-fermion r* env))
   (define (q2c-store-xy attr* addr* value c-n f-n part t r* env)
     (let* ([c-n (ce-lookup-x env 'const c-n "Color count")]
	    [f-n (ce-lookup-x env 'const f-n "Fermion size")]
	    [f-lo (if (eq? part 'high) (/ f-n 2) 0)]
	    [f-hi (if (eq? part 'low) (/ f-n 2) f-n)])
       (define (c-store c f r* env)
	 (let-values* ([(x env) (q2c-rename env value t c f)]
		       [c-s (ce-lookup-x env 'size-of 'COMPLEX "complex size")]
		       [d-s (ce-lookup-x env 'size-of 'complex-double
					 "double size")]
		       [c-t (if (= c-s d-s) 'complex-double 'complex-float)]
		       [off (* c-s (+ f (* f-n c)))])
		      (values (cons (make-qa0-store attr*
				      c-t
				      (append addr*
					      (list (make-c-expr-number off)))
				      (make-reg x))
				    r*)
			      env)))
       (let c-loop ([c 0] [r* r*] [env env])
	 (cond
	  [(= c c-n) (values r* env)]
	  [else (let f-loop ([f f-lo] [r* r*] [env env])
		  (cond
		   [(= f f-hi) (c-loop (+ c 1) r* env)]
		   [else (let-values* ([(r* e) (c-store c f r* env)])
		           (f-loop (+ f 1) r* e))]))]))))
   (define q2c-store*
     `((qcd-su-n               . ,q2c-store-su-n)
       (qcd-fermion            . ,q2c-store-fermion)
       (qcd-fermion-lo         . ,q2c-store-fermion-lo)
       (qcd-fermion-hi         . ,q2c-store-fermion-hi)
       (qcd-projected-fermion  . ,q2c-store-projected-fermion)))
   (define (q2c-store c attr* type addr* value r* env)
     (cond
      [(assq type q2c-store*)
       => (lambda (n&f) ((cdr n&f) attr* addr* value r* env))]
      [else (values (cons c r*) env)]))
   (define (q2c-load-su-n attr* output addr* r* env)
     (q2c-load-xy attr* output addr* '*colors* '*colors*
		  'all 'gauge r* env))
   (define (q2c-load-fermion attr* output addr* r* env)
     (q2c-load-xy attr* output addr* '*colors* '*fermion-dim*
		  'all 'fermion r* env))
   (define (q2c-load-fermion-lo attr* output addr* r* env)
     (q2c-load-xy attr* output addr* '*colors* '*fermion-dim*
		  'low 'fermion r* env))
   (define (q2c-load-fermion-hi attr* output addr* r* env)
     (q2c-load-xy attr* output addr* '*colors* '*fermion-dim*
		  'high 'fermion r* env))
   (define (q2c-load-projected-fermion attr* output addr* r* env)
     (q2c-load-xy attr* output addr* '*colors* '*projected-fermion-dim*
		  'all 'projected-fermion r* env))
   (define (q2c-load-xy attr* output addr* c-n f-n part t r* env)
     (let* ([c-n (ce-lookup-x env 'const c-n "Color count")]
	    [f-n (ce-lookup-x env 'const f-n "Fermion size")]
	    [f-lo (if (eq? part 'high) (/ f-n 2) 0)]
	    [f-hi (if (eq? part 'low) (/ f-n 2) f-n)])
       (define (c-load c f r* env)
	 (let-values* ([(x env) (q2c-rename env output t c f)]
		       [c-s (ce-lookup-x env 'size-of 'COMPLEX "complex size")]
		       [d-s (ce-lookup-x env 'size-of 'complex-double
					 "double size")]
		       [c-t (if (= c-s d-s) 'complex-double 'complex-float)]
		       [off (* c-s (+ f (* f-n c)))])
           (values (cons (make-qa0-load attr* c-t
					(make-reg x)
					(append addr*
						(list (make-c-expr-number
						       off))))
			 r*)
		   env)))
       (let c-loop ([c 0] [r* r*] [env env])
	 (cond
	  [(= c c-n) (values r* env)]
	  [else (let f-loop ([f f-lo] [r* r*] [env env])
		  (cond
		   [(= f f-hi) (c-loop (+ c 1) r* env)]
		   [else (let-values* ([(r* e) (c-load c f r* env)])
				      (f-loop (+ f 1) r* e))]))]))))
   (define q2c-load*
     `((qcd-su-n                . ,q2c-load-su-n)
       (qcd-projected-fermion   . ,q2c-load-projected-fermion)
       (qcd-fermion             . ,q2c-load-fermion)
       (qcd-fermion-lo          . ,q2c-load-fermion-lo)
       (qcd-fermion-hi          . ,q2c-load-fermion-hi)))
   (define (q2c-load c attr* type output addr* r* env)
     (cond
      [(assq type q2c-load*)
       => (lambda (n&f) ((cdr n&f) attr* output addr* r* env))]
      [else (values (cons c r*) env)]))
   (define (q2c-su-n-offset attr* output* input* r* env)
     (q2c-offset attr* output* input* '*colors* '*colors*
		 'gauge r* env))
   (define (q2c-fermion-offset attr* output* input* r* env)
     (q2c-offset attr* output* input* '*colors* '*fermion-dim*
		 'fermion r* env))
   (define (q2c-projected-fermion-offset attr* output* input* r* env)
     (q2c-offset attr* output* input* '*colors* '*projected-fermion-dim*
		 'projected-fermion r* env))
   (define (q2c-addu attr* output* input* r* env)
     (q2c-addx attr* output* input* '*colors*
	       'gauge r* env))
   (define (q2c-addf attr* output* input* r* env)
     (q2c-addx attr* output* input* '*fermion-dim*
	       'fermion r* env))
   (define (q2c-addh attr* output* input* r* env)
     (q2c-addx attr* output* input* '*projected-fermion-dim*
	       'projected-fermion r* env))
   (define (q2c-mulf attr* output* input* r* env)
     (q2c-mulx attr* output* input* '*fermion-dim*
	       'fermion r* env))
   (define (q2c-mulh attr* output* input* r* env)
     (q2c-mulx attr* output* input* '*projected-fermion-dim*
	       'projected-fermion r* env))
   (define (q2c-mulf-conj attr* output* input* r* env)
     (q2c-mulx-conj attr* output* input* '*fermion-dim*
		    'fermion r* env))
   (define (q2c-mulh-conj attr* output* input* r* env)
     (q2c-mulx-conj attr* output* input* '*projected-fermion-dim*
		    'projected-fermion r* env))
   (define (q2c-maddf attr* output* input* r* env)
     (q2c-maddx attr* output* input* '*fermion-dim*
		'all 'fermion r* env))
   (define (q2c-maddf-lo attr* output* input* r* env)
     (q2c-maddx attr* output* input* '*fermion-dim*
		'low 'fermion r* env))
   (define (q2c-maddf-hi attr* output* input* r* env)
     (q2c-maddx attr* output* input* '*fermion-dim*
		'high 'fermion r* env))
   (define (q2c-maddh attr* output* input* r* env)
     (q2c-maddx attr* output* input* '*projected-fermion-dim*
		'all 'projected-fermion r* env))
   (define (q2c-scaleu attr* output* input* r* env)
     (q2c-scalex attr* output* input* '*colors*
		 'all 'gauge r* env))
   (define (q2c-scalef attr* output* input* r* env)
     (q2c-scalex attr* output* input* '*fermion-dim*
		 'all 'fermion r* env))
   (define (q2c-scalef-lo attr* output* input* r* env)
     (q2c-scalex attr* output* input* '*fermion-dim*
		 'low 'fermion r* env))
   (define (q2c-scalef-hi attr* output* input* r* env)
     (q2c-scalex attr* output* input* '*fermion-dim*
		 'high 'fermion r* env))
   (define (q2c-scaleh attr* output* input* r* env)
     (q2c-scalex attr* output* input* '*projected-fermion-dim*
		 'all 'projected-fermion r* env))
   (define (q2c-fnorm-init attr* output* input* r* env)
     (q2c-check-list output* = 1 "QCD fermion norm init outputs")
     (q2c-check-list input* = 0 "QCD fermion norm init inputs")
     (values (cons (make-qa0-operation attr*
				       'complex-norm-init
				       output*
				       input*)
		   r*)
	     env))
   (define (q2c-fnorm-add attr* output* input* r* env)
     (define (complex-norm c f r* env)
       (let-values* ([(a env) (q2c-rename env (cadr input*) 'fermion c f)])
		    (values (cons (make-qa0-operation attr*
				    'complex-norm-add
				    output*
				    (list (car input*) (make-reg a)))
				  r*)
			    env)))
     (q2c-check-list output* = 1 "QCD fermion norm add outputs")
     (q2c-check-list input* = 2 "QCD fermion norm add inputs")
     (let ([c-n (ce-lookup-x env 'const '*colors* "Color count")]
	   [f-n (ce-lookup-x env 'const '*fermion-dim* "Fermion dimension")])
       (let c-loop ([c 0] [r* r*] [env env])
	 (cond
	  [(= c c-n) (values r* env)]
	  [else (let f-loop ([f 0] [r* r*] [env env])
		  (cond
		   [(= f f-n) (c-loop (+ c 1) r* env)]
		   [else (let-values* ([(r* env) (complex-norm c f r* env)])
			   (f-loop (+ f 1) r* env))]))]))))
   (define (q2c-fnorm-fini attr* output* input* r* env)
     (q2c-check-list output* = 1 "QCD fermion norm fini outputs")
     (q2c-check-list input* = 1 "QCD fermion norm fini inputs")
     (values (cons (make-qa0-operation attr* 'complex-norm-fini output* input*)
		   r*)
	     env))
   (define (q2c-fdot-init attr* output* input* r* env)
     (q2c-check-list output* = 1 "QCD fermion dot init outputs")
     (q2c-check-list input* = 0 "QCD fermion dot init inputs")
     (values (cons (make-qa0-operation attr* 'complex-dot-init output* input*)
		   r*)
	     env))
   (define (q2c-fdot-add attr* output* input* r* env)
     (define (complex-dot c f r* env)
       (let-values* ([(a env) (q2c-rename env (cadr input*) 'fermion c f)]
		     [(b env) (q2c-rename env (caddr input*) 'fermion c f)])
	 (values (cons (make-qa0-operation attr*
			 'complex-cmadd
			 output* (list (car input*) (make-reg a) (make-reg b)))
		       r*)
		 env)))
     (q2c-check-list output* = 1 "QCD fermion dot add outputs")
     (q2c-check-list input* = 3 "QCD fermion dot add inputs")
     (let ([c-n (ce-lookup-x env 'const '*colors* "Color count")]
	   [f-n (ce-lookup-x env 'const '*fermion-dim* "Fermion dimension")])
       (let c-loop ([c 0] [r* r*] [env env])
	 (cond
	  [(= c c-n) (values r* env)]
	  [else (let f-loop ([f 0] [r* r*] [env env])
		  (cond
		   [(= f f-n) (c-loop (+ c 1) r* env)]
		   [else (let-values* ([(r* env) (complex-dot c f r* env)])
		           (f-loop (+ f 1) r* env))]))]))))
   (define (q2c-fdot-fini attr* output* input* r* env)
     (q2c-check-list output* = 1 "QCD fermion dot fini outputs")
     (q2c-check-list input* = 1 "QCD fermion dot fini inputs")
     (values (cons (make-qa0-operation attr* 'complex-dot-fini output* input*)
		   r*)
	     env))
   (define (q2c-offset attr* output* input* c-n f-n t r* env)
     (let* ([c-n (ce-lookup-x env 'const c-n "Color count")]
	    [f-n (ce-lookup-x env 'const f-n "Fermion size")]
	    [d-n (ce-lookup-x env 'size-of 'COMPLEX "complex size")]
	    [c   (car input*)] [f (cadr input*)] [r (car output*)]
	    [r0  (make-reg (new-reg))] [m0 (* d-n f-n)]
	    [r1  (make-reg (new-reg))] [m1 d-n])
       (values `(,(make-qa0-operation '() 'int-add output* (list r0 r1))
		 ,(make-qa0-operation '() 'int-mul (list r1)
				      (list f (make-c-expr-number m1)))
		 ,(make-qa0-operation '() 'int-mul (list r0)
				      (list c (make-c-expr-number m0)))
		 ,@r*)
	       env)))
   (define (q2c-addx attr* output* input* f-n t r* env)
     (define (complex-add c f r* env)
       (let-values* ([(a env) (q2c-rename env (car input*) t c f)]
		     [(b env) (q2c-rename env (cadr input*) t c f)]
		     [(d env) (q2c-rename env (car output*) t c f)])
		    (values (cons (make-qa0-operation attr*
				    'complex-add
				    (list (make-reg d))
				    (list (make-reg a) (make-reg b)))
				  r*)
			    env)))
     (q2c-check-list output* = 1 "QCD add outputs")
     (q2c-check-list input* = 2 "QCD add inputs")
     (let ([c-n (ce-lookup-x env 'const '*colors* "Color count")]
	   [f-n (ce-lookup-x env 'const f-n "Field dimension")])
       (let c-loop ([c 0] [r* r*] [env env])
	 (cond
	  [(= c c-n) (values r* env)]
	  [else (let f-loop ([f 0] [r* r*] [env env])
		  (cond
		   [(= f f-n) (c-loop (+ c 1) r* env)]
		   [else (let-values* ([(r* env) (complex-add c f r* env)])
		           (f-loop (+ f 1) r* env))]))]))))
   (define (q2c-scalex attr* output* input* f-n part t r* env)
     (define (complex-scale c f r* env)
       (let-values* ([a (car input*)]
		     [(b env) (q2c-rename env (cadr input*) t c f)]
		     [(r env) (q2c-rename env (car output*) t c f)])
         (values (cons (make-qa0-operation attr*
			 'complex-rmul
			 (list (make-reg r))
			 (list a (make-reg b)))
		       r*)
		 env)))
     (q2c-check-list output* = 1 "QCD add outputs")
     (q2c-check-list input* = 2 "QCD add inputs")
     (let* ([c-n (ce-lookup-x env 'const '*colors* "Color count")]
	    [f-n (ce-lookup-x env 'const f-n "Field dimension")]
	    [f-lo (if (eq? part 'high) (/ f-n 2) 0)]
	    [f-hi (if (eq? part 'low) (/ f-n 2) f-n)])
       (let c-loop ([c 0] [r* r*] [env env])
	 (cond
	  [(= c c-n) (values r* env)]
	  [else (let f-loop ([f f-lo] [r* r*] [env env])
		  (cond
		   [(= f f-hi) (c-loop (+ c 1) r* env)]
		   [else (let-values* ([(r* env) (complex-scale c f r* env)])
		           (f-loop (+ f 1) r* env))]))]))))
   (define (q2c-madd-lohi attr* output* input* r* env)
     (let ([t 'fermion])
       (define (complex-madd c f s alpha a r* env)
	 (let-values* ([(s env) (q2c-rename env s t c f)]
		       [(a env) (q2c-rename env a t c f)]
		       [(x env) (q2c-rename env (car output*) t c f)])
	   (values (cons (make-qa0-operation attr*
			   'complex-rmadd
			   (list (make-reg x))
			   (list (make-reg s) alpha (make-reg a)))
			 r*)
		   env)))
       (q2c-check-list output* = 1 "QCD madd-lohi outputs")
       (q2c-check-list input* = 5 "QCD madd-lohi inputs")
       (let* ([c-n (ce-lookup-x env 'const '*colors* "Color count")]
	      [f-n (ce-lookup-x env 'const '*fermion-dim* "Field dimension")]
	      [f-m (/ f-n 2)]
	      [s  (car input*)]
	      [alpha (cadr input*)] [a (caddr input*)]
	      [beta (cadddr input*)] [b (car (cddddr input*))])
	 (let c-loop ([c 0] [r* r*] [env env])
	   (cond
	    [(= c c-n) (values r* env)]
	    [else
	     (let f-loop ([f 0] [r* r*] [env env])
	       (cond
		[(= f f-n) (c-loop (+ c 1) r* env)]
		[else
		 (let-values*
		     ([(r* env) (if (< f f-m)
				    (complex-madd c f s alpha a r* env)
				    (complex-madd c f s beta b r* env))])
		   (f-loop (+ f 1) r* env))]))])))))
   (define (q2c-maddx attr* output* input* f-n part t r* env)
     (define (complex-madd c f s alpha a r* env)
       (let-values* ([(s env) (q2c-rename env s t c f)]
		     [(a env) (q2c-rename env a t c f)]
		     [(x env) (q2c-rename env (car output*) t c f)])
         (values (cons (make-qa0-operation attr*
	  	         'complex-rmadd
			 (list (make-reg x))
			 (list (make-reg a) alpha (make-reg s)))
		       r*)
		 env)))
     (q2c-check-list output* = 1 "QCD madd outputs")
     (q2c-check-list input* = 3 "QCD madd inputs")
     (let* ([c-n (ce-lookup-x env 'const '*colors* "Color count")]
	    [f-n (ce-lookup-x env 'const f-n "Field dimension")]
	    [f-lo (if (eq? part 'high) (/ f-n 2) 0)]
	    [f-hi (if (eq? part 'low) (/ f-n 2) f-n)]
	    [a  (car input*)]
	    [alpha (cadr input*)]
	    [s (caddr input*)])
       (let c-loop ([c 0] [r* r*] [env env])
	 (cond
	  [(= c c-n) (values r* env)]
	  [else (let f-loop ([f f-lo] [r* r*] [env env])
		  (cond
		   [(= f f-hi) (c-loop (+ c 1) r* env)]
		   [else (let-values*
			     ([(r* env) (complex-madd c f s alpha a r* env)])
			   (f-loop (+ f 1) r* env))]))]))))
   (define (q2c-mul-g attr* output* input* f-n t op-0 op-k u-get r* env)
     (q2c-check-list output* = 1 "QCD mul outputs")
     (q2c-check-list input* = 2 "QCD mul inputs")
     (let ([c-n (ce-lookup-x env 'const '*colors* "Color count")]
	   [f-n (ce-lookup-x env 'const f-n "Field dimension")]
	   [r-a (car input*)]
	   [r-b (cadr input*)]
	   [r-r (car output*)])
       (define (s-mul-z r* env)
	 (let-values* ([(a env) (u-get r-a 0 0 env)])
	   (let loop ([f 0] [r* r*] [env env])
	     (cond
	      [(= f f-n) (values r* env)]
	      [else (let-values* ([(b env) (q2c-rename env r-b t 0 f)]
				  [(z env) (q2c-rename env r-r t 0 f)])
	              (loop (+ f 1)
			    (cons (make-qa0-operation attr*
				    op-0
				    (list (make-reg z))
				    (list (make-reg a) (make-reg b)))
				  r*)
			    env))]))))
       (define (s-mul-1 r* env)
	 (let ([q (new-reg)])
	   (let x-loop ([x 0] [r* r*] [env env])
	     (cond
	      [(= x c-n) (values r* env q)]
	      [else
	       (let y-loop ([y 0] [r* r*] [env env])
		 (cond
		  [(= y f-n) (x-loop (+ x 1) r* env)]
		  [else
		   (let-values* ([(q-v env) (q2c-rename env q t x y)]
				 [(a-v env) (u-get r-a x 0 env)]
				 [(b-v env) (q2c-rename env r-b t 0 y)])
                     (y-loop (+ y 1)
			     (cons (make-qa0-operation attr*
				     op-0
				     (list (make-reg q-v))
				     (list (make-reg a-v) (make-reg b-v)))
				   r*)
			     env))]))]))))
       (define (s-madd-x r-x r-r c r* env)
	 (let x-loop ([x 0] [r* r*] [env env])
	   (cond
	    [(= x c-n) (values r* env r-r)]
	    [else
	     (let y-loop ([y 0] [r* r*] [env env])
	       (cond
		[(= y f-n) (x-loop (+ x 1) r* env)]
		[else
		 (let-values* ([(q-v env) (q2c-rename env r-x t x y)]
			       [(r-v env) (q2c-rename env r-r t x y)]
			       [(a-v env) (u-get r-a x c env)]
			       [(b-v env) (q2c-rename env r-b t c y)])
		   (y-loop (+ y 1)
			   (cons (make-qa0-operation attr*
                                   op-k
				   (list (make-reg r-v))
				   (list (make-reg a-v)
					 (make-reg b-v)
					 (make-reg q-v)))
				 r*)
			   env))]))])))
       (if (= c-n 1) (s-mul-z r* env)
	   (let-values* ([(r* env r-x) (s-mul-1 r* env)])
	     (let loop ([c 1] [r* r*] [env env] [r-x r-x])
	       (cond
		[(= c (- c-n 1))
		  (let-values* ([(r* env r-x) (s-madd-x r-x r-r c r* env)])
		    (values r* env))]
		[else
                  (let-values* ([(r* env r-x) (s-madd-x r-x (new-reg) c
							r* env)])
		    (loop (+ c 1) r* env r-x))]))))))
   (define (q2c-mulx attr* output* input* f-n t r* env)
     (q2c-mul-g attr* output* input* f-n t
		'complex-mul 'complex-madd
		(lambda (r i j env) (q2c-rename env r 'gauge i j))
		r* env))
   (define (q2c-mulx-conj attr* output* input* f-n t r* env)
     (q2c-mul-g attr* output* input* f-n t
		'complex-cmul 'complex-cmadd
		(lambda (r i j env) (q2c-rename env r 'gauge j i))
		r* env))
   (define (q2c-check-list x* op size msg)
     (if (not (op (length x*) size))
	 (error 'qcd->complex "ERROR: ~a" msg)))
   (define (q2c-rename env base type i-a i-b)
     (let ([key (list 'qcd->complex base type i-a i-b)])
       (ce-search env key
		  (lambda (val)
		    (values val env))
		  (lambda ()
		    (let* ([r (new-reg)]
			   [env (ce-bind env key r)])
		      (values r env))))))
   (define (q2c-project attr* output* input* r* env)
     (q2c-check-list output* = 1 "QCD gamma projection result")
     (q2c-check-list input* = 1 "QCD gamma projection source")
     (let* ([kind (attr-lookup attr* 'project "qcd-project")]
	    [op* (ce-lookup env `(project ,@kind)
			    "project op-table for ~a" kind)]
	    [c-n (ce-lookup-x env 'const '*colors* "Color count")]
	    [h-n (/ (ce-lookup-x env 'const '*fermion-dim* "Fermion dim") 2)]
	    [r-r (car output*)]
	    [r-a (car input*)])
       (define (proj c h op r* env)
	 (q2c-check-list op = 4 "Projection operation")
	 (let-values* ([op-0 (car op)] [(f-0) (cadr op)]
		       [(a-0 env) (q2c-rename env r-a 'fermion c f-0)]
		       [op-1 (caddr op)] [(f-1) (cadddr op)]
		       [(a-1 env) (q2c-rename env r-a 'fermion c f-1)]
		       [cmd (binary-cmd op-0 op-1)]
		       [(r env) (q2c-rename env r-r 'projected-fermion c h)])
           (values (cons (make-qa0-operation attr*
                           cmd
			   (list (make-reg r))
			   (list (make-reg a-0) (make-reg a-1)))
			 r*)
		   env)))
       (define (binary-cmd op-0 op-1)
	 (case op-0
	   [(plus-one)
	    (case op-1
	      [(plus-one) 'complex-add]
	      [(minus-one) 'complex-sub]
	      [(plus-i) 'complex-add-i]
	      [(minus-i) 'complex-sub-i]
	      [else (error 'qa0 "Unknown second factor (~a ~a)" op-0 op-1)])]
	   [else (error 'qa0 "Unknown first factor (~a ~a)" op-0 op-1)]))
       (let c-loop ([c 0] [r* r*] [env env])
	 (cond
	  [(= c c-n) (values r* env)]
	  [else (let h-loop ([h 0] [op* op*] [r* r*] [env env])
		  (cond
		   [(null? op*) (c-loop (+ c 1) r* env)]
		   [else (let-values* ([(r* env) (proj c h (car op*) r* env)])
			   (h-loop (+ h 1) (cdr op*) r* env))]))]))))
   (define (q2c-unproject attr* output* input* r* env)
     (q2c-check-list output* = 1 "QCD gamma unprojection result")
     (q2c-check-list input* = 1 "QCD gamma unprojection source")
     (let* ([kind (attr-lookup attr* 'unproject "qcd-unproject")]
	    [op* (ce-lookup env `(unproject ,@kind)
			    "unproj op-table for ~a" kind)]
	    [c-n (ce-lookup-x env 'const '*colors* "Color count")]
	    [f-n (ce-lookup-x env 'const '*fermion-dim* "Fermion dim")]
	    [r-r (car output*)]
	    [r-a (car input*)])
       (define (unproj c f op r* env)
	 (q2c-check-list op = 2 "Unprojection operation")
	 (let-values* ([opcode (car op)] [(component) (cadr op)]
		       [(a env) (q2c-rename env r-a
					    'projected-fermion c component)]
		       [cmd (unary-cmd opcode)]
		       [(r env) (q2c-rename env r-r 'fermion c f)])
	   (values (cons (make-qa0-operation attr*
                           cmd (list (make-reg r)) (list (make-reg a)))
			 r*)
		   env)))
       (define (unary-cmd name)
	 (case name
	   [(plus-one)   'complex-move]
	   [(minus-one)  'complex-neg]
	   [(plus-i)     'complex-times-plus-i]
	   [(minus-i)    'complex-times-minus-i]
	   [else (error 'qa0 "Unknown unproject operation ~a" name)]))
       (let c-loop ([c 0] [r* r*] [env env])
	 (cond
	  [(= c c-n) (values r* env)]
	  [else (let f-loop ([f 0] [op* op*] [r* r*] [env env])
		  (cond
		   [(null? op*) (c-loop (+ c 1) r* env)]
		   [else (let-values* ([(r* env) (unproj c f (car op*) r* env)])
		           (f-loop (+ f 1) (cdr op*) r* env))]))]))))

   (define (q2c-unproject-add attr* output* input* r* env)
     (q2c-check-list output* = 1 "QCD gamma unproject-add result")
     (q2c-check-list input* = 2 "QCD gamma unproject-add sources")
     (let* ([kind (attr-lookup attr* 'unproject "qcd-unproject-add")]
	    [op* (ce-lookup env `(unproject ,@kind)
			    "unproj op-table for ~a" kind)]
	    [c-n (ce-lookup-x env 'const '*colors* "Color count")]
	    [r-r (car output*)]
	    [r-a (car input*)]
	    [r-b (cadr input*)])
       (define (upa c f op r* env)
	 (q2c-check-list op = 2 "Unprojection operation")
	 (let-values* ([(a env) (q2c-rename env r-a 'fermion c f)]
		       [opcode (car op)] [(component) (cadr op)]
		       [(b env) (q2c-rename env r-b 'projected-fermion
					    c component)]
		       [cmd (add-cmd opcode)]
		       [(r env) (q2c-rename env r-r 'fermion c f)])
	   (values (cons (make-qa0-operation attr*
                           cmd
			   (list (make-reg r))
			   (list (make-reg a) (make-reg b)))
			 r*)
		   env)))
       (define (add-cmd name)
	 (case name
	   [(plus-one)   'complex-add]
	   [(minus-one)  'complex-sub]
	   [(plus-i)     'complex-add-i]
	   [(minus-i)    'complex-sub-i]
	   [else (error 'qa0 "Unknown unproject operation ~a" name)]))
       (let c-loop ([c 0] [r* r*] [env env])
	 (cond
	  [(= c c-n) (values r* env)]
	  [else (let f-loop ([f 0] [op* op*] [r* r*] [env env])
		  (cond
		   [(null? op*) (c-loop (+ c 1) r* env)]
		   [else (let-values* ([(r* env) (upa c f (car op*) r* env)])
			   (f-loop (+ f 1) (cdr op*) r* env))]))]))))
   (define q2c-op*
     `((qcd-project                  . ,q2c-project)
       (qcd-unproject                . ,q2c-unproject)
       (qcd-unproject-add            . ,q2c-unproject-add)
       (qcd-mulf                     . ,q2c-mulf)
       (qcd-mulh                     . ,q2c-mulh)
       (qcd-mulf-conj                . ,q2c-mulf-conj)
       (qcd-mulh-conj                . ,q2c-mulh-conj)
       (qcd-scaleu                   . ,q2c-scaleu)
       (qcd-scalef                   . ,q2c-scalef)
       (qcd-scalef-lo                . ,q2c-scalef-lo)
       (qcd-scalef-hi                . ,q2c-scalef-hi)
       (qcd-scaleh                   . ,q2c-scaleh)
       (qcd-addu                     . ,q2c-addu)
       (qcd-addf                     . ,q2c-addf)
       (qcd-addh                     . ,q2c-addh)
       (qcd-maddf                    . ,q2c-maddf)
       (qcd-maddf-lo                 . ,q2c-maddf-lo)
       (qcd-maddf-hi                 . ,q2c-maddf-hi)
       (qcd-maddh                    . ,q2c-maddh)
       (qcd-madd-lohi                . ,q2c-madd-lohi)
       (qcd-fnorm-init               . ,q2c-fnorm-init)
       (qcd-fnorm-add                . ,q2c-fnorm-add)
       (qcd-fnorm-fini               . ,q2c-fnorm-fini)
       (qcd-fdot-init                . ,q2c-fdot-init)
       (qcd-fdot-add                 . ,q2c-fdot-add)
       (qcd-fdot-fini                . ,q2c-fdot-fini)
       (qcd-su-n-offset              . ,q2c-su-n-offset)
       (qcd-fermion-offset           . ,q2c-fermion-offset)
       (qcd-projected-fermion-offset . ,q2c-projected-fermion-offset)))
   (define (q2c-operation c attr* name output* input* r* env)
     (cond
      [(assq name q2c-op*)
       => (lambda (n&f) ((cdr n&f) attr* output* input* r* env))]
      [else (values (cons c r*) env)]))
   (define (qcd->complex ast env)
     (variant-case ast
       [qa0-top (decl*)
	 (let loop ([r* '()] [decl* decl*] [env env])
	   (cond
	    [(null? decl*) (values (make-qa0-top (reverse r*)) env)]
	    [else (let-values* ([(d e) (q2c-decl (car decl*) env)])
		    (loop (cons d r*) (cdr decl*) e))]))])))

