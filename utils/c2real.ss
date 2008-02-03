(module c2real
        mzscheme
   (require "common.ss")
   (require "ast.ss")
   (require "cenv.ss")
   (require "attr.ss")
   (provide complex->real)

   (define new-reg
     (let ([n 0])
       (lambda () (let ([r (gen-reg 'r n)])
		    (set! n (+ n 1))
		    r))))
   (define *reg-count* 0)
   (define (fresh-reg)
     (let ([x (gen-reg 'cr *reg-count*)])
       (set! *reg-count* (+ *reg-count* 1))
       (make-reg x)))
   (define (c2r-rename env base part)
     (let ([key (list 'complex->real base part)])
       (ce-search env key
		  (lambda (val)
		    (values val env))
		  (lambda ()
		    (let* ([r (new-reg)]
			   [env (ce-bind env key r)])
		      (values r env))))))
   (define (c2r-split c2r-xxx attr* output* input* r* env)
     (let-values* ([g* (list (fresh-reg))]
		   [(r* env) (c2r-xxx attr* g* input* r* env)])
       (c2r-complex-move attr* output* g* r* env)))
   (define (c2r-decl decl env)
     (variant-case decl
       [qa0-proc (attr* name arg-name* arg-type* arg-c-name* arg-c-type* code*)
         (c2r-proc attr* name arg-name* arg-type* arg-c-name* arg-c-type*
		   code* env)]         
       [else (values decl env)]))
   (define (c2r-proc attr* name arg-name* arg-type* arg-c-name* arg-c-type*
		     code* env)
     (let-values* ([(c* env) (c2r-code* code* env)])
       (values (make-qa0-proc attr* name arg-name* arg-type*
			      arg-c-name* arg-c-type* c*)
	       env)))
   (define (c2r-code* code* env)
     (let loop ([r* '()] [code* code*] [env env])
       (cond
	[(null? code*) (values (reverse r*) env)]
	[else (let-values* ([(r* env) (c2r-code (car code*) r* env)])
		(loop r* (cdr code*) env))])))
   (define (c2r-code c r* env)
     (variant-case c
       [qa0-operation (attr* name output* input*)
	 (c2r-operation c attr* name output* input* r* env)]
       [qa0-load (attr* type output addr*)
         (c2r-load c attr* type output addr* r* env)]
       [qa0-store (attr* type addr* value)
         (c2r-store c attr* type addr* value r* env)]
       [qa0-loop (attr* var low high code*)
         (c2r-loop c attr* var low high code* r* env)]
       [qa0-if (var true-code* false-code*)
	 (c2r-if c var true-code* false-code* r* env)]))
   (define (c2r-if c var true-code* false-code* r* env)
     (let-values* ([(t* env) (c2r-code* true-code* env)]
		   [(f* env) (c2r-code* false-code* env)])
       (values (cons (make-qa0-if var t* f*) r*)
	       env)))
   (define (c2r-loop c attr* var low high code* r* env)
     (let-values* ([(c* env) (c2r-code* code* env)])
       (values (cons (make-qa0-loop attr* var low high c*) r*)
	       env)))
   (define (c2r-store-x attr* addr* value r-type msg r* env)
     (let-values* ([m-s (ce-lookup-x env 'size-of r-type msg)]
		   [off-im (list (make-c-expr-number m-s))]
		   [(x-re env) (c2r-rename env value 'real)]
		   [(x-im env) (c2r-rename env value 'imag)])
       (values `(,(make-qa0-store attr* r-type (append addr* off-im)
				  (make-reg x-im))
		 ,(make-qa0-store attr* r-type addr*
				  (make-reg x-re))
		 ,@r*)
	       env)))
   (define (c2r-store-complex attr* type addr* value r* env)
     (c2r-store-x attr* addr* value 'REAL "REAL size" r* env))
   (define (c2r-store-complex-double attr* type addr* value r* env)
     (c2r-store-x attr* addr* value 'double "double size" r* env))
   (define (c2r-store-complex-float attr* type addr* value r* env)
     (c2r-store-x attr* addr* value 'float "float size" r* env))
   (define c2r-store*
     `((COMPLEX        . ,c2r-store-complex)
       (complex-double . ,c2r-store-complex-double)
       (complex-float  . ,c2r-store-complex-float)))
   (define (c2r-store c attr* type addr* value r* env)
     (cond
      [(assq type c2r-store*)
       => (lambda (n&f) ((cdr n&f) attr* type addr* value r* env))]
      [else (values (cons c r*) env)]))
   (define (c2r-load-x attr* type output addr* real-type msg r* env)
     (let-values* ([m-s (ce-lookup-x env 'size-of real-type msg)]
		   [off-im (list (make-c-expr-number m-s))]
		   [(x-re env) (c2r-rename env output 'real)]
		   [(x-im env) (c2r-rename env output 'imag)])
       (values `(,(make-qa0-load attr* real-type (make-reg x-im)
				 (append addr* off-im))
		 ,(make-qa0-load attr* real-type (make-reg x-re)
				 addr*)
		 ,@r*)
	       env)))
   (define (c2r-load-complex attr* type output addr* r* env)
     (c2r-load-x attr* type output addr* 'REAL "REAL size" r* env))
   (define (c2r-load-complex-double attr* type output addr* r* env)
     (c2r-load-x attr* type output addr* 'double "double size" r* env))
   (define (c2r-load-complex-float attr* type output addr* r* env)
     (c2r-load-x attr* type output addr* 'float "float size" r* env))
   (define c2r-load*
     `((COMPLEX        . ,c2r-load-complex)
       (complex-double . ,c2r-load-complex-double)
       (complex-float  . ,c2r-load-complex-float)))
   (define (c2r-load c attr* type output addr* r* env)
     (cond
      [(assq type c2r-load*)
       => (lambda (n&f) ((cdr n&f) attr* type output addr* r* env))]
      [else (values (cons c r*) env)]))
   (define (c2r-complex-move attr* output* input* r* env)
     (let-values* ([(a-re env) (c2r-rename env (car input*) 'real)]
		   [(a-im env) (c2r-rename env (car input*) 'imag)]
		   [(b-re env) (c2r-rename env (car output*) 'real)]
		   [(b-im env) (c2r-rename env (car output*) 'imag)])
       (values `(,(make-qa0-operation attr* 'double-move
				      (list (make-reg b-im))
				      (list (make-reg a-im)))
		 ,(make-qa0-operation attr* 'double-move
				      (list (make-reg b-re))
				      (list (make-reg a-re)))
		 ,@r*)
	       env)))
   (define (c2r-complex attr* output* input* r* env)
     (let-values* ([(b-re env) (c2r-rename env (car output*) 'real)]
		   [(b-im env) (c2r-rename env (car output*) 'imag)])
       (values `(,(make-qa0-operation attr* 'double-move
				      (list (make-reg b-im))
				      (list (cadr input*)))
		 ,(make-qa0-operation attr* 'double-move
				      (list (make-reg b-re))
				      (list (car input*)))
		 ,@r*)
	       env)))
   (define (c2r-complex-real attr* output* input* r* env)
     (let-values* ([(a-re env) (c2r-rename env (car input*) 'real)])
       (values `(,(make-qa0-operation attr* 'double-move
				      output*
				      (list (make-reg a-re)))
		 ,@r*)
	       env)))
   (define (c2r-complex-imag attr* output* input* r* env)
     (let-values* ([(a-im env) (c2r-rename env (car input*) 'imag)])
       (values `(,(make-qa0-operation attr* 'double-move
				      output*
				      (list (make-reg a-im)))
		 ,@r*)
	       env)))
   (define (c2r-complex-neg attr* output* input* r* env)
     (let-values* ([(a-re env) (c2r-rename env (car input*) 'real)]
		   [(a-im env) (c2r-rename env (car input*) 'imag)]
		   [(b-re env) (c2r-rename env (car output*) 'real)]
		   [(b-im env) (c2r-rename env (car output*) 'imag)])
       (values `(,(make-qa0-operation attr* 'double-neg
				      (list (make-reg b-im))
				      (list (make-reg a-im)))
		 ,(make-qa0-operation attr* 'double-neg
				      (list (make-reg b-re))
				      (list (make-reg a-re)))
		 ,@r*)
	       env)))
   (define (c2r-complex-add attr* output* input* r* env)
     (let-values* ([(a-re env) (c2r-rename env (car input*) 'real)]
		   [(a-im env) (c2r-rename env (car input*) 'imag)]
		   [(b-re env) (c2r-rename env (cadr input*) 'real)]
		   [(b-im env) (c2r-rename env (cadr input*) 'imag)]
		   [(c-re env) (c2r-rename env (car output*) 'real)]
		   [(c-im env) (c2r-rename env (car output*) 'imag)])
       (values `(,(make-qa0-operation attr* 'double-add
				      (list (make-reg c-im))
				      (list (make-reg a-im) (make-reg b-im)))
		 ,(make-qa0-operation attr* 'double-add
				      (list (make-reg c-re))
				      (list (make-reg a-re) (make-reg b-re)))
		 ,@r*)
	       env)))
   (define (c2r-complex-sub attr* output* input* r* env)
     (let-values* ([(a-re env) (c2r-rename env (car input*) 'real)]
		   [(a-im env) (c2r-rename env (car input*) 'imag)]
		   [(b-re env) (c2r-rename env (cadr input*) 'real)]
		   [(b-im env) (c2r-rename env (cadr input*) 'imag)]
		   [(c-re env) (c2r-rename env (car output*) 'real)]
		   [(c-im env) (c2r-rename env (car output*) 'imag)])
       (values `(,(make-qa0-operation attr* 'double-sub
				      (list (make-reg c-im))
				      (list (make-reg a-im) (make-reg b-im)))
		 ,(make-qa0-operation attr* 'double-sub
				      (list (make-reg c-re))
				      (list (make-reg a-re) (make-reg b-re)))
		 ,@r*)
	       env)))
   (define (c2r-complex-rmul attr* output* input* r* env)
     (let-values* ([a (car input*)]
		   [(b-re env) (c2r-rename env (cadr input*) 'real)]
		   [(b-im env) (c2r-rename env (cadr input*) 'imag)]
		   [(c-re env) (c2r-rename env (car output*) 'real)]
		   [(c-im env) (c2r-rename env (car output*) 'imag)])
       (values `(,(make-qa0-operation attr* 'double-mul
				      (list (make-reg c-im))
				      (list a (make-reg b-im)))
		 ,(make-qa0-operation attr* 'double-mul
				      (list (make-reg c-re))
				      (list a (make-reg b-re)))
		 ,@r*)
	       env)))
   (define (c2r-complex-rmadd attr* output* input* r* env)
     (let-values* ([(a-re env) (c2r-rename env (car input*) 'real)]
		   [(a-im env) (c2r-rename env (car input*) 'imag)]
		   [b (cadr input*)]
		   [(c-re env) (c2r-rename env (caddr input*) 'real)]
		   [(c-im env) (c2r-rename env (caddr input*) 'imag)]
		   [(d-re env) (c2r-rename env (car output*) 'real)]
		   [(d-im env) (c2r-rename env (car output*) 'imag)])
       (values `(,(make-qa0-operation attr* 'double-madd
				      (list (make-reg d-re))
				      (list (make-reg a-re) b (make-reg c-re)))
		 ,(make-qa0-operation attr* 'double-madd
				      (list (make-reg d-im))
				      (list (make-reg a-im) b (make-reg c-im)))
		 ,@r*)
	       env)))
   (define (c2r-complex-rmsub attr* output* input* r* env)
     (let-values* ([(a-re env) (c2r-rename env (car input*) 'real)]
		   [(a-im env) (c2r-rename env (car input*) 'imag)]
		   [b (cadr input*)]
		   [(c-re env) (c2r-rename env (caddr input*) 'real)]
		   [(c-im env) (c2r-rename env (caddr input*) 'imag)]
		   [(d-re env) (c2r-rename env (car output*) 'real)]
		   [(d-im env) (c2r-rename env (car output*) 'imag)])
       (values `(,(make-qa0-operation attr* 'double-msub
				      (list (make-reg d-im))
				      (list (make-reg a-im) b (make-reg c-im)))
		 ,(make-qa0-operation attr* 'double-msub
				      (list (make-reg d-re))
				      (list (make-reg a-re) b (make-reg c-re)))
		 ,@r*)
	       env)))
   (define (c2r-complex-norm-init attr* output* input* r* env)
     (values `(,(make-qa0-operation attr* 'double-zero output* '()) ,@r*)
	     env))
   (define (c2r-complex-norm-add attr* output* input* r* env)
     (let-values* ([a (car input*)]
		   [(b-re env) (c2r-rename env (cadr input*) 'real)]
		   [(b-im env) (c2r-rename env (cadr input*) 'imag)])
       (values `(,(make-qa0-operation attr* 'double-madd
				      output*
				      (list a (make-reg b-im) (make-reg b-im)))
		 ,(make-qa0-operation attr* 'double-madd
				      output*
				      (list a (make-reg b-re) (make-reg b-re)))
		 ,@r*)
	       env)))
   (define (c2r-complex-norm-fini attr* output* input* r* env)
     (values `(,(make-qa0-operation attr* 'double-move output* input*) ,@r*)
	     env))
   (define (c2r-complex-times-plus-i attr* output* input* r* env)
     (if (member (car output*) input*)
	 (c2r-split c2r-complex-times-plus-i attr* output* input* r* env)
	 (let-values* ([(a-re env) (c2r-rename env (car output*) 'real)]
		       [(a-im env) (c2r-rename env (car output*) 'imag)]
		       [(b-re env) (c2r-rename env (car input*) 'real)]
		       [(b-im env) (c2r-rename env (car input*) 'imag)])
	   (values `(,(make-qa0-operation attr* 'double-move
					  (list (make-reg a-im))
					  (list (make-reg b-re)))
		     ,(make-qa0-operation attr* 'double-neg
					  (list (make-reg a-re))
					  (list (make-reg b-im)))
		     ,@r*)
		   env))))
   (define (c2r-complex-times-minus-i attr* output* input* r* env)
     (if (member (car output*) input*)
	 (c2r-split c2r-complex-times-plus-i attr* output* input* r* env)
	 (let-values* ([(a-re env) (c2r-rename env (car output*) 'real)]
		       [(a-im env) (c2r-rename env (car output*) 'imag)]
		       [(b-re env) (c2r-rename env (car input*) 'real)]
		       [(b-im env) (c2r-rename env (car input*) 'imag)])
	   (values `(,(make-qa0-operation attr* 'double-move
					  (list (make-reg a-re))
					  (list (make-reg b-im)))
		     ,(make-qa0-operation attr* 'double-neg
					  (list (make-reg a-im))
					  (list (make-reg b-re)))
		     ,@r*)
		   env))))
   (define (c2r-complex-mul attr* output* input* r* env)
     (if (member (car output*) input*)
	 (c2r-split c2r-complex-mul attr* output* input* r* env)
	 (let-values* ([(a-re env) (c2r-rename env (car output*) 'real)]
		       [(a-im env) (c2r-rename env (car output*) 'imag)]
		       [(b-re env) (c2r-rename env (car input*) 'real)]
		       [(b-im env) (c2r-rename env (car input*) 'imag)]
		       [(c-re env) (c2r-rename env (cadr input*) 'real)]
		       [(c-im env) (c2r-rename env (cadr input*) 'imag)])
	   (values `(,(make-qa0-operation attr* 'double-madd
					  (list (make-reg a-im))
					  (list (make-reg a-im)
						(make-reg b-im)
						(make-reg c-re)))
		     ,(make-qa0-operation attr* 'double-msub
					  (list (make-reg a-re))
					  (list (make-reg a-re)
						(make-reg b-im)
						(make-reg c-im)))
		     ,(make-qa0-operation attr* 'double-mul
					  (list (make-reg a-im))
					  (list (make-reg b-re)
						(make-reg c-im)))
		     ,(make-qa0-operation attr* 'double-mul
					  (list (make-reg a-re))
					  (list (make-reg b-re)
						(make-reg c-re)))
		     ,@r*)
		   env))))
   (define (c2r-complex-cmul attr* output* input* r* env)
     (if (member (car output*) input*)
	 (c2r-split c2r-complex-cmul attr* output* input* r* env)
	 (let-values* ([(a-re env) (c2r-rename env (car output*) 'real)]
		       [(a-im env) (c2r-rename env (car output*) 'imag)]
		       [(b-re env) (c2r-rename env (car input*) 'real)]
		       [(b-im env) (c2r-rename env (car input*) 'imag)]
		       [(c-re env) (c2r-rename env (cadr input*) 'real)]
		       [(c-im env) (c2r-rename env (cadr input*) 'imag)])
	   (values `(,(make-qa0-operation attr* 'double-msub
					  (list (make-reg a-im))
					  (list (make-reg a-im)
						(make-reg b-im)
						(make-reg c-re)))
		     ,(make-qa0-operation attr* 'double-madd
					  (list (make-reg a-re))
					  (list (make-reg a-re)
						(make-reg b-im)
						(make-reg c-im)))
		     ,(make-qa0-operation attr* 'double-mul
					  (list (make-reg a-im))
					  (list (make-reg b-re)
						(make-reg c-im)))
		     ,(make-qa0-operation attr* 'double-mul
					  (list (make-reg a-re))
					  (list (make-reg b-re)
						(make-reg c-re)))
		     ,@r*)
		   env))))
   (define (c2r-complex-madd attr* output* input* r* env)
     (if (member (car output*) (cdr input*))
	 (c2r-split c2r-complex-madd attr* output* input* r* env)
	 (let-values* ([(a-re env) (c2r-rename env (car output*) 'real)]
		       [(a-im env) (c2r-rename env (car output*) 'imag)]
		       [(x-re env) (c2r-rename env (car input*) 'real)]
		       [(x-im env) (c2r-rename env (car input*) 'imag)]
		       [(b-re env) (c2r-rename env (cadr input*) 'real)]
		       [(b-im env) (c2r-rename env (cadr input*) 'imag)]
		       [(c-re env) (c2r-rename env (caddr input*) 'real)]
		       [(c-im env) (c2r-rename env (caddr input*) 'imag)])
	   (values `(,(make-qa0-operation attr* 'double-madd
					  (list (make-reg a-im))
					  (list (make-reg a-im)
						(make-reg b-im)
						(make-reg c-re)))
		     ,(make-qa0-operation attr* 'double-msub
					  (list (make-reg a-re))
					  (list (make-reg a-re)
						(make-reg b-im)
						(make-reg c-im)))
		     ,(make-qa0-operation attr* 'double-madd
					  (list (make-reg a-im))
					  (list (make-reg x-im)
						(make-reg b-re)
						(make-reg c-im)))
		     ,(make-qa0-operation attr* 'double-madd
					  (list (make-reg a-re))
					  (list (make-reg x-re)
					        (make-reg b-re)
						(make-reg c-re)))
		     ,@r*)
		   env))))
   (define (c2r-complex-cmadd attr* output* input* r* env)
     (if (member (car output*) (cdr input*))
	 (c2r-split c2r-complex-cmadd attr* output* input* r* env)
	 (let-values* ([(a-re env) (c2r-rename env (car output*) 'real)]
		       [(a-im env) (c2r-rename env (car output*) 'imag)]
		       [(x-re env) (c2r-rename env (car input*) 'real)]
		       [(x-im env) (c2r-rename env (car input*) 'imag)]
		       [(b-re env) (c2r-rename env (cadr input*) 'real)]
		       [(b-im env) (c2r-rename env (cadr input*) 'imag)]
		       [(c-re env) (c2r-rename env (caddr input*) 'real)]
		       [(c-im env) (c2r-rename env (caddr input*) 'imag)])
	   (values `(,(make-qa0-operation attr* 'double-msub
					  (list (make-reg a-im))
					  (list (make-reg a-im)
						(make-reg b-im)
						(make-reg c-re)))
		     ,(make-qa0-operation attr* 'double-madd
					  (list (make-reg a-re))
					  (list (make-reg a-re)
						(make-reg b-im)
						(make-reg c-im)))
		     ,(make-qa0-operation attr* 'double-madd
					  (list (make-reg a-im))
					  (list (make-reg x-im)
						(make-reg b-re)
						(make-reg c-im)))
		     ,(make-qa0-operation attr* 'double-madd
					  (list (make-reg a-re))
					  (list (make-reg x-re)
					        (make-reg b-re)
						(make-reg c-re)))
		     ,@r*)
		   env))))
   (define (c2r-complex-add-i attr* output* input* r* env)
     (if (member (car output*) (cdr input*))
	 (c2r-split c2r-complex-add-i attr* output* input* r* env)
	 (let-values* ([(a-re env) (c2r-rename env (car output*) 'real)]
		       [(a-im env) (c2r-rename env (car output*) 'imag)]
		       [(b-re env) (c2r-rename env (car input*) 'real)]
		       [(b-im env) (c2r-rename env (car input*) 'imag)]
		       [(c-re env) (c2r-rename env (cadr input*) 'real)]
		       [(c-im env) (c2r-rename env (cadr input*) 'imag)])
	   (values `(,(make-qa0-operation attr* 'double-add
					  (list (make-reg a-im))
					  (list (make-reg b-im)
						(make-reg c-re)))
		     ,(make-qa0-operation attr* 'double-sub
					  (list (make-reg a-re))
					  (list (make-reg b-re)
						(make-reg c-im)))
		     ,@r*)
		   env))))
   (define (c2r-complex-sub-i attr* output* input* r* env)
     (if (member (car output*) (cdr input*))
	 (c2r-split c2r-complex-sub-i attr* output* input* r* env)
	 (let-values* ([(a-re env) (c2r-rename env (car output*) 'real)]
		       [(a-im env) (c2r-rename env (car output*) 'imag)]
		       [(b-re env) (c2r-rename env (car input*) 'real)]
		       [(b-im env) (c2r-rename env (car input*) 'imag)]
		       [(c-re env) (c2r-rename env (cadr input*) 'real)]
		       [(c-im env) (c2r-rename env (cadr input*) 'imag)])
	   (values `(,(make-qa0-operation attr* 'double-sub
					  (list (make-reg a-im))
					  (list (make-reg b-im)
						(make-reg c-re)))
		     ,(make-qa0-operation attr* 'double-add
					  (list (make-reg a-re))
					  (list (make-reg b-re)
						(make-reg c-im)))
		     ,@r*)
		   env))))
   (define (c2r-complex-dot-init attr* output* input* r* env)
     (let-values* ([(a-re env) (c2r-rename env (car output*) 'real)]
		   [(a-im env) (c2r-rename env (car output*) 'imag)])
     (values `(,(make-qa0-operation attr* 'double-zero
				    (list (make-reg a-im)) '())
	       ,(make-qa0-operation attr* 'double-zero
				    (list (make-reg a-re)) '())
	       ,@r*)
	     env)))
   (define (c2r-complex-dot-add attr* output* input* r* env)
     (c2r-complex-cmadd attr* output* input* r* env))
   (define (c2r-complex-dot-fini attr* output* input* r* env)
     (let-values* ([(a-re env) (c2r-rename env (car output*) 'real)]
		   [(a-im env) (c2r-rename env (car output*) 'imag)]
		   [(b-re env) (c2r-rename env (car input*) 'real)]
		   [(b-im env) (c2r-rename env (car input*) 'imag)])
     (values `(,(make-qa0-operation attr* 'double-move
				    (list (make-reg a-im)) 
				    (list (make-reg b-im)))
	       ,(make-qa0-operation attr* 'double-move
				    (list (make-reg a-re))
				    (list (make-reg b-re)))
	       ,@r*)
	     env)))
   (define c2r-op*
     `((complex-move          . ,c2r-complex-move)
       (complex               . ,c2r-complex)
       (complex-real          . ,c2r-complex-real)
       (complex-imag          . ,c2r-complex-imag)
       (complex-neg           . ,c2r-complex-neg)
       (complex-times-plus-i  . ,c2r-complex-times-plus-i)
       (complex-times-minus-i . ,c2r-complex-times-minus-i)
       (complex-add           . ,c2r-complex-add)
       (complex-sub           . ,c2r-complex-sub)
       (complex-rmul          . ,c2r-complex-rmul)
       (complex-mul           . ,c2r-complex-mul)
       (complex-madd          . ,c2r-complex-madd)
       (complex-rmadd         . ,c2r-complex-rmadd)
       (complex-rmsub         . ,c2r-complex-rmsub)
       (complex-cmul          . ,c2r-complex-cmul)
       (complex-cmadd         . ,c2r-complex-cmadd)
       (complex-add-i         . ,c2r-complex-add-i)
       (complex-sub-i         . ,c2r-complex-sub-i)
       (complex-norm-init     . ,c2r-complex-norm-init)
       (complex-norm-add      . ,c2r-complex-norm-add)
       (complex-norm-fini     . ,c2r-complex-norm-fini)
       (complex-dot-init      . ,c2r-complex-dot-init)
       (complex-dot-add       . ,c2r-complex-dot-add)
       (complex-dot-fini      . ,c2r-complex-dot-fini)))
   (define (c2r-operation c attr* name output* input* r* env)
     (cond
      [(assq name c2r-op*)
       => (lambda (n&f) ((cdr n&f) attr* output* input* r* env))]
      [else (values (cons c r*) env)]))
   (define (complex->real ast env)
     (variant-case ast
       [qa0-top (decl*)
	 (let loop ([r* '()] [decl* decl*] [env env])
	   (cond
	    [(null? decl*) (values (make-qa0-top (reverse r*)) env)]
	    [else (let-values* ([(d e) (c2r-decl (car decl*) env)])
		    (loop (cons d r*) (cdr decl*) e))]))])))

