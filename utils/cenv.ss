(module cenv
	mzscheme
   (require "basis.ss")
   (require "common.ss")
   (require "ast.ss")
   (provide ce-empty-env
	    ce-search
	    ce-search-x
	    ce-lookup
	    ce-lookup-x
	    ce-bind
	    ce-bind-x
	    ce-add-const
	    ce-add-type
	    ce-add-array
	    ce-add-struct
	    ce-add-alias
	    ce-add-qcd-type
            ce-for-each
	    ce-bgl)
   (define (ce-empty-env) '())
   (define (ce-search env key k-found k-missed)
     (let ([x (assoc key env)])
       (if x (k-found (cdr x))
	   (k-missed))))
   (define (ce-search-x env type key k-found k-missed)
     (ce-search env `(,type ,key) k-found k-missed))
   (define (ce-lookup env key msg . args)
     (ce-search env key (lambda (x) x)
		(lambda () (apply error 'qa0 msg args))))
   (define (ce-lookup-x env type key msg . args)
     (ce-search env `(,type ,key) (lambda (x) x)
		(lambda () (apply error 'qa0 msg args))))
   (define (ce-bind env k v) `((,k . ,v) ,@env))
   (define (ce-bind-x env t k v) `(((,t ,k) . ,v) ,@env))
   (define (ce-add-const env name value)
     (let* ([t (list 'type name)]
	    [x (assoc t env)])
       (cond
	[x (error 'qa0 "Rebinding ~a to ~a is not allowed, old binding ~a"
		  name value (cdr x))]
	[else (let* ([env (ce-bind env t 'const)]
		     [env (ce-bind-x env 'const name value)])
		env)])))
   (define (ce-add-type env name c-name size align)
     (let* ([t (list 'type name)]
	    [x (assoc t env)])
       (cond
	[x (error 'qa0 "Redefining type ~a is not allowed" name)]
	[else (let* ([env (ce-bind env t 'type)]
		     [env (ce-bind-x env 'size-of name size)]
		     [env (ce-bind-x env 'align-of name align)]
		     [env (ce-bind-x env 'components name '())]
		     [env (ce-bind-x env 'aliased-to name name)]
		     [env (ce-bind-x env 'name-of name c-name)])
		env)])))
   (define (ce-add-array env name c-name base size)
     (let* ([bs (ce-lookup-x env 'size-of base "Size of array base ~a" base)]
	    [ba (ce-lookup-x env 'align-of base
			     "Alignment of array base ~a" base)]
	    [t (list 'type name)]
	    [x (assoc t env)])
       (cond
	[x (error 'qa0 "Redefining array ~a is not allowed" name)]
	[else (let* ([env (ce-bind env t 'array)]
		     [env (ce-bind-x env 'size-of name (* size bs))]
		     [env (ce-bind-x env 'align-of name ba)]
		     [env (ce-bind-x env 'components name '())]
		     [env (ce-bind-x env 'name-of name c-name)])
		env)])))
   (define (ce-add-struct env name c-name field* type*)
     (let* ([t (list 'type name)]
	    [x (assoc t env)])
       (cond
	[x (error 'qa0 "Redefining structure ~a is not allowed" name)]
	[else (let loop ([env env] [f* field*] [t* type*]
			 [size 0] [align 1])
		(cond
                 [(null? f*)
		  (let* ([env (ce-bind-x env 'type name 'struct)]
			 [env (ce-bind-x env 'size-of name size)]
			 [env (ce-bind-x env 'align-of name align)]
			 [env (ce-bind-x env 'components name field*)]
			 [env (ce-bind-x env 'name-of name c-name)])
		    env)]
                 [else
		  (let* ([f (car f*)] [t (car t*)]
			 [a-f (ce-lookup-x env 'align-of t
					   "Alignment of ~a.~a" name f)]
			 [s-f (ce-lookup-x env 'size-of t
					   "Size of field ~a.~a" name f)]
			 [start (* a-f (quotient (+ size a-f -1) a-f))]
			 [align (max a-f align)]
			 [env (ce-bind env (list 'offset-of name f) start)])
		    (loop env (cdr f*) (cdr t*) (+ start s-f) align))]))])))
   (define (ce-add-alias env new old)
     (let ([t (ce-lookup-x env 'type old "type of ~a" old)])
       (case t
	 [(const) (let* ([v (ce-lookup-x env 'const old "value of ~a" old)]
			 [env (ce-bind-x env 'type new t)]
			 [env (ce-bind-x env 'const new v)])
		    env)]
	 [(type array struct)
	  (let* ([c (ce-lookup-x env 'components old
				 "components of ~a" old)]
		 [s (ce-lookup-x env 'size-of old "size of ~a" old)]
		 [a (ce-lookup-x env 'align-of old "align of ~a" old)]
		 [x (ce-lookup-x env 'name-of old "C name of ~a" old)]
                 [tn (ce-lookup-x env 'aliased-to old "True name of ~a" old)]
		 [env (ce-bind-x env 'type new t)]
		 [env (ce-bind-x env 'size-of new s)]
		 [env (ce-bind-x env 'align-of new a)]
		 [env (ce-bind-x env 'components new c)]
		 [env (ce-bind-x env 'aliased-to new tn)]
		 [env (ce-bind-x env 'name-of new x)])
	    (let loop ([c c] [env env])
	      (cond
	       [(null? c) env]
	       [else (let ([o (ce-lookup env (list 'offset-of old (car c))
					 "offset of ~a in ~a"
					 (car c) old)])
		       (loop (cdr c)
			     (ce-bind env (list 'offset-of new (car c))
				      o)))])))]
	 [else (error 'qa0 "Internal error in ce-add-alias")])))
   (define (ce-add-qcd-type env name c-name a-dim b-dim)
     (let ([c-size (ce-lookup-x env 'size-of 'COMPLEX "(size-of COMPLEX)")]
	   [c-align (ce-lookup-x env 'align-of 'COMPLEX "(align-of COMPLEX)")]
	   [a-size (ce-lookup-x env 'const a-dim "(const ~a)" a-dim)]
	   [b-size (ce-lookup-x env 'const b-dim "(const ~a)" b-dim)])
       (ce-add-type env name c-name (* a-size b-size c-size) c-align)))
   (define (ce-for-each env predicate? proc)
     (let loop ([env env])
       (cond
	[(null? env) #t]
	[(predicate? (caar env) (cdar env)) (proc (caar env) (cdar env))
	 (loop (cdr env))]
	[else (loop (cdr env))])))
   (define (ce-bgl env)
     (let* ([env (ce-add-type env 'int            "int"              4  4)]
	    [env (ce-add-type env 'pointer        "void *"           4  4)]
	    [env (ce-add-type env 'float          "float"            4  4)]
	    [env (ce-add-type env 'double         "double"           8  8)]
	    [env (ce-add-type env 'vector-float   "vector float"     8  8)]
	    [env (ce-add-type env 'vector-double  "vector double"   16 16)]
	    [env (ce-add-type env 'complex-float  "float _Complex"   8  8)]
	    [env (ce-add-type env 'complex-double "double _Complex" 16 16)]
	    [env (ce-add-const env '*colors*            3)]
	    [env (ce-add-const env '*dim*               4)]
	    [env (ce-add-const env '*fermion-dim*       4)]
	    [env (ce-add-const env '*projected-fermion-dim*  2)])
       env)))