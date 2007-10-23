(module cheader
        mzscheme
  (require "common.ss")
  (require "ast.ss")
  (require "cenv.ss")
  (require "attr.ss")
  (require "verbose.ss")

  (provide c-header
	   build-proc-name
	   build-proc-type)
  (define (c-header ast env)
    (define (gh-decl decl)
      (variant-case decl
        [qa0-array (name c-name base-name size)
          (gh-add-array name c-name base-name size)]
        [qa0-struct (name c-name field-type* field-c-name*)
          (gh-add-struct  name c-name field-type* field-c-name*)]
        [qa0-proc (attr* name arg-c-name* arg-c-type*)
          (gh-add-proc attr* name arg-c-name* arg-c-type*)]
	[qa0-verbose (target* data*) (process-verbose 'c-header target* data*)]
	[else #f]))
    (define (gh-add-array name c-name base size)
      (printf "typedef ~a ~a[~a];~%~%"
	      (ce-lookup-x env 'name-of base
			   "make-c-header: Internal error: no name for ~a" base)
	      c-name (variant-case size
                       [c-expr-number (number) number])))
    (define (gh-add-struct name c-name field-type* field-c-name*)
      (printf "struct ~a {~%" c-name)
      (for-each (lambda (n t) 
		  (printf "  ~a ~a;~%" (ce-lookup-x env 'name-of n
						    "Missing name for type ~a"
						    n) t))
		field-type* field-c-name*)
      (printf "};~%~%"))
    (define (gh-add-proc attr* name arg-c-name* arg-c-type*)
      (printf "~a ~a"
	      (build-proc-type name attr* env)
	      (build-proc-name name attr* env))
      (if (null? arg-c-name*) (printf "(void);~%")
	  (begin (printf "(")
		 (let loop ([name* arg-c-name*] [type* arg-c-type*])
		   (cond
		    [(null? name*)]
		    [else (printf "~a~a ~a" (if (eq? name* arg-c-name*) "" ", ")
				  (car type*) (car name*))
			  (loop (cdr name*) (cdr type*))]))
		 (printf ");~%"))))
    (variant-case ast
      [qa0-top (decl*) (let loop ([in* decl*])
			 (cond
			  [(null? in*) #t]
			  [else (gh-decl (car in*)) (loop (cdr in*))]))]))
  (define (build-proc-type name attr* env)
    (let ([f-count (attr-search attr* 'count-flops (lambda (v*) #t)
				(lambda () #f))]
	  [r-type (attr-search attr* 'return
			       (lambda (v*)
				 (if (< (length v*) 3)
				     (error 'proc-type: "short attribute ~a"
					    v*)
				     (caddr v*)))
			       (lambda () #f))])
      (cond
       [(and f-count r-type)
	(error 'proc-type "Both count-flops and return specified for ~a" name)]
       [f-count "unsigned int"]
       [r-type r-type]
       [else "void"])))
  (define (build-proc-name name attr* env)
    (let ([prefix (ce-lookup env 'proc-prefix "missing proc-prefix in CE")]
	  [colors (ce-lookup-x env 'const '*colors* "missing *colors* in CE")]
	  [prec   (ce-lookup env 'prec-letter "missing prec-letter in CE")]
	  [infix  (ce-lookup env 'proc-infix "missing proc-infix in CE")]
	  [suffix (ce-lookup env 'proc-suffix "missing proc-suffix in CE")]
	  [stem*  (attr-lookup attr* 'stem "missing stem for ~a" name)]
	  [has-df? (attr-search attr* 'prec&color
				(lambda (v*) #t) (lambda () #f))])
      (let loop ([r (if has-df?
			(format "~a_~a~a_~a_" prefix prec colors infix)
			(format "~a_~a_" prefix infix))]
		 [stem* stem*])
	(cond
	 [(null? stem*) (format "~a~a" r suffix)]
	 [else (loop (format "~a~a" r (car stem*)) (cdr stem*))])))))
