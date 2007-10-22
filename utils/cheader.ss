(module cheader
        mzscheme
  (require "common.ss")
  (require "ast.ss")
  (require "cenv.ss")
  (require "attr.ss")
  (require "parser.ss")

  (provide c-header)
  (define (c-header ast env)
    (define (gh-decl decl)
      (variant-case decl
        [qa0-array (name c-name base-name size)
          (gh-add-array name c-name base-name size)]
        [qa0-struct (name c-name field-type* field-c-name*)
          (gh-add-struct  name c-name field-type* field-c-name*)]
        [qa0-proc (attr* name arg-c-name* arg-c-type*)
          (gh-add-proc attr* name arg-c-name* arg-c-type*)]
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
      (printf "~a ~a" (attr-search attr* 'count-flops
                                   (lambda (v) "int")
                                   (lambda () "void"))
	      (build-proc-name (attr-lookup attr* 'stem
					    "missing stem in proc ~a" name)
			       env))
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
			  [else (gh-decl (car in*)) (loop (cdr in*))]))])))
