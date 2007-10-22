(module common
	mzscheme
  (provide dbg
	   assert
	   define-variant
	   variant-type
	   variant-case
	   let-values*)
	   
  (define *verbose?* #f)
  (define (progress-report fmt . args)
    (if *verbose?*
	(begin (printf "  ")
	       (apply printf format args)
	       (printf "...~%"))))
  (define (dbg . msg)
    (display "D:")
    (for-each (lambda (x) (display " ") (write x)) msg)
    (newline)
    #f)
  (define-syntax assert
    (lambda (x)
      (syntax-case x ()
	[(_ test fmt arg ...)
	 (syntax (if (not test) (error 'assert (string-append "FAILED: " fmt)
				       arg ...)))])))
  (define-syntax define-variant
    (lambda (x)
      (define (mk-name template-id . args)
	(datum->syntax-object template-id
			      (string->symbol (apply string-append
                                (map (lambda (x)
				       (cond
                                        [(string? x) x]
                                        [(symbol? x) (symbol->string x)]
                                        [else (symbol->string
					       (syntax-object->datum x))]))
				     args)))))
      (syntax-case x ()
	[(_ name (field ...))
	 (with-syntax
           ([constructor (mk-name (syntax name) 'make- (syntax name))]
	    [predicate (mk-name (syntax name) (syntax name) "?")]
	    [(reader ...) (map (lambda (fld)
				 (mk-name (syntax name) (syntax name) "->" fld))
			       (syntax->list (syntax (field ...))))]
	    [count (length (syntax->list (syntax (name field ...))))])
	   (with-syntax
             ([(index ...) (let f ([i 1])
			     (if (= i (syntax-object->datum (syntax count))) '()
				 (cons i (f (+ 1 i)))))])
	     (syntax
	      (begin
                (provide constructor
                         predicate
                         reader ...)
		(define constructor (lambda (field ...)
				      (vector 'name field ...)))
		(define predicate (lambda (object)
				    (and (vector? object)
					 (= (vector-length object) count)
					 (eq? (vector-ref object 0) 'name))))
		(define reader (lambda (object) (vector-ref object index)))
		...))))])))
  (define variant-type (lambda (variant) (vector-ref variant 0)))
  (define-syntax variant-case
    (lambda (x)
      (define (mk-name template-id . args)
	(datum->syntax-object template-id
          (string->symbol (apply string-append
				 (map (lambda (x)
					(cond
					 [(string? x) x]
					 [(symbol? x) (symbol->string x)]
					 [else (symbol->string
						(syntax-object->datum x))]))
				      args)))))
      (syntax-case x (else)
	[(_ exp clause ...) (not (identifier? (syntax exp)))
	 (syntax (let ([var exp]) (_ var clause ...)))]
	[(_ var) (syntax (error 'variant-case "no clause matches ~a" var))]
	[(_ var (else exp exp1 ...)) (syntax (begin exp exp1 ...))]
	[(_ var (name [(fname field) ...] exp1 exp2 ...) clause ...)
	 (with-syntax
	  ([predicate (mk-name (syntax name) (syntax name) "?")]
	   [(reader ...) (map (lambda (fld)
				(mk-name (syntax name) (syntax name) "->" fld))
			      (syntax->list (syntax (field ...))))])
	  (syntax (if (predicate var)
		      (let ([fname (reader var)] ...) exp1 exp2 ...)
		      (_ var clause ...))))]
	[(_ var (name (field ...) exp1 exp2 ...) clause ...)
	 (syntax (_ var (name ([field field] ...) exp1 exp2 ...)
		    clause ...))])))
  (define-syntax let-values*
    (lambda (x)
      (syntax-case x ()
	[(_ () e ...) (syntax (let () e ...))]
	[(_ ([(n ...) i ...]) e ...)
	 (syntax (call-with-values (lambda () (let () i ...))
		   (lambda (n ...) (let () e ...))))]
	[(_ ([n i ...]) e ...) (syntax (let ([n (let () i ...)])
					 (let () e ...)))]
	[(_ (b0 b1 ...) e ...) (syntax (_ (b0) (_ (b1 ...) e ...)))]))))
