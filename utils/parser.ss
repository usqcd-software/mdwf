(module parser
        mzscheme
   (require "common.ss")
   (require "ast.ss")
   (require "cenv.ss")
   (provide parse-qa0-file
	    user-reg)

   (define (user-reg name) (string->symbol (format "_~a" name)))
   (define (parse-qa0-file file-name)
     (define (check-list msg in cmp min-size)
       (if (not (list? in))
	   (error 'parse-qa0 "List expected in ~a, found~%  ~a~%"
		  msg in))
       (if (not (cmp (length in) min-size))
	   (error 'parse-qa0
		  "List of ~a ~a is expected in ~a, found:~%   ~a~%"
		  (if (eq? cmp =) "exactly" "at least")
		  min-size msg in)))
     (define (check-name msg form name)
       (if (not (symbol? name))
	   (error 'parse-qa0
		  "Expecting name for ~a, found ~a in:~%~a~%"
		  msg name form)))
     (define (check-sym-or-list-of-syms msg form name)
       (cond
	[(symbol? name) #t]
	[(and (list? name) (let loop ([name name])
			     (cond
			      [(null? name) #t]
			      [(symbol? (car name)) (loop (cdr name))]
			      [else #f]))) #t]
	[else (error 'parse-qa0
		     "Expecting name or list of names for ~a, found ~a in ~a~%"
		     msg name form)]))
     (define (check-string msg form string)
       (if (not (string? string))
	   (error 'parse-qa0
		  "Expecting string for ~a, found ~a in:~%~a~%"
		  msg string form)))
     (define (check-attrib* msg form attr*)
       (for-each (lambda (attr) (check-attrib msg form attr)) attr*))
     (define (check-attrib msg form attr)
       (if (not (or (symbol? attr)
		    (and (list? attr)
			 (>= (length attr) 1)
			 (symbol? (car attr)))))
	   (error 'parse-qa0 "Expecting attribute for ~a, found ~a in~%~a~%"
		  msg attr form)))
     (define (check-output* msg form out*)
       (cond
	[(null? out*) #t]
	[(not (symbol? (car out*)))
         (error 'parse-qa0 "Bad value for output of ~a, found ~a in~%~a~%"
		msg (car out*) form)]
	[else (check-output* msg form (cdr out*))]))
     (define (check-input* msg form in*)
       (for-each (lambda (in) (check-input msg form in)) in*))
     (define (check-input msg form in)
       (check-list msg in = 2)
       (if (case (car in)
	     [(reg) (not (symbol? (cadr in)))]
	     [(const) #f]
	     [else #t])
	   (error 'parse-qa0 "Bad value for input of ~a, found ~a in ~%~a~%"
		  msg in form)))
     (define (check-iterator msg form iter)
       (check-list msg iter >= 2)
       (check-name msg form (car iter))
       (case (length iter)
	 [(2) (check-list "iterator range" (cadr iter) >= 0)]
	 [(3) #t]
	 [else (error 'parse-qa0 "Bad iterator form of ~a, found ~a in~%~a~%"
		      msg iter form)]))
     (define (empty-ast) (make-qa0-top '()))
     (define (parse-top s-expr ast)
       (check-list "top level" s-expr >= 1)
       (variant-case ast
         [qa0-top (decl*)
	   (case (car s-expr)
	     [(alias)
	      (make-qa0-top (append decl* (list (parse-alias s-expr))))]
	     [(constant)
	      (make-qa0-top (append decl* (list (parse-constant s-expr))))]
	     [(array)
	      (make-qa0-top (append decl* (list (parse-array s-expr))))]
	     [(structure)
	      (make-qa0-top (append decl* (list (parse-structure s-expr))))]
	     [(repeat procedure)
	      (make-qa0-top (append decl* (parse-repeat s-expr)))]
	     [(include)
              (make-qa0-top (append decl* (parse-include s-expr)))]
	     [(verbose)
	      (make-qa0-top (append decl* (list (parse-verbose s-expr))))]
	     [else (error 'parse-qa0 "Unknown form:~%~a~%" s-expr)])]))
     (define (parse-include s-expr)
       (check-list "include form" s-expr = 2)
       (check-string "include file name" s-expr (cadr s-expr))
       (let ([f (open-input-file (cadr s-expr))])
	 (let loop ([r (read f)] [ast (empty-ast)])
	   (cond
	    [(eof-object? r) (close-input-port f) (qa0-top->decl* ast)]
	    [else (loop (read f) (parse-top r ast))]))))
     (define (parse-verbose s-expr)
       (check-list "verbose form" s-expr >= 1)
       (map (lambda (x)
	      (check-list "verbose case" x = 2)
	      (check-sym-or-list-of-syms "verbose target" s-expr (car x))
	      (check-string "verbose data" s-expr (cadr x)))
	    (cdr s-expr))
       (make-qa0-verbose (map car (cdr s-expr)) (map cadr (cdr s-expr))))
     (define (parse-alias s-expr)
       (check-list "top level alias" s-expr = 3)
       (check-name "alias new name" s-expr (cadr s-expr))
       (check-name "alias old name" s-expr (caddr s-expr))
       (make-qa0-alias (caddr s-expr) (cadr s-expr)))
     (define (parse-constant s-expr)
       (check-list "constant" s-expr = 3)
       (check-name "constant name" s-expr (cadr s-expr))
       (make-qa0-const (cadr s-expr) (parse-const-expr (caddr s-expr))))
     (define (parse-array s-expr)
       (check-list "array" s-expr = 5)
       (check-name "array name" s-expr (cadr s-expr))
       (check-string "array c-name" s-expr (caddr s-expr))
       (check-name "array base name" s-expr (cadddr s-expr))
       (make-qa0-array (cadr s-expr) (caddr s-expr) (cadddr s-expr)
		       (parse-const-expr (car (cddddr s-expr)))))
     (define (parse-structure s-expr)
       (check-list "structure" s-expr = 4)
       (let ([name (cadr s-expr)]
	     [c-name (caddr s-expr)]
	     [f* (cadddr s-expr)])
	 (check-name "struct name" s-expr name)
	 (check-string "struct c-name" s-expr c-name)
	 (map (lambda (f)
		(check-list "struct field" f = 3)
		(check-name "struct field name" s-expr (car f))
		(check-string "struct field c-name" s-expr (cadr f))
		(check-name "struct field type" s-expr (caddr f)))
	      f*)
	 (make-qa0-struct name c-name
			  (map car f*) (map caddr f*) (map cadr f*))))
     (define (parse-repeat s-expr)
       (check-list "top level repeat" s-expr >= 1)
       (case (car s-expr)
	 [(repeat)
	  (check-list "top level repeat" s-expr >= 3)
	  (let loop ([bind* (cadr s-expr)])
	    (cond
             [(null? bind*)
	      (let do-body ([p* '()] [b* (cddr s-expr)])
		(cond
		 [(null? b*) p*]
		 [else (do-body (append p* (parse-repeat (car b*)))
				(cdr b*))]))]
             [else
	      (check-iterator "top level repeat" s-expr (car bind*))
	      (list
	       (case (length (car bind*))
		 [(2) (make-qa0-repeat-list (caar bind*)
					    (map parse-a-value (cadar bind*))
					    (loop (cdr bind*)))]
		 [(3) (make-qa0-repeat-range (caar bind*)
					     (parse-const-expr (cadar bind*))
					     (parse-const-expr (caddar bind*))
					     (loop (cdr bind*)))]
		 [else (error 'qa0 "Internal error in parse-repeat")]))]))]
	 [(procedure) (list (parse-procedure s-expr))]
	 [else (error 'qa0 "Bad form at top level~%~a~%" s-expr)]))
     (define (parse-a-value expr)
       (if (or (number? expr)
	       (symbol? expr)
	       (string? expr)) expr
	       (error 'parse-qa0 "a-value is not ~a~%" expr)))
     (define (parse-procedure s-expr)
       (check-list "procedure" s-expr >= 5)
       (check-name "procedure name" s-expr (cadr s-expr))
       (check-attrib* "procedure attributes" s-expr (caddr s-expr))
       (check-list "procedure arguments" (cadddr s-expr) >= 0)
       (map (lambda (arg)
	      (check-list "procedure argument" arg = 4)
	      (check-name "argument name" arg (car arg))
	      (check-name "argument type" arg (cadr arg))
	      (check-string "argument c-name" arg (caddr arg))
	      (check-string "argument c-type" arg (cadddr arg)))
	    (cadddr s-expr))
       (parse-code* (cddddr s-expr)
		    (lambda (code*)
		      (make-qa0-proc (map parse-attr (caddr s-expr))
				     (cadr s-expr)
				     (map user-reg (map car (cadddr s-expr)))
				     (map cadr (cadddr s-expr))
				     (map caddr (cadddr s-expr))
				     (map cadddr (cadddr s-expr))
				     code*))))
     (define (parse-code* code* k)
       (if (null? code*) (k '())
	   (parse-code* (cdr code*)
			(lambda (op*)
			  (let ([code (car code*)])
			    (check-list "operation" code >= 1)
			    (case (car code)
			      [(begin)
			       (parse-code* (cdr code)
					    (lambda (op2*)
					      (k (append op2* op*))))]
			      [(macro) (k (parse-macro code op*))]
			      [else (k (cons (parse-code code) op*))]))))))
     (define (parse-attr f)
       (cond
	[(symbol? f) (make-qa0-attr f '())]
	[else (check-list "attribute" f >= 1)
	      (make-qa0-attr (car f) (map parse-a-value (cdr f)))]))
     (define (parse-code f)
       (check-list "code" f >= 1)
       (case (car f)
	 [(op) (parse-op f)]
	 [(load) (parse-load f)]
	 [(store) (parse-store f)]
	 [(loop) (parse-loop f)]
	 [(if) (parse-if f)]
	 [(if-else) (parse-if-else f)]
	 [else (error 'parse-qa0 "Unexpected code~%~a~%" f)]))
     (define (parse-op f)
       (check-list "code op" f = 5)
       (let ([name (cadr f)]
	     [attr* (caddr f)]
	     [out* (cadddr f)]
	     [in* (car (cddddr f))])
	 (check-name "code op name" f name)
	 (check-attrib* "code op attirbutes" f attr*)
	 (check-output* "code op outputs" f out*)
	 (check-input* "code op inputs" f in*)
	 (make-qa0-operation (map parse-attr attr*)
			     name
			     (map parse-output out*)
			     (map parse-input in*))))
     (define (parse-output out) (parse-reg out))
     (define (parse-input in)
       (case (car in)
	 [(reg) (parse-reg (cadr in))]
	 [(const) (parse-const-expr in)]
	 [else (error 'parse-qa0 "Internal error in parse-input")]))
     (define (parse-reg name) (make-reg (user-reg name)))
     (define (parse-load f)
       (check-list "load" f = 5)
       (let ([type (cadr f)]
	     [attr* (caddr f)]
	     [out (cadddr f)]
	     [addr (car (cddddr f))])
	 (check-attrib* "load" f attr*)
	 (check-name "load" f type)
	 (check-output* "load" f (list out))
	 (check-input* "load" f addr)
	 (make-qa0-load (map parse-attr attr*)
			type
			(parse-output out)
			(map parse-input addr))))
     (define (parse-store f)
       (check-list "store" f = 5)
       (let ([type (cadr f)]
	     [attr* (caddr f)]
	     [addr (cadddr f)]
	     [value (car (cddddr f))])
	 (check-attrib* "store" f attr*)
	 (check-name "store" f type)
	 (check-input* "store" f addr)
	 (check-input "store" f value)
	 (make-qa0-store (map parse-attr attr*)
			 type
			 (map parse-input addr)
			 (parse-input value))))
     (define (parse-loop f)
       (check-list "loop" f >= 4)
       (let ([attr* (cadr f)]
	     [ctl (caddr f)]
	     [code* (cdddr f)])
	 (check-attrib* "loop" f attr*)
	 (check-list "loop control" ctl = 3)
	 (parse-code* code* (lambda (op*)
			      (make-qa0-loop (map parse-attr attr*)
					     (parse-output (car ctl))
					     (parse-input (cadr ctl))
					     (parse-input (caddr ctl))
					     op*)))))
     (define (parse-if f)
       (check-list "if" f = 3)
       (let ([c (cadr f)]
	     [t (caddr f)])
	 (check-input "if predicate" f c)
	 (check-list "if branch" t >= 1)
	 (make-qa0-if (parse-input c) (parse-seq t) '())))
     (define (parse-if-else form)
       (check-list "if-else" form = 4)
       (let ([c (cadr form)]
	     [t (caddr form)]
	     [f (cadddr form)])
	 (check-input "if-else predicate" form c)
	 (check-list "if-else true branch" t >= 1)
	 (check-list "if-else false branch" f >= 1)
	 (make-qa0-if (parse-input c) (parse-seq t) (parse-seq f))))
     (define (parse-seq s)
       (if (and (list s) (> (length s) 0) (equal? (car s) 'begin))
	   (parse-code* (cdr s) (lambda (p*) p*))
	   (list (parse-code s))))
     (define (parse-macro form op*)
       (check-list "macro" form >= 3)
       (let ([iterator* (cadr form)]
	     [code* (cddr form)])
	 (let loop ([i* iterator*] [t* op*])
	   (cond
	    [(null? i*) (parse-code* code* (lambda (c*) (append c* t*)))]
	    [else (let ([i (car i*)] [i* (cdr i*)])
		    (check-iterator "macro" form i)
		    (case (length i)
		      [(2) (let ([c* (loop i* '())]
				 [v* (map parse-a-value (cadr i))]
				 [id (car i)])
			     (cons (make-qa0-macro-list id v* c*) t*))]
		      [(3) (let ([c* (loop i* '())]
				 [id (car i)]
				 [lo (parse-const-expr (cadr i))]
				 [hi (parse-const-expr (caddr i))])
			     (cons (make-qa0-macro-range id lo hi c*) t*))]
		      [else (error 'qa0 "Internal error in parse-macro")]))]))))
     (define (parse-const-expr form)
       (check-list "constant expression" form = 2)
       (parse-c-expr (cadr form)))
     (define (parse-c-expr form)
       (cond
	[(number? form) (make-c-expr-number form)]
	[(symbol? form) (make-c-expr-id form)]
	[(string? form) (make-c-expr-string form)]
	[(and (= (length form) 2) (eq? (car form) 'quote))
	 (make-c-expr-quote (cadr form))]
	[(and (list? form) (>= (length form) 1) (symbol? (car form)))
	 (make-c-expr-op (car form) (map parse-c-expr (cdr form)))]
	[else (error 'parse-qa0 "constant expression is bad:~%~a~%" form)]))
     (let ([f (open-input-file file-name)])
       (let loop ([s-expr (read f)] [ast (empty-ast)])
	 (cond
	  [(eof-object? s-expr) ast]
	  [else (loop (read f) (parse-top s-expr ast))])))))
