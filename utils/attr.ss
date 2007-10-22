(module attr
        mzscheme
   (require "common.ss")
   (require "ast.ss")
   (require "cenv.ss")
   (provide attr-search
            attr-lookup)

   (define (attr-search attr* key found missed)
     (let loop ([attr* attr*])
       (cond
	[(null? attr*) (missed)]
	[else (variant-case (car attr*)
               [qa0-attr (name value*) (if (eq? key name) (found value*)
					   (loop (cdr attr*)))])])))
   (define (attr-lookup attr* key msg . arg*)
     (attr-search attr* key
		  (lambda (v*) v*)
		  (lambda () (apply error 'qa0 msg arg*)))))
