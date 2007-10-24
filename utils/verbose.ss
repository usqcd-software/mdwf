(module verbose
        mzscheme
   (provide process-verbose)

   (define (process-verbose key target* data*)
     (let loop ([target* target*] [data* data*])
       (cond
	[(null? target*) (printf "~%")]
	[(or (eq? key (car target*))
	     (and (list? (car target*))
		  (memq key (car target*))))
	 (printf "~a" (car data*))
	 (loop (cdr target*) (cdr data*))]
	[else (loop (cdr target*) (cdr data*))]))))
