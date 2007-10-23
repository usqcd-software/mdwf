(module verbose
        mzscheme
   (provide process-verbose)

   (define (process-verbose key target* data*)
     (let loop ([target* target*] [data* data*])
       (cond
	[(null? target*) (printf "~%")]
	[else (if (eq? key (car target*)) (printf "~a" (car data*)))
	      (loop (cdr target*) (cdr data*))]))))
