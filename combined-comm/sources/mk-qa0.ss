(require (lib "pretty.ss"))
(define *output* "qa0.ss")
(define *input* "qa0.xx")

(define (process-file name)
  (let ([f (open-input-file name)])
    (let loop ([r (read f)])
      (cond
        [(eof-object? r) (close-input-port f)]
        [else (pretty-print r) (loop (read f))]))))

(if (file-exists? *output*)
    (delete-file *output*))

(with-output-to-file *output*
  (lambda ()
    (printf "(~%")
    (let ([f (open-input-file *input*)])
      (let loop ([r (read f)])
	(cond
	 [(null? r) (close-input-port f)]
	 [(and (list? (car r)) (= (length (car r)) 2) (eq? (caar r) 'include))
	  (process-file (cadar r)) (loop (cdr r))]
	 [else (pretty-print (car r)) (loop (cdr r))])))
    (printf ")~%")))

      