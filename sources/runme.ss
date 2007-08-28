(load "common.ss")
(load "basis.ss")
(load "parser.ss")
(load "qa0print.ss")
(load "cfolding.ss")
(load "q2complex.ss")
(load "backends.ss")
(load "mk-cheader.ss")
;(load "gen-header.ss")

;(define (c99-emit . args) (printf "~%/* TODO: c99-emit */~%"))

(define (c-t name machine d/f)
  (let-values* ([ast (parse-qa0-file name)]
                [(ast env) (fold-constants/env ast machine d/f)]
                [(ast e-x) (qcd->complex ast env)]
;		[_ (printf "~%**** Tree before back end~%")]
;		[_ (print-tree ast)]
;		[_ (printf "~%**** Environment before back end~%")]
;		[_ (for-each (lambda (r)
;			       (printf "   ~s : ~s~%" (car r) (cdr r))) e-x)]
;		[_ (printf "~%*****~%~%")]
                [(ast env) (complex->back-end ast env)]
;		[_ (printf "~%**** Tree at the back end~%")]
;		[_ (print-tree ast)]
;		[_ (printf "~%**** Environment after back end~%")]
;		[_ (for-each (lambda (r)
;			       (printf "   ~s : ~s~%" (car r) (cdr r))) env)]
;		[_ (printf "~%*****~%~%")]
		)
    (emit-back-end ast env)
    (printf "~%~%~%HEADER BEGIN~%~%")
    (mk-c-header ast env)
    (printf "~%~%~%HEADER END~%~%")
))

                