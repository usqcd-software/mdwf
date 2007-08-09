(load "common.ss")
(load "basis.ss")
(load "parser.ss")
(load "qa0print.ss")
(load "cfolding.ss")
(load "q2complex.ss")

(define (c-t name d/f)
  (let-values* ([(ast) (parse-qa0-file name)]
                [(ast env) (fold-constants/env ast ce-bgl d/f)]
                [(ast) (qcd->complex ast env)])
    (printf "~%**** Tree~%")
    (print-tree ast)
    (printf "~%**** Environment~%")
    (for-each (lambda (r) (printf "   ~s : ~s~%" (car r) (cdr r))) env)))

                