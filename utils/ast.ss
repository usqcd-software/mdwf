(module ast
        mzscheme
   (require "common.ss")
   (define-variant qa0-top (decl*))
   (define-variant qa0-alias (old-name new-name))
   (define-variant qa0-const (name value))
   (define-variant qa0-struct (name c-name
			       field-name* field-type* field-c-name*))
   (define-variant qa0-array (name c-name base-name size))
   (define-variant qa0-proc (attr* name
			     arg-name* arg-type* arg-c-name* arg-c-type*
			     code*))
   (define-variant qa0-repeat-list (id value* body*))
   (define-variant qa0-repeat-range (id low high body*))
   (define-variant qa0-operation (attr* name output* input*))
   (define-variant qa0-load (attr* type output addr*))
   (define-variant qa0-store (attr* type addr* value))
   (define-variant qa0-loop (attr* var low high code*))
   (define-variant qa0-if (var true-code* false-code*))
   (define-variant qa0-macro-list (id value* code*))
   (define-variant qa0-macro-range (id low high code*))
   (define-variant qa0-attr (name value*))
   (define-variant c-expr-id (id))
   (define-variant c-expr-quote (literal))
   (define-variant c-expr-number (number))
   (define-variant c-expr-string (string))
   (define-variant c-expr-op (name c-expr*))
   (define-variant reg (name)))
