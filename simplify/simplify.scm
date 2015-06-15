#| -*-Scheme-*-

$Id$

Copyright (c) 2002 Massachusetts Institute of Technology

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
|#

;;;;       General Recursive Simplifier Maker

;;; Given a set of operations, this procedure makes a recursive
;;;  simplifier that simplifies expressions involving these
;;;  operations, treating other combinations as atomic.

;;; To break an expression up into manipulable and nonmanipulable
;;; parts with respect to a set of algebraic operators.  This is done
;;; by the introduction of auxiliary variables.

;;; For example, the equation
;;;    I = Is (exp((V2 - V3)/Vt) - 1) ; I, V2, V3
;;; can be broken into three equations
;;;    I + Is = Is*X                  ; I, X
;;;    V2/Vt - V3/Vt = Y              ; V2, V3, Y
;;;    X = (exp Y)                    ; X, Y

;;; where X and Y are new variables.  The first two parts contain only
;;; addition, subtraction, multiplication, and division and the third
;;; is not expressible in terms of those operations.


(declare (usual-integrations))

(define (make-analyzer ->expression expression-> known-operators)
  (let ((auxiliary-variable-table '())
	(uorder '())
	(priority '()))

    ;; Default simplifier
    (define (simplify expr)
      (new-analysis)
      (simplify-expression expr))


    ;; Simplify relative to existing tables
    (define (simplify-expression expr)	
      (backsubstitute (analyze-expression expr)))


    ;; Analyze relative to existing tables
    (define (analyze-expression expr)
      (fluid-let ((incremental-simplifier #f))
	(base-simplify (analyze expr))))


    ;; Set up new analysis
    (define (new-analysis)		
      (set! auxiliary-variable-table '())
      (set! uorder '())
      (set! priority '())
      'done)


    ;; Define ordering of variables
    (define (set-priority! . exprs)
      (set! priority (map add-symbol! exprs))
      priority)


    ;; Get kernel table
    (define (get-auxiliary-variable-defs)
      (map (lambda (entry)
	     (list (cdr entry) (car entry)))
	   auxiliary-variable-table))

    ;; Implementation -----------------------

    (define (analyze expr)
      (let ((vars (sort (variables-in expr) alphaless?)))
	(set! uorder
	      (append (map add-symbol! priority)
		      vars)))
      (ianalyze expr))

    (define (ianalyze expr)
      (if (and (pair? expr) (not (eq? (car expr) 'quote)))
	  (let ((sexpr (map ianalyze expr)))
	    ;; At this point all subexpressions are canonical.
	    (if (and (memq (operator sexpr) known-operators)
		     (not (and (expt? sexpr)
			       (not (integer?
				     (second (operands sexpr)))))))
		sexpr
		(let ((as-seen (expression-seen sexpr)))
		  (if as-seen
		      as-seen
		      (new-kernels sexpr)))))
	  expr))

    (define (new-kernels expr)
      (let ((sexpr (map base-simplify expr)))
	(let ((v (hash-table/get symbolic-operator-table
				 (operator sexpr)
				 #f)))
	  (if v
	      (let ((w (apply v (operands sexpr))))
		(if (and (pair? w) (eq? (operator w) (operator sexpr)))
		    (add-symbols! w)
		    (ianalyze w)))		      
	      (add-symbols! sexpr)))))

    (define (base-simplify expr)
      (if (and (pair? expr) (not (eq? (car expr) 'quote)))
	  (expression-> expr ->expression vless?)
	  expr))

    (define (backsubstitute expr)
      (cond ((pair? expr) (map backsubstitute expr))
	      ((symbol? expr)
	       (let ((v (rassq expr auxiliary-variable-table)))
		 (if v
		     (backsubstitute (car v))
		     ;Can improve by doing only once.
		     ; Could make new table on entry to back.
		     expr)))
	      (else expr)))

    (define (add-symbols! expr)
      (let ((new (map add-symbol! expr)))
	(add-symbol! new)))

    (define (add-symbol! expr)
      (if (and (pair? expr) (not (eq? (car expr) 'quote)))
	  (let ((as-seen (expression-seen expr)))
	    (if as-seen
		as-seen
		(let ((newvar
		       (generate-uninterned-symbol "kernel")))
		  (set! auxiliary-variable-table
			(cons (cons expr newvar) auxiliary-variable-table))
		  newvar)))
	  expr))

    (define (expression-seen expr)
      (let ((entry (assoc expr auxiliary-variable-table)))
	(and entry (cdr entry))))

    (define (vless? var1 var2)
      (let ((in (memq var1 uorder)))
	(cond (in
	       (cond ((memq var2 in) true)
		     ((memq var2 uorder) false)
		     (else true)))
	      ((memq var2 uorder) false)
	      (else
	       (let ((in1 (rassq var1 auxiliary-variable-table))
		     (in2 (rassq var2 auxiliary-variable-table)))
		 (cond ((and in1 in2)
			(memq in2 (memq in1 auxiliary-variable-table)))
		       (in1 true)
		       (in2 false)
		       (else
			(alphaless? var1 var2))))))))

    (vector simplify
	    simplify-expression
	    new-analysis
	    set-priority!
	    analyze-expression
	    get-auxiliary-variable-defs)))


(define (default-simplifier analyzer) (vector-ref analyzer 0))

(define (expression-simplifier analyzer) (vector-ref analyzer 1))

(define (initializer analyzer) (vector-ref analyzer 2))

(define (priority-setter analyzer) (vector-ref analyzer 3))

(define (expression-analyzer analyzer) (vector-ref analyzer 4))

(define (auxiliary-variable-fetcher analyzer) (vector-ref analyzer 5))

(define fpf:analyzer
  (make-analyzer fpf:->expression fpf:expression-> fpf:operators-known))

(define fpf:simplify (default-simplifier fpf:analyzer))


(define pcf:analyzer
  (make-analyzer poly:->expression poly:expression-> poly:operators-known))

(define pcf:simplify (default-simplifier pcf:analyzer))


(define rcf:analyzer
  (make-analyzer rcf:->expression rcf:expression-> rcf:operators-known))

(define rcf:simplify (default-simplifier rcf:analyzer))
#|
((initializer rcf:analyzer))

(pp ((expression-analyzer rcf:analyzer)
     '(- i (* Is (- (exp (/ (- v2 v3) Vt)) 1)))))
(+ (* (+ 1 (* -1 kernel17)) Is) i)

(pp ((auxiliary-variable-fetcher rcf:analyzer)))
((kernel17 (exp kernel16))
 (kernel16 (/ (+ v2 (* -1 v3)) Vt)))

(pp ((expression-analyzer rcf:analyzer)
     '(exp (/ (- v3 v2) (- Vt)))))
kernel17

(pp ((expression-simplifier rcf:analyzer)
     '(- i (* Is (- (exp (/ (- v2 v3) Vt)) 1)))))
(+ (* (+ 1 (* -1 (exp (/ (+ v2 (* -1 v3)) Vt)))) Is) i)
;Unspecified return value
|#