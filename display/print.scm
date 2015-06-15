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

;;; Hamiltonians look better if we divide them out.

(define *divide-out-terms* #t)

(define *fully-divide-out-terms* #t)

(define (ham:simplify hexp)
  (cond ((and (quotient? hexp) *divide-out-terms*)
	 (cond ((sum? (symb:numerator hexp))
		(let ((d (symb:denominator hexp)))
		  (a-reduce symb:+
			    (map (lambda (n)
				   (simplify (symb:/ n d)))
				 (operands (symb:numerator hexp))))))
	       (else hexp)))
	((and (compound-data-constructor? hexp) *fully-divide-out-terms*)
	 (cons (operator hexp) (map ham:simplify (operands hexp))))
	(else hexp)))


;;; Equations are often prettier if we get rid of the denominators,
;;; but watch out for singularities.

(define (eqn:simplify hexp)
  (cond ((quotient? hexp)
	 (symb:numerator hexp))
	((matrix? hexp)
	 ((m:elementwise eqn:simplify) hexp))
	((vector? hexp)
	 ((v:elementwise eqn:simplify) hexp))
	(else hexp)))

(define (flush-derivative expr)
  (substitute derivative-symbol
	      'derivative
	      expr))


;;; Preserves row/col distinction.
(define *flushing-row/col* #f)

(define (flush-row/col expr)
  (if *flushing-row/col*
      (let lp ((expr expr))
	(if (pair? expr)
	    (if (memq (car expr) '(row column))
		(cons 'vector
		      (map lp (cdr expr)))
		(map lp expr))
	    expr))
      expr))

(define (simplify exp)
  (flush-derivative
   (flush-row/col
    (ham:simplify
     (easy-simplify			;from rule-environment
      (expression exp))))))


(define *only-printing* #f)
(define *last-expression-printed*)

(define (canonicalize-numbers expr)
  (cond ((pair? expr)
	 (cons (canonicalize-numbers (operator expr))
	       (map canonicalize-numbers (operands expr))))
	((number? expr)
	 (heuristic-canonicalize-complex expr))
	(else
	 expr)))

(define (prepare-for-printing expr simplifier)
  (set! *last-expression-printed* 
	(if (memq expr '(#t #f))
	    expr
	    (let ((nexpr
		   (canonicalize-numbers (expression expr))))
	      (simplifier nexpr))))
  *last-expression-printed*)



(define (show-expression expr #!optional simplifier)
  (if (default-object? simplifier)
      (set! simplifier simplify))
  (prepare-for-printing expr simplifier)
  (cond (*only-printing*
	 (pp *last-expression-printed*))
	(else
	 (internal-show-expression *last-expression-printed*))))

(define (print-expression expr #!optional simplifier)
  (if (default-object? simplifier)
      (set! simplifier simplify))
  (prepare-for-printing expr simplifier)
  (pp *last-expression-printed*))

(define pe print-expression)
(define se show-expression)
