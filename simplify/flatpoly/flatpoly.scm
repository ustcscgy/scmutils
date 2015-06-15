#| -*-Scheme-*-

$Id: copyright.scm,v 1.5 2005/09/25 01:28:17 cph Exp $

Copyright 2005 Massachusetts Institute of Technology

This file is part of MIT/GNU Scheme.

MIT/GNU Scheme is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

MIT/GNU Scheme is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with MIT/GNU Scheme; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

|#

;;;;      Flat Polynomial Form, for Commutative Rings

(declare (usual-integrations))

;;; An fpf is a sorted list of terms.  
;;;  Each term has exponents and a coefficient. 

(define (fpf:make terms)
  (cond ((null? terms) :zero)
	((and (null? (cdr terms))
	      (fpf:constant-term? (car terms)))
	 (fpf:coefficient (car terms)))
	(else
	 (cons '*fpf* terms))))

(define (explicit-fpf? x)
  (and (pair? x)
       (eq? (car x) '*fpf*)))

(define (fpf:arity fpf)
  (if (number? fpf)
      0
      (length (fpf:exponents (car (fpf:terms fpf))))))

(define (fpf:constant-term? term)
  (all-zeros? (fpf:exponents term)))

(define (all-zeros? exponents)
  (or (null? exponents)
      (and (fix:= 0 (car exponents))
	   (all-zeros? (cdr exponents)))))

(define (fpf:terms fpf)
  (if (number? fpf)
      (if (zero? fpf)
	  '()
	  (list (fpf:make-term '() fpf)))
      (cdr fpf)))


(define (fpf:make-term exponents coeff)
  (if (all-zeros? exponents)
      (cons '() coeff)
      (cons exponents coeff)))

(define (fpf:exponents term)
  (car term))

(define (fpf:coefficient term)
  (cdr term))

(define (fpf:same-exponents? fs1 fs2)
  (equal? fs1 fs2))


(define (fpf:>exponents? fs1 fs2)
  (and (pair? fs1)
       (or (not (pair? fs2))
	   (fix:> (car fs1) (car fs2))
	   (fpf:>exponents? (cdr fs1) (cdr fs2)))))
       
(define fpf:zero :zero)
(define fpf:one  :one)
(define fpf:-one :-one)

(define fpf:identity
  (fpf:make (list (fpf:make-term (list 1) :one))))

(define (fpf:new-variables n)
  (make-initialized-list n
    (lambda (i)
      (fpf:make (list (fpf:make-term
		       (make-initialized-list n
			 (lambda (j)
			   (if (fix:= i j) 1 0)))
		       :one))))))


;;; Operations +, -, *, /

(define (fpf:+ a1 a2)
  (fpf:make (fpf:add-terms (fpf:terms a1) (fpf:terms a2))))

(define (fpf:add-terms xlist ylist)
  (cond ((null? xlist) ylist)
	((null? ylist) xlist)
	(else
	 (let ((f1 (fpf:exponents (car xlist)))
	       (f2 (fpf:exponents (car ylist))))
	   (cond ((fpf:same-exponents? f1 f2)
		  (let ((ncoeff
			 (+ (fpf:coefficient (car xlist))
			    (fpf:coefficient (car ylist)))))
		    (if (zero? ncoeff)
			(fpf:add-terms (cdr xlist) (cdr ylist))
			(cons (fpf:make-term f1 ncoeff)
			      (fpf:add-terms (cdr xlist)
					     (cdr ylist))))))
		 ((fpf:>exponents? f1 f2)
		  (cons (car xlist)
			(fpf:add-terms (cdr xlist) ylist)))
		 (else
		  (cons (car ylist)
			(fpf:add-terms xlist (cdr ylist)))))))))

(define (fpf:- minuend subtrahend)
  (fpf:+ minuend (fpf:* fpf:-one subtrahend)))

(define (fpf:negate x)
  (fpf:* fpf:-one x))

(define (fpf:* m1 m2)
  (fpf:make (fpf:mul-terms (fpf:terms m1) (fpf:terms m2))))

(define (fpf:mul-terms xlist ylist)
  (if (null? xlist)
      '()
      (fpf:add-terms (fpf:term*terms (car xlist) ylist)
		     (fpf:mul-terms (cdr xlist) ylist))))

(define (fpf:term*terms term terms)
  (let ((exponents (fpf:exponents term))
	(coeff (fpf:coefficient term)))
    (let lp ((terms terms))
      (if (null? terms)
	  '()
	  (cons (fpf:make-term
		 (fpf:combine-exponents exponents
					(fpf:exponents (car terms)))
		 (* coeff (fpf:coefficient (car terms))))
		(lp (cdr terms)))))))

(define (fpf:combine-exponents exponents1 exponents2)
  (cond ((null? exponents1) exponents2)
	((null? exponents2) exponents1)
	(else
	 (map fix:+ exponents1 exponents2))))

(define (fpf:square p)
  (fpf:* p p))

(define (fpf:expt base exponent)
  (define (expt-iter x count answer)
    (if (int:zero? count)
	answer
	(if (even? count)
	    (expt-iter (fpf:square x) (int:quotient count 2) answer)
	    (expt-iter x (int:- count 1) (fpf:* x answer)))))
  (cond ((number? base) (expt base exponent))
	((not (exact-nonnegative-integer? exponent))
	 (error "No inverse (FPF:EXPT):" base exponent))
	((zero? exponent) :one)
	(else
	 (expt-iter base exponent :one))))

;;; Converting between flat polynomials and other kinds of expressions

(define (fpf:->expression p vars)
  (cond ((number? p) p)
	((explicit-fpf? p)
	 (a-reduce symb:+
		   (map (lambda (term)
			  (symb:* (fpf:coefficient term)
				  (a-reduce symb:*
					    (map (lambda (exponent var)
						   (symb:expt var exponent))
						 (fpf:exponents term)
						 vars))))
			(fpf:terms p))))
	(else
	 (error "Bad fpf -- ->EXPRESSION" p vars))))


(define (fpf:expression-> expr cont #!optional less?)
  ;; cont = (lambda (poly vars) ... )
  (let ((evars
	 (sort (list-difference (variables-in expr)
				fpf:operators-known)
		(if (default-object? less?) alphaless? less?))))
    (cont ((expression-walker
	    (pair-up evars
		     (fpf:new-variables (length evars))
		     fpf:operator-table))
	   expr)
	  evars)))


(define +$fpf (accumulation fpf:+ fpf:zero))
(define -$fpf (inverse-accumulation fpf:- fpf:+ fpf:negate fpf:zero))
(define *$fpf (accumulation fpf:* fpf:one))

(define fpf:operator-table
  `((+        ,+$fpf)
    (-        ,-$fpf)
    (*        ,*$fpf)
    (negate   ,fpf:negate)
    (square   ,fpf:square)
    (expt     ,fpf:expt)))

(define fpf:operators-known (map car fpf:operator-table))
