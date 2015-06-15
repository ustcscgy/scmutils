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

;;; GIGO:  Operators assume that arguments are in canonical form.
;;; Here, any object that is not either a sum, a product, or a power,
;;;  is considered a literal.

;;; An FPF is    a term, or
;;;              a sum of terms.
;;; A  term is   a number, or
;;;              a factor, or
;;;              a product of factors, or
;;;              a product of a number and factors.
;;; A  factor is a non-fpf, or
;;;              a positive integer power of a non-fpf.

#|
;;; Defined in Scheme

number?

;;; Defined in NUMSYMB

operator
operands
sum?
symb:+
product?
symb:*
difference?
symb:-
quotient?
symb:/
expt?
symb:expt
|#

(define (fpf? expr)
  (or (fpf:term? expr) (sum? expr)))

(define (fpf:term? expr)
  (or (number? expr) (fpf:factor? expr) (product? expr)))

(define (fpf:factor? expr)
  (or (non-fpf? expr) (fpf:expt? expr)))

(define (non-fpf? expr)
  (not (or (number? expr)
	   (fpf:expt? expr)
	   (product? expr)
	   (sum? expr))))

(define (fpf:expt? x)
  (and (expt? x) (exact-integer? (cadr (operands x)))))

(define (fpf:terms p)
  (cond ((sum? p) (operands p))
	(else (list p))))

(define (fpf:make terms)  
  (cond ((null? terms) :zero)
	((null? (cdr terms)) (car terms))
	(else (a-reduce symb:+ terms))))


(define (fpf:make-term factors coeff)
  (cond ((null? factors) coeff)
	((one? coeff)
	 (a-reduce symb:* factors))
	(else
	 (symb:* coeff (a-reduce symb:* factors)))))

(define (fpf:coefficient p)
  (cond ((number? p) p)
	((product? p)
	 (if (number? (car (operands p)))
	     (car (operands p))
	     :one))
	(else :one)))

(define (fpf:factors p)
  (cond ((number? p) '())
	((product? p)
	 (if (number? (car (operands p)))
	     (cdr (operands p))
	     (operands p)))
	(else (list p))))

(define (fpf:same-factors? fs1 fs2)
  (equal? fs1 fs2))

(define (fpf:<factors? fs1 fs2)
  (cond ((null? fs1) (not (null? fs2)))
	((null? fs2) #f)
	(else
	 (or (fix:< (length fs1) (length fs2))
	     (fpf:<base? (fpf:base (car fs1))
			 (fpf:base (car fs2)))
	     (and (fpf:same-base? (fpf:base (car fs1))
				  (fpf:base (car fs2)))
		  (fpf:<factors? (cdr fs1) (cdr fs2)))))))

(define (fpf:make-factor exponent base)
  (symb:expt base exponent))

(define (fpf:base factor)
  (cond ((fpf:expt? factor) (car (operands factor)))
	(else factor)))

(define (fpf:exponent factor)
  (cond ((fpf:expt? factor) (cadr (operands factor)))
	(else :one)))


(define (fpf:same-base? b1 b2)
  (equal? b1 b2))

(define (fpf:<base? b1 b2)
  (expr:< b1 b2))


;;; Operations +, -, *, /

(define (fpf:+ a1 a2)
  (fpf:make (fpf:add-terms (fpf:terms a1) (fpf:terms a2))))

(define (fpf:add-terms xlist ylist)
  (cond ((null? xlist) ylist)
	((null? ylist) xlist)
	(else
	 (let ((f1 (fpf:factors (car xlist)))
	       (f2 (fpf:factors (car ylist))))
	   (cond ((fpf:same-factors? f1 f2)
		  (let ((ncoeff
			 (+ (fpf:coefficient (car xlist))
			    (fpf:coefficient (car ylist)))))
		    (if (zero? ncoeff)
			(fpf:add-terms (cdr xlist) (cdr ylist))
			(cons (fpf:make-term f1 ncoeff)
			      (fpf:add-terms (cdr xlist)
					     (cdr ylist))))))
		 ((fpf:<factors? f1 f2)
		  (cons (car xlist)
			(fpf:add-terms (cdr xlist) ylist)))
		 (else
		  (cons (car ylist)
			(fpf:add-terms xlist (cdr ylist)))))))))

(define (fpf:- minuend subtrahend)
  (fpf:+ minuend (fpf:* :-one subtrahend)))

(define (fpf:negate x)
  (fpf:* :-one x))

(define (fpf:* m1 m2)
  (fpf:make (fpf:mul-terms (fpf:terms m1) (fpf:terms m2))))

(define (fpf:mul-terms xlist ylist)
  (if (null? xlist)
      '()
      (fpf:add-terms (fpf:term*terms (car xlist) ylist)
		     (fpf:mul-terms (cdr xlist) ylist))))

(define (fpf:term*terms term terms)
  (let ((factors (fpf:factors term))
	(coeff (fpf:coefficient term)))
    (let lp ((terms terms))
      (if (null? terms)
	  '()
	  (cons (fpf:make-term
		 (fpf:combine-factors factors (fpf:factors (car terms)))
		 (* coeff (fpf:coefficient (car terms))))
		(lp (cdr terms)))))))

(define (fpf:combine-factors factors1 factors2)
  (cond ((null? factors1) factors2)
	((null? factors2) factors1)
	(else
	 (let ((b1 (fpf:base (car factors1)))
	       (b2 (fpf:base (car factors2))))
	   (cond ((fpf:same-base? b1 b2)
		  (let ((nexpt
			 (+ (fpf:exponent (car factors1))
			    (fpf:exponent (car factors2)))))
		    (cond ((zero? nexpt)
			   (fpf:combine-factors (cdr factors1)
						(cdr factors2)))
			  (else
			   (cons (fpf:make-factor nexpt b1)
				 (fpf:combine-factors (cdr factors1)
						      (cdr factors2)))))))
		 ((fpf:<base? b1 b2)
		  (cons (car factors1)
			(fpf:combine-factors (cdr factors1) factors2)))
		 (else
		  (cons (car factors2)
			(fpf:combine-factors factors1 (cdr factors2)))))))))

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
  p)

(define (fpf:expression-> expr cont #!optional less?)
  ;; cont = (lambda (poly vars) ... )
  (let ((evars
	 (sort (list-difference (variables-in expr)
				fpf:operators-known)
		(if (default-object? less?) alphaless? less?))))
    (cont ((expression-walker
	    (pair-up evars
		     evars
		     fpf:operator-table))
	   expr)
	  evars)))


(define +$fpf (accumulation fpf:+ :zero))
(define -$fpf (inverse-accumulation fpf:- fpf:+ fpf:negate :zero))
(define *$fpf (accumulation fpf:* :one))

(define fpf:operator-table
  `((+        ,+$fpf)
    (-        ,-$fpf)
    (*        ,*$fpf)
    (negate   ,fpf:negate)
    (square   ,fpf:square)
    (expt     ,fpf:expt)))

(define fpf:operators-known (map car fpf:operator-table))
