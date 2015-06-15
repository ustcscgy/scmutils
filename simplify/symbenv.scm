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

;;;; Symbolic environment for simplification

(declare (usual-integrations))

(define (symbolic-operator operator-symbol)
  (let ((v (hash-table/get symbolic-operator-table operator-symbol #f)))
    (if v
	v
	(error "Undefined symbolic operator" operator-symbol))))


;;; N-ary Operator extensions

(define symbolic:=:bin (symbolic-operator '=))

(define (symbolic:=:n args)
  (cond ((null? args) #t)
	((null? (cdr args)) #t)
	(else
	 (let lp ((args (cddr args))
		  (larg (cadr args))
		  (ans (symbolic:=:bin (car args) (cadr args))))
	   (if (null? args)
	       ans
	       (lp (cdr args)
		   (car args)
		   (and ans (symbolic:=:bin larg (car args)))))))))


(define symbolic:+:bin (symbolic-operator '+))

(define (symbolic:+:n args)
  (cond ((null? args) :zero)
	((null? (cdr args)) (car args))
	(else
	 (let lp ((args (cddr args))
		  (ans (symbolic:+:bin (car args) (cadr args))))
	   (if (null? args)
	       ans
	       (lp (cdr args)
		   (symbolic:+:bin ans (car args))))))))

(define symbolic:*:bin (symbolic-operator '*))

(define (symbolic:*:n args)
  (cond ((null? args) :one)
	((null? (cdr args)) (car args))
	(else
	 (let lp ((args (cddr args))
		  (ans (symbolic:*:bin (car args) (cadr args))))
	   (if (null? args)
	       ans
	       (lp (cdr args)
		   (symbolic:*:bin ans (car args))))))))


(define symbolic:-:bin (symbolic-operator '-))
(define symbolic:negate (symbolic-operator 'negate))

(define (symbolic:-:n args)
  (cond ((null? args) :zero)
	((null? (cdr args)) (symbolic:negate (car args)))
	(else
	 (symbolic:-:bin (car args)
			 (symbolic:+:n (cdr args))))))


(define symbolic:/:bin (symbolic-operator '/))
(define symbolic:invert (symbolic-operator 'invert))

(define (symbolic:/:n args)
  (cond ((null? args) :one)
	((null? (cdr args)) (symbolic:invert (car args)))
	(else
	 (symbolic:/:bin (car args)
			 (symbolic:*:n (cdr args))))))


(define (symbolic-environment-maker)
  (let ((e (extend-ic-environment scmutils-base-environment)))
    (let ((d
	   (lambda (name value)
	     (local-assignment e name value))))

      (d '*environment* 'symbolic-environment)

      ;; Unary operators from generic.scm
      #|
      (d 'type (symbolic-operator 'type))
      (d 'type-predicate (symbolic-operator 'type-predicate))
      (d 'arity (symbolic-operator 'arity))

      (d 'inexact? (symbolic-operator 'inexact?))

      (d 'zero-like (symbolic-operator 'zero-like))
      (d 'one-like (symbolic-operator 'one-like))
      (d 'identity-like (symbolic-operator 'identity-like))
      |#
      (d 'zero? (symbolic-operator 'zero?))
      (d 'one? (symbolic-operator 'one?))
      ;;	(d 'identity? (symbolic-operator 'identity?))

      (d 'negate (symbolic-operator 'negate))
      (d 'invert (symbolic-operator 'invert))

      (d 'square (symbolic-operator 'square))
      (d 'cube   (symbolic-operator 'cube))

      (d 'sqrt (symbolic-operator 'sqrt))

      (d 'exp (symbolic-operator 'exp))
      (d 'log (symbolic-operator 'log))
      #|
      (d 'exp2  (symbolic-operator 'exp2))
      (d 'exp10 (symbolic-operator 'exp10))
      (d 'log2  (symbolic-operator 'log2))
      (d 'log10 (symbolic-operator 'log10))
      |#
      (d 'sin (symbolic-operator 'sin))
      (d 'cos (symbolic-operator 'cos))
      (d 'tan (symbolic-operator 'tan))
      (d 'sec (symbolic-operator 'sec))
      (d 'csc (symbolic-operator 'csc))

      (d 'asin (symbolic-operator 'asin))
      (d 'acos (symbolic-operator 'acos))

      (d 'sinh (symbolic-operator 'sinh))
      (d 'cosh (symbolic-operator 'cosh))
      #|
      (d 'tanh (symbolic-operator 'tanh))
      (d 'sech (symbolic-operator 'sech))
      (d 'csch (symbolic-operator 'csch))
      |#
      (d 'abs (symbolic-operator 'abs))

      ;; (d 'derivative (symbolic-operator 'derivative))
      
      ;; Binary (and nary) operators from generic.scm

      (d 'expt (symbolic-operator 'expt))
      ;; (d 'gcd (symbolic-operator 'gcd))


      ;; Complex operators from generic.scm

      (d 'make-rectangular (symbolic-operator 'make-rectangular))
      (d 'make-polar (symbolic-operator 'make-polar))

      (d 'real-part (symbolic-operator 'real-part))
      (d 'imag-part (symbolic-operator 'imag-part))
      (d 'magnitude (symbolic-operator 'magnitude))
      (d 'angle (symbolic-operator 'angle))

      (d 'conjugate (symbolic-operator 'conjugate))


      ;; Wierd operators from generic.scm

      (d 'atan (symbolic-operator 'atan))
      #|
      (d 'partial-derivative (symbolic-operator 'partial-derivative))
      (d 'partial (symbolic-operator 'partial))

      (d 'apply (symbolic-operator 'apply))


      ;; Compound operators from mathutil.scm

      (d 'arg-scale (symbolic-operator 'arg-scale))
      (d 'arg-shift (symbolic-operator 'arg-shift))

      (d 'sigma (symbolic-operator 'sigma))

      (d 'compose (symbolic-operator 'compose))
      |#


      (d '= (named-lambda (= . args) (symbolic:=:n args)))

      (d '+ (named-lambda (+ . args) (symbolic:+:n args)))

      (d '* (named-lambda (* . args) (symbolic:*:n args)))

      (d '- (named-lambda (- . args) (symbolic:-:n args)))

      (d '/ (named-lambda (/ . args) (symbolic:/:n args))))
    e))

(define symbolic-environment (symbolic-environment-maker))
