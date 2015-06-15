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

;;;; Primitive Generic Operation Declarations 

(declare (usual-integrations))

;;; Unary Operators 

(define (g:type x) (generic-apply 'type x))
(define (g:type-predicate x) (generic-apply 'type-predicate x))
(define (g:arity x) (generic-apply 'arity x))

;;; The default arity is #f.  Objects must define their arity.

(assign-operation 'arity       none?        any?)


(define (g:inexact? x)
  (if (number? x)
      (inexact? x)
      (generic-apply 'inexact? x)))

(define (g:zero-like x) (generic-apply-default :zero 'zero-like x))
(define (g:one-like x) (generic-apply-default :one 'one-like x))
(define (g:identity-like x) (generic-apply 'identity-like x))

;;; Generic tests are conservative.  
;;; They will return #f unless the answer is known true.

(define (g:zero? x)
  (if (number? x)
      (zero? x)
      (generic-predicate 'zero? x)))
(define (g:one? x)
  (if (number? x)
      (one? x)
      (generic-predicate 'one? x)))
(define (g:identity? x)
  (if (number? x)
      (one? x)
      (generic-predicate 'identity? x)))

(define (g:negate x)
  (if (number? x)
      (- x)
      (generic-apply 'negate x)))
(define (g:invert x)
  (if (number? x)
      (/ x)
      (generic-apply 'invert x)))

(define (g:square x)
  (if (number? x)
      (* x x)
      (generic-apply-default-operation
       (lambda (args) (g:* (car args) (car args)))
       'square x)))

(define (g:dot-product x y)
  (cond ((and (number? x) (number? y)) (* x y))
	(else (generic-apply 'dot-product x y))))

(define (g:sqrt x)
  (if (number? x)
      (sqrt x)
      (generic-apply 'sqrt x)))

(define (g:exp x)
  (if (number? x)
      (exp x)
      (generic-apply 'exp x)))
(define (g:log x)
  (if (number? x)
      (log x)
      (generic-apply 'log x)))

(define (g:sin x)
  (if (number? x)
      (sin x)
      (generic-apply 'sin x)))
(define (g:cos x)
  (if (number? x)
      (cos x)
      (generic-apply 'cos x)))

(define (g:asin x)
  (if (number? x)
      (asin x)
      (generic-apply 'asin x)))
(define (g:acos x)
  (if (number? x)
      (acos x)
      (generic-apply 'acos x)))

(define (g:sinh x)
  (if (number? x)
      (sinh x)
      (generic-apply 'sinh x)))
(define (g:cosh x)
  (if (number? x)
      (cosh x)
      (generic-apply 'cosh x)))

(define (g:abs x)
  (if (number? x)
      (abs x)
      (generic-apply 'abs x)))

(define (g:determinant x)
  (if (number? x)
      x
      (generic-apply 'determinant x)))

(define (g:trace x)
  (if (number? x)
      x
      (generic-apply 'trace x)))


;;; Duplicate of text in OPERATOR.SCM, except that the explicit type
;;; tag is here rather than the variable operator-type-tag.  This is
;;; necessary because of a problem of load order.

(define (make-operator p #!optional name subtype arity #!rest opts)
  (if (default-object? name) (set! name #f))
  (if (default-object? subtype) (set! subtype #f))
  (if (default-object? arity) (set! arity (procedure-arity p)))
  (make-apply-hook p `(*operator* ,subtype ,name ,arity ,@opts)))

#| ;;;Old
(define (make-operator p #!optional name subtype #!rest opts)
  (if (default-object? name) (set! name #f))
  (if (default-object? subtype) (set! subtype #f))
  (make-apply-hook p `(*operator* ,subtype ,name ,@opts)))

(define (make-operator p . opt)
  (make-apply-hook (lambda (x) (p x))
		   `(*operator* ,p . ,opt)))
|#

(define g:derivative
  (make-operator
   (lambda (f)
     (generic-apply 'partial-derivative f '()))
   'derivative))

;;; Binary Operators

(define (g:=:bin x y)
  (if (and (number? x) (number? y))
      (= x y)
      (generic-predicate '= x y)))

(define (g:+:bin x y)
  (cond ((and (number? x) (number? y)) (+ x y))
	((g:zero? x) y)
	((g:zero? y) x)
	(else (generic-apply '+ x y))))

(define (g:-:bin x y)
  (cond ((and (number? x) (number? y)) (- x y))
	((g:zero? y) x)
	(else (generic-apply '- x y))))

(define (g:*:bin x y)
  (cond ((and (number? x) (number? y)) (* x y))
	((exact-zero? x) (g:zero-like y))
	((exact-zero? y) (g:zero-like x))
	((g:one? x) y)
	((g:one? y) x)
	(else (generic-apply '* x y))))

;;; In g:*:bin we test for exact (numerical) zero 
;;; because it is possible to produce a wrong-type 
;;; zero here, as follows:

;;;		  |0|             |0|
;;;	  |a b c| |0|   |0|       |0|
;;;	  |d e f| |0| = |0|, not  |0|

;;; We are less worried about the zero? below,
;;; because any invertible matrix is square.

(define (g:/:bin x y)
  (cond ((and (number? x) (number? y)) (/ x y))
	((g:zero? x) (g:zero-like y))
	((g:one? y) x)
	(else (generic-apply '/ x y))))

(define (g:expt x y)
  (cond ((and (number? x) (number? y)) (expt x y))
	;;((g:zero? x) x) ;No! consider 0^{-1}
	((g:one? x) x)
	((g:zero? y) (g:one-like x))
	((g:one? y) x)
	(else (generic-apply 'expt x y))))

(define (g:gcd:bin x y)
  (if (and (number? x) (number? y))
      (gcd x y)
      (generic-apply 'gcd x y)))

;;; Complex Operators

(define (g:make-rectangular real imag)
  (generic-apply 'make-rectangular real imag))

(define (g:make-polar mag ang)
  (generic-apply 'make-polar mag ang))

(define (g:real-part z) (generic-apply 'real-part z))
(define (g:imag-part z) (generic-apply 'imag-part z))
(define (g:magnitude z) (generic-apply 'magnitude z))
(define (g:angle z) (generic-apply 'angle z))

(define (g:conjugate x) (generic-apply 'conjugate x))


;;; Weird operators

(define (g:atan y #!optional x)
  (if (default-object? x)
      (g:atan1 y)
      (g:atan2 y x)))

(define (g:atan1 y) (generic-apply 'atan1 y))
(define (g:atan2 y x) (generic-apply 'atan2 y x))

(define (g:partial-derivative f . varspecs)
  (generic-apply 'partial-derivative f varspecs))

(define (g:partial . varspecs)
  (make-operator
   (lambda (f)
     (generic-apply 'partial-derivative f varspecs))
   `(partial ,@varspecs)))


(define *enable-literal-apply* #f)

(define (g:apply f . apply-args)
  (define (collapse l)
    (if (null? (cdr l))
	(car l)
	(cons (car l)
	      (collapse (cdr l)))))
  (if (null? apply-args)
      (error "No argument list for G:APPLY")
      (let ((args (collapse apply-args)))
	(cond ((procedure? f)
	       (apply f args))
	      ((and (symbol? f) *enable-literal-apply*)
	       (apply
		(literal-function f
				  (default-function-type (length args)))
		args))
	      (else
	       (generic-apply 'apply f args))))))

(define (with-literal-apply-enabled thunk)
  (fluid-let ((*enable-literal-apply* #t))
    (thunk)))

;;; N-ary Operator extensions

(define (g:= . args)
  (g:=:n args))

(define (g:=:n args)
  (cond ((null? args) #t)
	((null? (cdr args)) #t)
	(else
	 (let lp ((args (cddr args))
		  (larg (cadr args))
		  (ans (g:=:bin (car args) (cadr args))))
	   (if (null? args)
	       ans
	       (lp (cdr args)
		   (car args)
		   (and ans (g:=:bin larg (car args)))))))))


(define (g:+ . args)
  (g:+:n args))

(define (g:+:n args)
  (cond ((null? args) :zero)
	((null? (cdr args)) (car args))
	(else
	 (let lp ((args (cddr args))
		  (ans (g:+:bin (car args) (cadr args))))
	   (if (null? args)
	       ans
	       (lp (cdr args)
		   (g:+:bin ans (car args))))))))

(define (g:* . args)
  (g:*:n args))

(define (g:*:n args)
  (cond ((null? args) :one)
	((null? (cdr args)) (car args))
	(else
	 (let lp ((args (cddr args))
		  (ans (g:*:bin (car args) (cadr args))))
	   (if (null? args)
	       ans
	       (lp (cdr args)
		   (g:*:bin ans (car args))))))))

(define (g:- . args)
  (g:-:n args))

(define (g:-:n args)
  (cond ((null? args) :zero)
	((null? (cdr args)) (g:negate (car args)))
	(else
	 (g:-:bin (car args)
		  (g:+:n (cdr args))))))

(define (g:/ . args)
  (g:/:n args))

(define (g:/:n args)
  (cond ((null? args) :one)
	((null? (cdr args)) (g:invert (car args)))
	(else
	 (g:/:bin (car args)
		  (g:*:n (cdr args))))))


(define (g:gcd . args)
  (g:gcd:n args))

(define (g:gcd:n args)
  (cond ((null? args) :zero)
	((null? (cdr args)) (car args))
	(else
	 (let lp
	     ((as (cddr args))
	      (ans (g:gcd:bin (car args) (cadr args))))
	   (cond ((null? as) ans)
		 ((g:one? ans) ans)
		 (else
		  (lp (cdr as) (g:gcd:bin ans (car as)))))))))

