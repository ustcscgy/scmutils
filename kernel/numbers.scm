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

;;;; Generic Numerical Arithmetic

(declare (usual-integrations))

(define (n:type m) number-type-tag)
(define (n:type-predicate m) numerical-quantity?)

;;; If this variable is set to #f then numbers are not allowed
;;; to be applied to arguments.  Otherwise they are considered 
;;; constant functions of any number of arguments.

(define *numbers-are-constant-functions* #f)

(define (n:arity m) *at-least-zero*)

(define (n:self n . args)
  (if *numbers-are-constant-functions*
      n
      (error "Application of a number not allowed" n args)))

(define (n:deriv n . args) :zero)

;;; We rename all of the standard numerical functions to make
;;;  it possible to reliably access them, even in a generic
;;;  environment.  See NUMERIC.SCM for other such renamings.

(define n:inexact? inexact?)

(define (n:zero-like x)
  (if (exact? x) 0 0.0))

(define (n:one-like x)
  (if (exact? x) 1 1.0))

(define n:zero? zero?)
(define n:one? one?)

(define n:negate negate)
(define n:invert invert)

(define (n:sqrt x)
  (if (exact-rational? x)
      (/ (sqrt (numerator x))
	 (sqrt (denominator x)))
      (sqrt x)))

(define n:exp exp)
(define n:log log)

(define n:sin sin)
(define n:cos cos)
(define n:tan tan)
(define n:sec sec)
(define n:csc csc)

(define n:asin asin)
(define n:acos acos)
(define n:atan atan)

(define n:sinh sinh)
(define n:cosh cosh)
(define n:tanh tanh)
(define n:sech sech)
(define n:csch csch)

(define n:= =)

(define n:+ +)
(define n:- -)
(define n:* *)
(define n:/ /)

(define n:expt expt)
(define n:gcd gcd)

(define n:make-rectangular make-rectangular)
(define n:make-polar make-polar)
(define n:real-part real-part)
(define n:imag-part imag-part)
(define n:magnitude magnitude)
(define n:angle angle)

(define n:conjugate conjugate)

(define n:log10 log10)
(define n:log2 log2)
(define n:exp10 exp10)
(define n:exp2 exp2)

(define n:square square)
(define n:cube cube)

(define n:sigma sigma)

;;; Here we assign the primitive numerical generic operators.

(assign-operation 'type             n:type             number?)
(assign-operation 'type-predicate   n:type-predicate   number?)
(assign-operation 'arity            n:arity            number?)

(assign-operation 'inexact?         n:inexact?         number?)

(assign-operation 'zero-like        n:zero-like        number?)
(assign-operation 'one-like         n:one-like         number?)
(assign-operation 'identity-like    n:one-like         number?)

(assign-operation 'zero?            n:zero?            number?)
(assign-operation 'one?             n:one?             number?)

(assign-operation 'negate           n:negate           number?)
(assign-operation 'invert           n:invert           number?)

(assign-operation 'sqrt             n:sqrt             number?)

(assign-operation 'exp              n:exp              number?)
(assign-operation 'log              n:log              number?)

(assign-operation 'sin              n:sin              number?)
(assign-operation 'cos              n:cos              number?)

(assign-operation 'asin             n:asin             number?)
(assign-operation 'acos             n:acos             number?)

(assign-operation 'sinh             n:sinh             number?)
(assign-operation 'cosh             n:cosh             number?)

(assign-operation '=          n:=            number? number?)

(assign-operation '+          n:+            number? number?)
(assign-operation '-          n:-            number? number?)
(assign-operation '*          n:*            number? number?)
(assign-operation '/          n:/            number? number?)

(assign-operation 'dot-product          n:*            number? number?)

(assign-operation 'expt       n:expt         number? number?)
(assign-operation 'gcd        n:gcd          number? number?)

(assign-operation 'make-rectangular    n:make-rectangular real? real?)
(assign-operation 'make-polar          n:make-polar       real? real?)
(assign-operation 'real-part           n:real-part        number?)
(assign-operation 'imag-part           n:imag-part        number?)
(assign-operation 'magnitude           n:magnitude        number?)
(assign-operation 'angle               n:angle            number?)

(assign-operation 'conjugate           n:conjugate        number?)

(assign-operation 'atan1               n:atan             number?)
(assign-operation 'atan2               n:atan             number? number?)

;(assign-operation 'partial-derivative  n:deriv            number? any?)

(assign-operation 'apply               n:self             number? any?)

;;;; Abstract numbers generalize numerical quantities.

(define (literal-number symbol)
  (make-literal number-type-tag symbol)) 

(define (make-numerical-combination operator #!optional reverse?)
  (if (default-object? reverse?)
      (lambda operands 
	(make-numsymb-expression operator operands))
      (lambda operands 
	(make-numsymb-expression operator (reverse operands)))))

(define (an:zero-like n) :zero)
(define (an:one-like n) :one)

(assign-operation 'type             n:type             abstract-number?)
(assign-operation 'type-predicate   n:type-predicate   abstract-number?)
(assign-operation 'arity            n:arity            abstract-number?)

(assign-operation 'inexact?   (has-property? 'inexact)         abstract-number?)

(assign-operation 'zero-like        an:zero-like        abstract-number?)
(assign-operation 'one-like         an:one-like         abstract-number?)
(assign-operation 'identity-like    an:one-like         abstract-number?)

(assign-operation 'negate  (make-numerical-combination 'negate)  abstract-number?)
(assign-operation 'invert  (make-numerical-combination 'invert)  abstract-number?)

(assign-operation 'sqrt    (make-numerical-combination 'sqrt)    abstract-number?)

(assign-operation 'exp     (make-numerical-combination 'exp)     abstract-number?)
(assign-operation 'log     (make-numerical-combination 'log)     abstract-number?)

(assign-operation 'sin     (make-numerical-combination 'sin)     abstract-number?)
(assign-operation 'cos     (make-numerical-combination 'cos)     abstract-number?)

(assign-operation 'asin    (make-numerical-combination 'asin)    abstract-number?)
(assign-operation 'acos    (make-numerical-combination 'acos)    abstract-number?)

(assign-operation 'sinh    (make-numerical-combination 'sinh)    abstract-number?)
(assign-operation 'cosh    (make-numerical-combination 'cosh)    abstract-number?)

;(assign-operation 'derivative       n:deriv            abstract-number?)

(assign-operation '+   (make-numerical-combination '+)     abstract-number? abstract-number?)
(assign-operation '+   (make-numerical-combination '+)     number?          abstract-number?)
(assign-operation '+   (make-numerical-combination '+ 'r)  abstract-number? number?)

(assign-operation '-   (make-numerical-combination '-)     abstract-number? abstract-number?)
(assign-operation '-   (make-numerical-combination '-)     number?          abstract-number?)
(assign-operation '-   (make-numerical-combination '-)     abstract-number? number?)

(assign-operation '*   (make-numerical-combination '*)     abstract-number? abstract-number?)
(assign-operation '*   (make-numerical-combination '*)     number?          abstract-number?)
(assign-operation '*   (make-numerical-combination '* 'r)  abstract-number? number?)

(assign-operation 'dot-product   (make-numerical-combination '*)     abstract-number? abstract-number?)
(assign-operation 'dot-product   (make-numerical-combination '*)     number?          abstract-number?)
(assign-operation 'dot-product   (make-numerical-combination '* 'r)  abstract-number? number?)

(assign-operation '/   (make-numerical-combination '/)     abstract-number? abstract-number?)
(assign-operation '/   (make-numerical-combination '/)     number?          abstract-number?)
(assign-operation '/   (make-numerical-combination '/)     abstract-number? number?)

(assign-operation 'expt (make-numerical-combination 'expt) abstract-number? abstract-number?)
(assign-operation 'expt (make-numerical-combination 'expt) number?          abstract-number?)
(assign-operation 'expt (make-numerical-combination 'expt) abstract-number? number?)

(assign-operation 'gcd (make-numerical-combination 'gcd)   abstract-number? abstract-number?)
(assign-operation 'gcd (make-numerical-combination 'gcd)   number?          abstract-number?)
(assign-operation 'gcd (make-numerical-combination 'gcd 'r) abstract-number? number?)

(assign-operation 'make-rectangular
		  (make-numerical-combination 'make-rectangular)
		  abstract-number? abstract-number?)
(assign-operation 'make-rectangular
		  (make-numerical-combination 'make-rectangular)
		  number? abstract-number?)
(assign-operation 'make-rectangular
		  (make-numerical-combination 'make-rectangular)
		  abstract-number? number?)

(assign-operation 'make-polar
		  (make-numerical-combination 'make-polar)
		  abstract-number? abstract-number?)
(assign-operation 'make-polar
		  (make-numerical-combination 'make-polar)
		  number? abstract-number?)
(assign-operation 'make-polar
		  (make-numerical-combination 'make-polar)
		  abstract-number? number?)

(assign-operation 'real-part (make-numerical-combination 'real-part) abstract-number?)
(assign-operation 'imag-part (make-numerical-combination 'imag-part) abstract-number?)
(assign-operation 'magnitude (make-numerical-combination 'magnitude) abstract-number?)
(assign-operation 'angle     (make-numerical-combination 'angle)     abstract-number?)

(assign-operation 'conjugate (make-numerical-combination 'conjugate) abstract-number?)

(assign-operation 'atan1     (make-numerical-combination 'atan)      abstract-number?)
(assign-operation 'atan2     (make-numerical-combination 'atan)
		                                    abstract-number? abstract-number?)
(assign-operation 'atan2     (make-numerical-combination 'atan) number? abstract-number?)
(assign-operation 'atan2     (make-numerical-combination 'atan) abstract-number? number?)

;(assign-operation 'partial-derivative  n:deriv            abstract-number? any?)
(assign-operation 'apply               n:self             abstract-number? any?)


;;; Conservative tests...  These tests will return TRUE only if the
;;; default simplifier can prove that the answer is TRUE.  This is a
;;; bad idea to make generic, because it puts a simplification burden
;;; on lots of parts of the system. -- GJS

(define (an:= x y)
  (exact-zero? (simplify (g:- x y))))

(define (an:zero? x)
  (exact-zero? (simplify x)))

(define (an:one? x)
  (exact-one? (simplify x)))

;;; Sigh.
;;; (assign-operation '=              an:=            abstract-number? abstract-number?)
;;; (assign-operation 'zero?          an:zero?        abstract-number?)
;;; (assign-operation 'one?           an:one?         abstract-number?)
