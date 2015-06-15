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

;;;; Rule systems for simplification

;;; Default is simplifier lives dangerously.

;;; Allows (log (exp x)) => x 
;;;  Can confuse x=(x0+n*2pi)i with x0

(define log-exp-simplify? true)

;;; Allows (sqrt (square x)) => x
;;;  Not good if x is negative.
(define sqrt-expt-simplify? true)

;;; Allows (asin (sin x)) => x, etc
;;;  Loses multivalue info, as in log-exp
(define inverse-simplify? true)

;;; Wierd case: ((D magnitude) (square x)) => 1
(define ignore-zero? true)

;;; Allows commutation of partial derivatives.
;;;  Only OK if components selected by partials are unstructured (e.g. Real)
(define commute-partials? true)


;;; However, we have control over the defaults

(define (log-exp-simplify doit?)
  (assert (boolean? doit?) "Argument must be a boolean.")
  (clear-memoizer-tables)
  (set! log-exp-simplify? doit?))

(define (sqrt-expt-simplify doit?)
  (assert (boolean? doit?) "Argument must be a boolean.")
  (clear-memoizer-tables)
  (set! sqrt-expt-simplify? doit?))

(define (inverse-simplify doit?)
  (assert (boolean? doit?) "Argument must be a boolean.")
  (clear-memoizer-tables)
  (set! inverse-simplify? doit?))

(define (ignore-zero-simplify doit?)
  (assert (boolean? doit?) "Argument must be a boolean.")
  (clear-memoizer-tables)
  (set! ignore-zero? doit?))

(define (commute-partials-simplify doit?)
  (assert (boolean? doit?) "Argument must be a boolean.")
  (clear-memoizer-tables)
  (set! commute-partials? doit?))


(define (non-integer? x)
  (not (integer? x)))

(define (even-integer? x)
  (and (integer? x) (even? x) (fix:> x 1)))

(define (odd-integer? x)
  (and (integer? x) (odd? x) (fix:> x 1)))

(define universal-reductions
  (rule-system
   ( (EXP (* (? n integer?) (LOG (? x))))
     none
     (EXPT (: x) (: n)) )
     
   ( (EXP (LOG (? x))) none (: x) )
   ( (LOG (EXP (? x))) log-exp-simplify? (: x) )

   ( (EXPT (SQRT (? x)) (? n even-integer?))
     none
     (EXPT (: x) (: (quotient n 2))) )
   ( (SQRT (EXPT (? x) (? n even-integer?)))
     sqrt-expt-simplify?
     (EXPT (: x) (: (quotient n 2))) )

   ( (EXPT (SQRT (? x)) (? n odd-integer?))
     none
     (* (SQRT (: x)) (EXPT (: x) (: (quotient (fix:- n 1) 2)))) )
   ( (SQRT (EXPT (? x) (? n odd-integer?)))
     sqrt-expt-simplify?
     (* (SQRT (: x)) (EXPT (: x) (: (quotient (fix:- n 1) 2)))) )

   ( (SQRT (EXP (? x)))
     sqrt-expt-simplify?
     (EXP (/ (: x) 2)) )

   ( (LOG (SQRT (? x)))
     none
     (* 1/2 (LOG (: x))) )

   
   ( (/ (? x) (SQRT (? x)))
     none
     (SQRT (: x)) )

   ( (/ (SQRT (? x)) (? x))
     none
     (/ 1 (SQRT (: x))) )

   ( (/ (* (?? u) (? x) (?? v)) (SQRT (? x)))
     none
     (* (:: u) (SQRT (: x)) (:: v)) )

   ( (/ (* (?? u) (SQRT (? x)) (?? v)) (? x))
     none
     (/ (* (:: u) (:: v)) (SQRT (: x))) )

   ( (/ (? x) (* (?? u) (SQRT (? x)) (?? v)))
     none
     (/  (SQRT (: x)) (* (:: u) (:: v))) )

   ( (/ (SQRT (? x)) (* (?? u) (? x) (?? v)))
     none
     (/ 1 (* (:: u) (SQRT (: x)) (:: v))) )

   ( (/ (* (?? p) (? x) (?? q))
	(* (?? u) (SQRT (? x)) (?? v)))
     none
     (/ (* (:: p) (SQRT (: x)) (:: q))
	(* (:: u) (:: v))) )

   ( (/ (* (?? p) (SQRT (? x)) (?? q))
	(* (?? u) (? x) (?? v)))
     none
     (/ (* (:: p) (:: q))
	(* (:: u) (SQRT (: x)) (:: v))) )


   ( (expt (expt (? x) (? p/q)) (? m*q))
     (integer? (* p/q m*q))
     (expt (: x) (: (* p/q m*q))) )
     
   ( (SIN (ASIN (? x))) none (: x) )
   ( (COS (ACOS (? x))) none (: x) )
   ( (TAN (ATAN (? x))) none (: x) )
   ( (SIN (ACOS (? x))) none (SQRT (- 1 (EXPT (: x) 2))) )
   ( (COS (ASIN (? y))) none (SQRT (- 1 (EXPT (: y) 2))) )
   ( (TAN (ASIN (? y))) none (/ (: y) (SQRT (- 1 (EXPT (: y) 2)))) )
   ( (TAN (ACOS (? x))) none (/ (SQRT (- 1 (EXPT (: x) 2))) (: x)) )

   ;; sin atan, cos atan below

   ( (ASIN (SIN (? x))) inverse-simplify? (: x) )
   ( (ASIN (COS (? x))) inverse-simplify? (- pi/2 (: x)) )
   ( (ACOS (COS (? x))) inverse-simplify? (: x) )
   ( (ACOS (SIN (? x))) inverse-simplify? (- pi/2 (: x)) )
   ( (ATAN (TAN (? x))) inverse-simplify? (: x) )

   #|
   ( (ATAN (? y) (? x)) none (ATAN (/ (: y) (: x))) )

   ( (ATAN (/ (SIN (? x)) (COS (? x)))) inverse-simplify? (: x) )

   ( (SIN (ATAN (/ (? a) (? b))))
     none
     (/ (: a) (SQRT (+ (EXPT (: a) 2) (EXPT (: b) 2)))) )

    ( (COS (ATAN (/ (? a) (? b))))
     none
     (/ (: b) (SQRT (+ (EXPT (: a) 2) (EXPT (: b) 2)))) )

   ( (SIN (ATAN (? a)))
     none
     (/ (: a) (SQRT (+ 1 (EXPT (: a) 2)))) )

   ( (COS (ATAN (? a)))
     none
     (/ 1 (SQRT (+ 1 (EXPT (: a) 2)))) )
   |#

   ( (ATAN (/ (? y) (? x))) none (ATAN (: y) (: x)) )

   ( (ATAN (? y)) none (ATAN (: y) 1) )

   ( (ATAN (SIN (? x)) (COS (? x))) inverse-simplify? (: x) )

   ( (SIN (ATAN (? a) (? b)))
     none
     (/ (: a) (SQRT (+ (EXPT (: a) 2) (EXPT (: b) 2)))) )

    ( (COS (ATAN (? a) (? b)))
     none
     (/ (: b) (SQRT (+ (EXPT (: a) 2) (EXPT (: b) 2)))) )


   ( (MAGNITUDE (* (? x) (? y) (?? ys)))
     none
     (* (MAGNITUDE (: x)) (MAGNITUDE (* (: y) (:: ys)))) )

   ( (MAGNITUDE (EXPT (? x) (? n even-integer?)))
     none
     (EXPT (: x) (: n)) )

   ( ((DERIVATIVE MAGNITUDE) (EXPT (? x) (? n even-integer?)))
     ignore-zero?
     1 )
   ))

(define sqrt-expand
  (rule-system
		     
   ( (SQRT (* (? x) (? y)))
     none
     (* (SQRT (: x)) (SQRT (: y))) )

   ( (SQRT (* (? x) (? y) (?? ys)))
     none
     (* (SQRT (: x)) (SQRT (* (: y) (:: ys)))) )

   ( (SQRT (/ (? x) (? y)))
     none
     (/ (SQRT (: x)) (SQRT (: y))) )

   ( (SQRT (/ (? x) (? y) (?? ys)))
     none
     (/ (SQRT (: x)) (SQRT (* (: y) (:: ys)))) )
   ))

(define sqrt-contract
  (rule-system
   ( (* (?? a) (SQRT (? x)) (?? b) (SQRT (? y)) (?? c))
     none
     (* (:: a) (:: b) (:: c) (SQRT (* (: x) (: y)))) )
   ))

(define specfun->logexp
  (rule-system
   ( (SQRT (? x)) none (EXP (* 1/2 (LOG (: x)))) )

   ( (ATAN (? z))
     none
     (/ (- (LOG (+ 1 (* +i (: z)))) (LOG (- 1 (* +i (: z))))) +2i) )

   ( (ASIN (? z))
     none
     (* -i (LOG (+ (* +i (: z)) (SQRT (- 1 (EXPT (: Z) 2)))))) )

   ( (ACOS (? z))
     none
     (* -i (LOG (+ (: z) (* +i (SQRT (- 1 (EXPT (: Z) 2))))))) )

   ( (SINH (? u)) none (/ (- (EXP (: u)) (EXP (* -1 (: u)))) 2) )

   ( (COSH (? u)) none (/ (+ (EXP (: u)) (EXP (* -1 (: u)))) 2) )

   ( (EXPT (? x) (? y non-integer?)) none (EXP (* (: y) (LOG (: x)))) )
   ))

(define logexp->specfun
  (rule-system
     ( (EXP (* -1 (LOG (? x)))) none (EXPT (: x) -1) )

     ( (EXP (* 1/2 (LOG (? x1)))) none (SQRT (: x1)) )

     ( (EXP (* -1/2 (LOG (? x1)))) none (/ 1 (SQRT (: x1))) )

     ( (EXP (* 3/2 (LOG (? x1)))) none (EXPT (SQRT (: x1)) 3) )

     ( (EXP (* -3/2 (LOG (? x1)))) none (EXPT (SQRT (: x1)) -3) )
     ))

(define log-contract
  (rule-system
   ( (+ (?? x1) (LOG (? x2)) (?? x3) (LOG (? x4)) (?? x5))
     none
     (+ (:: x1) (:: x3) (:: x5) (LOG (* (: x2) (: x4)))) )

   ( (* (? n integer?) (?? f1) (LOG (? x)) (?? f2))
     none
     (* (:: f1) (LOG (EXPT (: x) (: n))) (:: f2)) )

   ( (+ (?? x1)
	(* (?? f1) (LOG (? x)) (?? f2))
	(?? x2)
	(* (?? f3) (LOG (? y)) (?? f4))
	(?? x3))
     (exact-zero?
      (rcf:simplify
       (- (apply * (append f1 f2))
	  (apply * (append f3 f4)))))
     (+ (:: x1)
	(:: x2)
	(:: x3)
	(* (:: f1) (LOG (* (: x) (: y))) (:: f2))) )
   ))

(define log-expand
  (rule-system
   ( (LOG (* (? x1) (? x2) (?? xs)))
     none
     (+ (LOG (: x1)) (LOG (* (: x2) (:: xs)))) )

   ( (LOG (EXPT (? x) (? e))) none (* (: e) (LOG (: x))) )
   ))


(define (list< l1 l2)
  (cond ((null? l1) (not (null? l2)))
	((null? l2) #f)
	((< (car l1) (car l2)) #t)
	((> (car l1) (car l2)) #f)
	(else (list< (cdr l1) (cdr l2)))))

(define canonicalize-partials
  (rule-system
   ( ( (partial (?? i)) ( (partial (?? j)) (? f) ) )
     (and commute-partials? (list< j i))
     ( (partial (:: j)) ( (partial (:: i)) (: f) ) ) )))


;;;; Trigonometry

;;; The following rules are used to convert all trig expressions to
;;; ones involving only SIN and COS functions.

(define trig->sincos
  (rule-system
   ( (TAN (? x)) none (/ (SIN (: x)) (COS (: x))) )

   ( (COT (? x)) none (/ (COS (: x)) (SIN (: x))) )

   ( (SEC (? x)) none (/ 1 (COS (: x))) )

   ( (CSC (? x)) none (/ 1 (SIN (: x))) )
   ))


;;; Sometimes we want to express combinations of SIN and COS in terms
;;; of other functions.

(define sincos->trig
  (rule-system
   ( (/ (SIN (? x)) (COS (? x))) none (TAN (: x)) )

   ( (/ (* (?? n1) (SIN (? x)) (?? n2)) (COS (? x)))
     none
     (* (:: n1) (TAN (: x)) (:: n2)) )

     
   ( (/ (SIN (? x)) (* (?? d1) (COS (? x)) (?? d2)))
     none
     (/ (TAN (: x)) (* (:: d1) (:: d2))) )
     

   ( (/ (* (?? n1) (SIN (? x)) (?? n2))
	(* (?? d1) (COS (? x)) (?? d2)))
     none
     (/ (* (:: n1) (TAN (: x)) (:: n2))
	(* (:: d1) (:: d2))) )

;   ( (/ (COS (? x)) (SIN (? x))) none (COT (: x)) )

;   ( (/ (* (?? n1) (COS (? x)) (?? n2)) (SIN (? x)))
;     none
;     (* (:: n1) (COT (: x)) (:: n2)) )

     
;   ( (/ (COS (? x)) (* (?? d1) (SIN (? x)) (?? d2)))
;     none
;     (/ (COT (: x)) (* (:: d1) (:: d2))) )
     
;   ( (/ (* (?? n1) (COS (? x)) (?? n2))
;	(* (?? d1) (SIN (? x)) (?? d2)))
;     none
;     (/ (* (:: n1) (COT (: x)) (:: n2))
;	(* (:: d1) (:: d2))) )
   ))

;;; SIN is odd, and COS is even.  We canonicalize by moving the sign
;;; out of the first term of the argument.

(define angular-parity
  (rule-system
   ( (COS (? n negative-number?))
     none
     (COS (: (- n))) )

   ( (COS (* (? n negative-number?) (?? x)))
     none
     (COS (* (: (- n)) (:: x))) )

   ( (COS (+ (* (? n negative-number?) (?? x)) (?? y)))
     none
     (COS (- (* (: (- n)) (:: x)) (:: y))) )

   ( (SIN (? n negative-number?))
     none
     (- (SIN (: (- n)))) )

   ( (SIN (* (? n negative-number?) (?? x)))
     none
     (- (SIN (* (: (- n)) (:: x)))) )

   ( (SIN (+ (* (? n negative-number?) (?? x)) (?? y)))
     none
     (- (SIN (- (* (: (- n)) (:: x)) (:: y)))) )
   ))

(define expand-multiangle
  (rule-system
   ( (SIN (* 2 (? x)))
     none
     (* 2 (SIN (: x)) (COS (: x))) )

   ( (COS (* 2 (? x)))
     none
     (- (* 2 (EXPT (COS (: x)) 2)) 1) )

   ( (SIN (* (? n exact-integer?) (? f) (?? fs))) ;at least one f
     (> n 1)
     (+ (* (SIN (* (: (- n 1)) (: f) (:: fs))) (COS (* (: f) (:: fs))))
	(* (COS (* (: (- n 1)) (: f) (:: fs))) (SIN (* (: f) (:: fs))))) )

   ( (SIN (+ (? x) (? y) (?? ys)))	;at least one y
     none
     (+ (* (SIN (: x)) (COS (* (: y) (:: ys))))
	(* (COS (: x)) (SIN (* (: y) (:: ys))))) )

   ( (COS (* (? n exact-integer?) (? f) (?? fs))) ;at least one f
     (> n 1)
     (- (* (COS (* (: (- n 1)) (: f) (:: fs))) (COS (* (: f) (:: fs))))
	(* (SIN (* (: (- n 1)) (: f) (:: fs))) (SIN (* (: f) (:: fs))))) )

   ( (COS (+ (? x) (? y) (?? ys)))	;at least one y
     none
     (- (* (COS (: x)) (COS (* (: y) (:: ys))))
	(* (SIN (: x)) (SIN (* (: y) (:: ys))))) )
   ))

(define contract-multiangle
  (rule-system
   ( (* (?? u) (SIN (? x)) (?? v) (SIN (? y)) (?? w))
     none
     (* 1/2 (- (COS (- (: x) (: y))) (COS (+ (: x) (: y)))) (:: u) (:: v) (:: w)) )

   ( (* (?? u) (COS (? x)) (?? v) (COS (? y)) (?? w))
     none
     (* 1/2 (+ (COS (- (: x) (: y))) (COS (+ (: x) (: y)))) (:: u) (:: v) (:: w)) )

   ( (* (?? u) (SIN (? x)) (?? v) (COS (? y)) (?? w))
     none
     (* 1/2 (+ (SIN (+ (: x) (: y))) (SIN (- (: x) (: y)))) (:: u) (:: v) (:: w)) )

   ( (* (?? u) (COS (? y)) (?? v) (SIN (? x)) (?? w))
     none
     (* 1/2 (+ (SIN (+ (: x) (: y))) (SIN (- (: x) (: y)))) (:: u) (:: v) (:: w)) )
   ))

(define contract-expt-trig
  (rule-system
   ( (EXPT (SIN (? x)) (? n exact-integer?))
     (> n 1)
     (* 1/2 (- 1 (COS (* 2 (: x)))) (EXPT (SIN (: x)) (: (- n 2)))))

   ( (EXPT (COS (? x)) (? n exact-integer?))
     (> n 1)
     (* 1/2 (+ 1 (COS (* 2 (: x)))) (EXPT (COS (: x)) (: (- n 2)))))
   ))


;;; SINCOS.SCM has code for sin^2 x + cos^2 x => 1,
;;; however, sometimes we want a few other rules to help:

(define sin^2->cos^2
  (rule-system
   ( (expt (sin (? x)) (? n at-least-two?))
     none
     (evaluate
      (* (expt (sin (: x)) (- (: n) 2))
	 (- 1 (expt (cos (: x)) 2)))) )
   ))


(define cos^2->sin^2
  (rule-system
   ( (expt (cos (? x)) (? n at-least-two?))
     none
     (evaluate
      (* (expt (cos (: x)) (- (: n) 2))
	 (- 1 (expt (sin (: x)) 2)))) )
   ))

;;; Here are some residual rules.

(define sincos-random
  (rule-system

   ( (+ (?? a1) (? a) (?? a2) (EXPT (COS (? x)) (? n at-least-two?)) (?? a3))
     (exact-zero? (rcf:simplify (+ a (expt (cos x) (- n 2)))))
     (+ (:: a1) (:: a2) (:: a3) (* (EXPT (SIN (: x)) 2) (: a))) )

   ( (+ (?? a1) (EXPT (COS (? x)) (? n at-least-two?)) (?? a2) (? a) (?? a3))
     (exact-zero? (rcf:simplify (+ a (expt (cos x) (- n 2)))))
     (+ (:: a1) (:: a2) (:: a3) (* (EXPT (SIN (: x)) 2) (: a))) )

   ( (+ (?? a1) (? a) (?? a2) (EXPT (SIN (? x)) (? n at-least-two?)) (?? a3))
     (exact-zero? (rcf:simplify (+ a (expt (sin x) (- n 2)))))
     (+ (:: a1) (:: a2) (:: a3) (* (EXPT (COS (: x)) 2) (: a))) )

   ( (+ (?? a1) (EXPT (SIN (? x)) (? n at-least-two?)) (?? a2) (? a) (?? a3))
     (exact-zero? (rcf:simplify (+ a (expt (sin x) (- n 2)))))
     (+ (:: a1) (:: a2) (:: a3) (* (EXPT (COS (: x)) 2) (: a))) )

   ( (+ (?? a1)
	(? a)
	(?? a2)
	(* (?? b1) (EXPT (COS (? x)) (? n at-least-two?)) (?? b2))
	(?? a3))
     (exact-zero?
      (rcf:simplify
       (+ (apply * (append b1 b2 (list (expt (cos x) (- n 2))))) a)))
     (+ (:: a1) (:: a2) (:: a3) (* (: a) (EXPT (SIN (: x)) 2))) )
     
   ( (+ (?? a1)
	(? a)
	(?? a2)
	(* (?? b1) (EXPT (SIN (? x)) (? n at-least-two?)) (?? b2))
	(?? a3))
     (exact-zero?
      (rcf:simplify
       (+ (apply * (append b1 b2 (list (expt (sin x) (- n 2))))) a)))
     (+ (:: a1) (:: a2) (:: a3) (* (: a) (EXPT (COS (: x)) 2))) )


   ( (+ (?? a1)
	(* (?? b1) (EXPT (COS (? x)) (? n at-least-two?)) (?? b2))
	(?? a2)
	(? a)
	(?? a3))
     (exact-zero?
      (rcf:simplify
       (+ (apply * (append b1 b2 (list (expt (cos x) (- n 2))))) a)))
     (+ (:: a1) (:: a2) (:: a3) (* (: a) (EXPT (SIN (: x)) 2))) )
     
   ( (+ (?? a1)
	(* (?? b1) (EXPT (SIN (? x)) (? n at-least-two?)) (?? b2))
	(?? a2)
	(? a)
	(?? a3))
     (exact-zero?
      (rcf:simplify
       (+ (apply * (append b1 b2 (list (expt (sin x) (- n 2))))) a)))
     (+ (:: a1) (:: a2) (:: a3) (* (: a) (EXPT (COS (: x)) 2))) )

   ))

;;; We can eliminate SIN and COS in favor of complex exponentials 

(define sincos->exp1
  (rule-system
   ( (SIN (? x)) 
     none
     (/ (- (EXP (* +i (: x))) (EXP (* -i (: x))))
	+2i) )
     
   ( (COS (? x)) 
     none
     (/ (+ (EXP (* +i (: x))) (EXP (* -i (: x))))
	2) )
   ))

(define sincos->exp2
  (rule-system
   ( (SIN (? x)) 
     none
     (/ (- (EXP (* +i (: x))) (/ 1 (EXP (* +i (: x)))))
	+2i) )
     
   ( (COS (? x)) 
     none
     (/ (+ (EXP (* +i (: x))) (/ 1 (EXP (* +i (: x)))))
	2) )
   ))


;;; Under favorable conditions, we can replace the trig functions.

(define exp->sincos
  (rule-system
     ( (EXP (? c1 imaginary-number?))
       (positive? (imag-part c1))
       (+ (COS (: (n:imag-part c1)))
	  (* +i (SIN (: (imag-part c1))))) )

     ( (EXP (? c1 imaginary-number?))
       (negative? (n:imag-part c1))
       (+ (COS (: (- (imag-part c1))))
	  (* -i (SIN (: (- (imag-part c1)))))) )

     ( (EXP (* (? c1 imaginary-number?) (?? f)))
       (positive? (n:imag-part c1))
       (+ (COS (* (: (imag-part c1)) (:: f)))
	  (* +i (SIN (* (: (imag-part c1)) (:: f))))) )

     ( (EXP (* (? c1 imaginary-number?) (?? f)))
       (negative? (n:imag-part c1))
       (+ (COS (* (: (- (imag-part c1))) (:: f)))
	  (* -i (SIN (* (: (- (imag-part c1))) (:: f))))) )
     ))

;;; The following predicates are used in trig rules.

(define (negative-number? x)
  (and (number? x) (negative? x)))

(define (complex-number? z)
  (and (complex? z)
       (not (n:zero? (n:real-part z)))
       (not (n:zero? (n:imag-part z)))))

(define (imaginary-number? z)
  (and (complex? z)
       (not (n:zero? z))
       (n:zero? (n:real-part z))))

(define (imaginary-integer? z)
  (and (complex? z)
       (not (n:zero? z))
       (zero? (n:real-part z))
       (exact-integer? (n:imag-part z))))

	    
(define exp-contract
  (rule-system
   ( (* (?? x1) (EXP (? x2)) (?? x3) (EXP (? x4)) (?? x5))
     none
     (* (:: x1) (:: x3) (:: x5) (EXP (+ (: x2) (: x4)))) )

   ( (EXPT (EXP (? x)) (? p)) none (EXP (* (: p) (: x))) )
   ))


(define exp-expand
  (rule-system
   ( (EXP (+ (? x1) (? x2) (?? xs)))
     none
     (* (EXP (: x1)) (EXP (+ (: x2) (:: xs)))) )

   ( (EXP (* (? x imaginary-integer?) (?? factors)))
     (> (n:imag-part x) 1)
     (EXPT (EXP (* +i (:: factors))) (: (n:imag-part x))) )

   ( (EXP (* (? x imaginary-integer?) (?? factors)))
     (< (n:imag-part x) -1)
     (EXPT (EXP (* -i (:: factors))) (: (- (n:imag-part x)))) )

   ( (EXP (* (? n exact-integer?) (?? factors)))
     (> n 1)
     (EXPT (EXP (* (:: factors))) (: n)) )

   ( (EXP (* (? n exact-integer?) (?? factors)))
     (< n -1)
     (EXPT (EXP (* -1 (:: factors))) (: (- n))) )

   ( (EXP (? x complex-number?))
     none
     (* (EXP (: (n:real-part x)))
	(EXP (: (n:* (n:imag-part x) +i)))) )

   ( (EXP (* (? x complex-number?) (?? factors)))
     none
     (* (EXP (* (: (n:real-part x)) (:: factors)))
	(EXP (* (: (n:* (n:imag-part x) +i)) (:: factors)))) )
   ))

(define exp-general
  (rule-system
   ( (EXP (- (? x1)))
     none
     (/ 1 (EXP (: x1))) )

   ( (EXP (- (? x1) (? x2)))
     none
     (/ (EXP (: x1)) (EXP (: x2))) )

   ))

;;;; Simplifiers defined using these rule sets

;;; Assuming that expression comes in canonical it goes out canonical

(define (simplify-until-stable rule-simplify canonicalize)
  (define (simp exp)
    (let ((newexp (rule-simplify exp)))
      (if (equal? exp newexp)
	  exp
	  (simp (canonicalize newexp)))))
  simp)


;;; The usual canonicalizer is

(define simplify-and-flatten
  (compose fpf:simplify rcf:simplify))

(define ->poisson-form
  (compose simplify-and-flatten
	   angular-parity
	   (simplify-until-stable contract-multiangle simplify-and-flatten)
	   (simplify-until-stable contract-expt-trig simplify-and-flatten)
	   simplify-and-flatten
	   trig->sincos))

(define (trigexpand exp)
  ((compose simplify-and-flatten
	    sincos->trig
	    (simplify-until-stable sincos-flush-ones simplify-and-flatten)
	    simplify-and-flatten
	    exp->sincos
	    (simplify-until-stable exp-expand simplify-and-flatten)	
	    (simplify-until-stable exp-contract simplify-and-flatten)
	    (simplify-until-stable exp-expand simplify-and-flatten)
	    simplify-and-flatten
	    sincos->exp1
	    trig->sincos)
   exp))

(define (trigcontract exp)
  ((compose simplify-and-flatten
	    sincos->trig
	    (simplify-until-stable sincos-flush-ones simplify-and-flatten)
	    simplify-and-flatten
	    exp->sincos
	    (simplify-until-stable exp-expand simplify-and-flatten)
	    (simplify-until-stable sincos-flush-ones simplify-and-flatten)
	    simplify-and-flatten
	    exp->sincos
	    (simplify-until-stable exp-contract simplify-and-flatten)
	    (simplify-until-stable exp-expand simplify-and-flatten)
	    simplify-and-flatten
	    sincos->exp1
	    trig->sincos)
   exp))

(define (full-simplify exp)
  ((compose rcf:simplify

	    (simplify-until-stable universal-reductions
				   rcf:simplify)
;	    (simplify-until-stable sqrt-contract
;				   rcf:simplify)
	    (simplify-until-stable sqrt-expand
				   rcf:simplify)
	    (simplify-until-stable sqrt-contract
				   rcf:simplify)
	    rcf:simplify
	    logexp->specfun
	    sincos->trig	    
	    (simplify-until-stable sincos-flush-ones
				   rcf:simplify)
	    rcf:simplify
	    exp->sincos
	    (simplify-until-stable (compose log-expand exp-expand)
				   rcf:simplify)	
	    (simplify-until-stable (compose log-contract exp-contract)
				   rcf:simplify)
	    (simplify-until-stable (compose log-expand exp-expand)
				   rcf:simplify)
	    rcf:simplify
	    sincos->exp1
	    trig->sincos
	    specfun->logexp
	    (simplify-until-stable (compose universal-reductions sqrt-expand)
				   rcf:simplify)
	    rcf:simplify
	    )
   exp))

(define (oe-simplify exp)
  ((compose (simplify-until-stable universal-reductions
				   simplify-and-flatten)
	    (simplify-until-stable sqrt-expand
				   simplify-and-flatten)
	    (simplify-until-stable sqrt-contract
				   simplify-and-flatten)
	    simplify-and-flatten
	    sincos->trig	   
	    (simplify-until-stable sincos-random
				   simplify-and-flatten)


	    simplify-and-flatten
	    sin^2->cos^2
	    (simplify-until-stable sincos-flush-ones
				   simplify-and-flatten)

	    (simplify-until-stable (compose log-expand exp-expand)
				   simplify-and-flatten)	
	    (simplify-until-stable (compose log-contract exp-contract)
				   simplify-and-flatten)
	    (simplify-until-stable (compose log-expand exp-expand)
				   simplify-and-flatten)
	    (simplify-until-stable angular-parity
				   simplify-and-flatten)
	    (simplify-until-stable (compose universal-reductions sqrt-expand)
				   simplify-and-flatten)
	    simplify-and-flatten
	    trig->sincos
	    canonicalize-partials
	    )
   exp))

(define (easy-simplify exp)
  ((compose (simplify-until-stable (compose universal-reductions sqrt-expand)
				   simplify-and-flatten)
	    simplify-and-flatten
	    root-out-squares
	    (simplify-until-stable sqrt-contract
				   simplify-and-flatten)

	    sincos->trig
	    (simplify-until-stable sincos-random
				   simplify-and-flatten)
	    simplify-and-flatten
	    sin^2->cos^2
	    (simplify-until-stable sincos-flush-ones
				   simplify-and-flatten)

	    (simplify-until-stable (compose log-expand exp-expand)
				   simplify-and-flatten)	
	    (simplify-until-stable (compose log-contract exp-contract)
				   simplify-and-flatten)

	    (simplify-until-stable (compose universal-reductions
					    angular-parity
					    log-expand
					    exp-expand
					    sqrt-expand)
				   simplify-and-flatten)
	    simplify-and-flatten
	    trig->sincos
	    canonicalize-partials
	    )
   exp))


