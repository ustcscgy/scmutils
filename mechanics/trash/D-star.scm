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

;;; Locally a path is described by a stream of derivatives.

(define (provably-zero? x)
  (zero? (simplify x)))

(define (function-provably-zero? f)
  (let ((a (literal-number (generate-uninterned-symbol 'a))))
    (provably-zero? (f a))))

(define (path->local path)
  (define (derivs f)
    (if (function-provably-zero? f)
	the-empty-stream
	(cons-stream f (derivs (D f)))))
  (lambda (time)
    (cons-stream time
      (map-stream (lambda (f) (f time))
		  (derivs path)))))

#|
(define (path->local path)
  (define derivs
    (cons-stream path
		 (map-stream D derivs)))
  (lambda (time)
    (cons-stream time
      (map-stream (lambda (f) (f time))
		  derivs))))
|#

(define (path->finite-local path n)
  (define derivs
    (cons-stream path
		 (map-stream D derivs)))
  (lambda (time)
    (cons-stream time
      (map-stream (lambda (f) (f time))
		  (list->stream (stream-head derivs (fix:+ n 1)))))))


;;; If we are given particular values...

(define (make-local t q . derivatives)
  (cons-stream t
	       (cons-stream q
			    (list->stream derivatives))))


(define (local->time local)
  (stream-ref local 0))

(define (local->coordinates local)
  (stream-ref local 1))

(define (local->velocities local)
  (stream-ref local 2))

(define (local->accelerations local)
  (stream-ref local 3))

(define (stream:subst val i s)
  (if (fix:= i 0)
      (cons-stream val (stream-cdr s))
      (cons-stream (stream-car s)
		   (stream:subst val (fix:- i 1) (stream-cdr s)))))

(define (partial-time F)
  (lambda (local)
    ((D (lambda (time)
	  (F (stream:subst time 0 local))))
     (stream-ref local 0))))

(define (partial-coordinates F)
  (lambda (local)
    ((D (lambda (coords)
	  (F (stream:subst coords 1 local))))
     (stream-ref local 1))))

(define (partial-velocities F)
  (lambda (local)
    ((D (lambda (velocities)
	  (F (stream:subst velocities 2 local))))
     (stream-ref local 2))))

(define ((partial-local n) F)
  (lambda (local)
    ((D (lambda (velocities)
	  (F (stream:subst velocities n local))))
     (stream-ref local n))))


;;; "total time derivative"

(define (D*-procedure F)
  (define (DF-on-path q)
    (D (compose F (path->local q))))
  (abstract-to-local-function DF-on-path))

(define D* (make-operator D*-procedure 'D*))

;;; We want to write a system with pull rather than push.
;;; How many derivatives should we include?  
;;; This depends on the number of derivatives we want to take...

(define ((abstract-to-local-function f) local)
  ((f (osculating-path local))
   (local->time local)))

(define (osculating-path local)
  ;; Here local must be a finite stream.
  (let ((t0 (stream-ref local 0))
	(q0 (stream-ref local 1))
	(derivatives (stream-tail local 2)))
    (define (the-path t)
      (let ((dt (- t t0)))
	(let loop ((to-go derivatives) (sum q0) (n 2) (dt^n/n! dt))
	  (if (null? to-go)
	      sum
	      (loop (stream-cdr to-go)
		    (+ sum (* (stream-car to-go) dt^n/n!))
		    (fix:+ n 1)
		    (/ (* dt^n/n! dt) n))))))
    the-path))

(define (Lagrange-equations Lagrangian)
  (- (D* (partial-velocities Lagrangian))
     (partial-coordinates Lagrangian)))


(define ((L-harmonic m k) local)
  (let ((q (local->coordinates local))
	(qdot (local->velocities local)))
    (- (* 1/2 m (square qdot))
       (* 1/2 k (square q)))))

#|
;;; So we can get Lagrange's equations

(print-expression
 ((Lagrange-equations (L-harmonic 'm 'k))
  (make-local 't 'x 'v 'a)))
(+ (* a m) (* k x))
;Unspecified return value


;;; Works for 2 degrees of freedom

(print-expression
 ((Lagrange-equations (L-harmonic 'm 'k))
  (make-local 't #(x y) #(vx vy) #(ax ay))))
(vector (+ (* ax m) (* k x)) (+ (* ay m) (* k y)))
;Unspecified return value



;;; Adding extra local components is harmless.

(print-expression
 ((Lagrange-equations (L-harmonic 'm 'k))
  (make-local 't 'x 'v 'a 'j)))
(+ (* a m) (* k x))
;Unspecified return value

;;; But watch out.  If not enuf local componenents
;;;  are specified the answer is wrong.

(print-expression
 ((Lagrange-equations (L-harmonic 'm 'k))
  (make-local 't 'x 'v)))
(* k x)
;Unspecified return value


|#

;;; exp(x) - initial segment of the power series for exp(x)
;;;  is a function that is zero and has zero derivatives at zero up to
;;;  some order.  All higher order derivatives are one.

(define ((contact n) x)
  (if (= n 0)
      (- (exp x) 1)
      (let lp ((i n) (val (/ x n)))
	(if (fix:= i 1)
	    (- (exp x) (+ val 1))
	    (lp (fix:- i 1)
		(/ (* x (+ 1 val)) (fix:- i 1)))))))

#|

(print-expression (((expt D 0) (contact 3)) 0))
0

(print-expression (((expt D 1) (contact 3)) 0))
0

(print-expression (((expt D 2) (contact 3)) 0))
0

(print-expression (((expt D 3) (contact 3)) 0))
0

(print-expression (((expt D 4) (contact 3)) 0))
1

(print-expression (((expt D 5) (contact 3)) 0))
1

(print-expression (((expt D 6) (contact 3)) 0))
1

|#

#|
;;; The following version of OSCULATING-PATH produces a nasty object
;;;  in the answer if the local stream is not big enough.  This uses
;;;  Jack's CONTACT function.

(define (osculating-path local)
  ;; Here local must be a finite stream.
  (let ((m (- (stream-length local) 2))	
	(t0 (stream-ref local 0))
	(q0 (stream-ref local 1))
	(derivatives (stream-tail local 2)))
    (if (negative? m)
	(error "Local description too small"))
    (let ((extra (* (contact m) q0 'error-in-osculating-path )))
      (define (the-path t)
	(let ((dt (- t t0)))
	  (let loop ((to-go derivatives) (sum q0) (n 2) (dt^n/n! dt))
	    (if (null? to-go)
		(+ sum (extra dt))
		(loop (stream-cdr to-go)
		      (+ sum (* (stream-car to-go) dt^n/n!))
		      (+ n 1)
		      (/ (* dt^n/n! dt) n))))))
      the-path)))
|#

#|
(print-expression
 ((Lagrange-equations (L-harmonic 'm 'k))
  ((path->finite-local (vector (literal-function 'x)
			       (literal-function 'y))
		       2)
   't)))
(vector (+ (* k (x t)) (* m (((expt D 2) x) t)))
	(+ (* k (y t)) (* m (((expt D 2) y) t))))




(print-expression
 ((Lagrange-equations (L-harmonic 'm 'k))
  ((path->finite-local (vector (literal-function 'x)
			       (literal-function 'y))
		       1)
   't)))
(vector (+ (* error-in-osculating-path m (x t)) (* k (x t)))
	(+ (* error-in-osculating-path m (y t)) (* k (y t))))


(define ((central-Lagrangian-polar m V) local)
  (let ((q (local->coordinates local))
        (qdot (local->velocities local)))
    (let ((r (vector-ref q 0))
          (phi (vector-ref q 1))
          (rdot (vector-ref qdot 0))
          (phidot (vector-ref qdot 1)))
      (- (* 1/2 m
           (+ (square rdot)
              (square (* r phidot))) )
         (V r)))))
;Value: central-Lagrangian-polar

(print-expression
 ((Lagrange-equations (central-Lagrangian-polar 'm (literal-function 'V)))
  ((path->finite-local
    (vector (literal-function 'r) (literal-function 'phi))
    2)
  't)))
(vector
 (+ (* m (((expt D 2) r) t))
    (* -1 m (r t) (expt ((D phi) t) 2))
    ((D V) (r t)))
 (+ (* 2 m ((D r) t) (r t) ((D phi) t))
    (* m (((expt D 2) phi) t) (expt (r t) 2))))
|#

(define (generalized-Lagrange-equations Lagrangian)
  (lambda (local)
    (let ((m (- (stream-length local) 2)))
      (let lp ((n m))
	(if (= n 1)
	    ((partial-coordinates Lagrangian) local)
	    (- (((expt D* (- n 1)) ((partial-local n) Lagrangian)) local)
	       (lp (- n 1))))))))

#|
;;; With provably-zero trick, and not using CONTACT... 
;;;   D* does not loop on simple problems.

(print-expression ((D* (D* (lambda (local) (local->coordinates local))))
		   (make-local 't 'x 'v 'a)))
a

(print-expression ((D* (D* (lambda (local)
			     (square (local->coordinates local)))))
		   (make-local 't 'x 'v 'a)))
(+ (* 2 a x) (* 2 (expt v 2)))
|#

