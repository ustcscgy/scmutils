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

;;; Jack's Extend idea is quite nice... (f t q qd)
#|
(define (Dt f)
  (lambda (t q qd qdd)
    (+ (((partial 0) f) t q qd)
       (* qd (((partial 1) f) t q qd))
       (* qdd (((partial 2) f) t q qd)))))

; Dt (partial_qdot L) - partial_q L = 0

(define (extend f)
  (lambda (t q qd qdd)
    (f t q qd)))

(define ((L k m) t q qd)
  (- (* 1/2 m (square qd))
     (* 1/2 k (square q))))

(define (LE L)
  (- (Dt ((partial 2) L))
     (extend ((partial 1) L))))

#|

(print-expression
 ((LE (L 'k  'm)) 't 'x 'v 'a))
(+ (* a m) (* k x))
|#
|#

;;; Dt is an operator on functions that take local descriptions.
;;;  (F t q qd qdd ...)
;;;  It produces new functions of local descriptions.

(define ((Dt F) t q . derivs)
  (+ (apply ((partial 0) F) t q derivs)
     (apply +
	    (map (lambda (deriv index)
		   (* (apply ((partial (fix:+ index 1)) F) t q derivs)
		      deriv))
		 derivs
		 (iota (length derivs))))))


;;; Problem here is "state" is restricted to t,q,qd by construction
;;; in lag.scm.   Need to loosen this to make things really work.

(define ((extend F) t q qdot . derivs)
  (if (vector? q)
      (F (vector-append (vector t) q qdot))
      (F (vector t q qdot))))

(define (LE L)
  (let ((LL (extend L)))
    (- (Dt ((partial 2) LL))
       ((partial 1) LL))))

#|
(define ((L-harmonic m k) state)
  (let ((q (state->q state))
	(qdot (state->qdot state)))
    (- (* 1/2 m (square qdot))
       (* 1/2 k (square q)))))

(print-expression
 ((LE (L-harmonic 'm 'k))
  't 'x 'xd 'xdd))
(+ (* k x) (* m xdd))
|#
#|
(define ((central-Lagrangian-polar m V) state)
  (let ((q (state->q state))
        (qdot (state->qdot state)))
    (let ((r (vector-ref q 0))
          (phi (vector-ref q 1))
          (rdot (vector-ref qdot 0))
          (phidot (vector-ref qdot 1)))
      (- (* 1/2 m
           (+ (square rdot)
              (square (* r phidot))) )
         (V r)))))

(print-expression
 ((LE (central-Lagrangian-polar 'm (literal-function 'V)))
  't #(r phi) #(rdot phidot) #(rdd phidd)))
(vector (+ (* -1 m (expt phidot 2) r) (* m rdd) ((D V) r))
        (+ (* m phidd (expt r 2)) (* 2 m phidot r rdot)))
|#

