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

;;; for chap6 continuous transformations examples

(define ((Hdo alpha omega omega0) state)
  (let ((p (state->p state))
	(t (state->t state))
	(q (state->q state)))
    (+ (* 1/2 (square p)) 
       (* 1/2 (square omega0) (square q))
       (* -1 alpha q (cos (* omega t))))))

(define (make-advancer alpha omega omega0)
  (let ((alphap (/ alpha (- (square omega0) (square omega)))))
    (lambda (delta-t)
      (lambda (state)
	(let ((t (state->t state))
	      (q (state->q state))
	      (p (state->p state)))
	  (let ((z1 (- q (* alphap (cos (* omega t)))))
		(z2 (/ (+ p (* alphap omega (sin (* omega t)))) omega0))
		(u1 (* alphap (cos (* omega (+ t delta-t)))))
		(u2 (* -1 alphap (/ omega omega0) (sin (* omega (+ t delta-t)))))
		(cosw0dt (cos (* omega0 delta-t)))
		(sinw0dt (sin (* omega0 delta-t))))
	    (->H-state (+ (state->t state) delta-t)
		       (+ (* cosw0dt z1) (* sinw0dt z2) u1)
		       (* omega0 (+ (* -1 sinw0dt z1) (* cosw0dt z2) u2)))))))))
(define ->C make-advancer)

(define (((solution advancer) state0) t)
  ((advancer (- t (state->t state0))) state0))


#|

(print-expression
 (simplify
  (let ((state0 (->H-state 't 'q 'p)))
    (((solution (->C 'alpha 'omega 'omega0)) state0) 't))))
;(vector t q p)


(print-expression
 (simplify
  (let ((state0 (->H-state 't 'q 'p)))
    ((derivative ((solution (->C 'alpha 'omega 'omega0)) state0)) 't))))
;(vector 1 p (+ (* alpha (cos (* omega t))) (* -1 (expt omega0 2) q)))

(print-expression
 (let ((state0 (->H-state 't0 'q0 'p0)))
     (((Hamilton-equations (Hdo 'alpha 'omega 'omega0)) 
       (compose state->q ((solution (->C 'alpha 'omega 'omega0)) state0))
       (compose state->p ((solution (->C 'alpha 'omega 'omega0)) state0))) 
      't)))
;(vector 0 0)

(print-expression
 (let ((state0 (->H-state 't0 'q0 'p0)))
   (let ((sigma ((solution (->C 'alpha 'omega 'omega0)) state0)))
     ((- (D (compose ((->C 'alpha 'omega 'omega0) 'delta-t) sigma))
	(+ (compose (partial_t ((->C 'alpha 'omega 'omega0) 'delta-t)) sigma)
	   (compose (Poisson-bracket ((->C 'alpha 'omega 'omega0) 'delta-t)
				     (Hdo 'alpha 'omega 'omega0))
		    sigma))) 
      't))))
(vector 0 0 0)

(print-expression
 (->poisson-form
  (expression
   (let ((state0 (->H-state 't0 'q0 'p0)))
     (let ((sigma ((solution (->C 'alpha 'omega 'omega0)) state0)))
       ((compose (partial_t ((->C 'alpha 'omega 'omega0) 'delta-t)) sigma) 
	't))))))
; a big mess

but (((C delta-t) sigma) t) = (sigma (+ t delta-t))

(print-expression
 (->poisson-form
  (expression
   (let ((state0 (->H-state 't0 'q0 'p0)))
     (let ((sigma ((solution (->C 'alpha 'omega 'omega0)) state0)))
       (- (((->C 'alpha 'omega 'omega0) 'delta-t) (sigma 't))
	  (sigma (+ 't 'delta-t))))))))
(vector 0 0 0)

so it should satisfy Hamilton's equations in Poisson bracket form
D(Q circ sigma) = {Q, H} circ sigma

or equally
D(Q circ sigma') = {Q, H} circ sigma'

where sigma' is C circ sigma
D(Q circ C circ sigma) 
  = { Q circ C, H } circ sigma
???
  = {Q, H} circ (C sigma)

so apparently 
{C, H} circ sigma + partialt C circ sigma = {Q, H} circ (C sigma)

check this
(print-expression
 (let ((state0 (->H-state 't0 'q0 'p0)))
   (let ((sigma ((solution (->C 'alpha 'omega 'omega0)) state0)))
     ((- (compose 
	  (Poisson-bracket (lambda (state) state) (Hdo 'alpha 'omega 'omega0))
	  (compose ((->C 'alpha 'omega 'omega0) 'delta-t) sigma))
	 (+ (compose (partial_t ((->C 'alpha 'omega 'omega0) 'delta-t)) sigma)
	    (compose (Poisson-bracket ((->C 'alpha 'omega 'omega0) 'delta-t)
				      (Hdo 'alpha 'omega 'omega0))
		     sigma))) 
      't))))
(vector -1 0 0)
not quite?

|#

(define ((shift-t delta-t) state)
  (->H-state (+ (state->t state) delta-t)
             (state->q state)
             (state->p state)))


#|

;;; check symplectic

(print-expression
 (symplectic-transform? ((->C 'alpha 'omega 'omega0) 'delta-t) (->H-state 't0 'q0 'p0)))
(matrix-by-rows (list 0 0) (list 0 0))

so it is symplectic

now Cp
(define ((advance->Cp delta-advance) delta-t) 
  (compose (shift-t (- delta-t)) (delta-advance delta-t)))

(print-expression
 (symplectic-transform? (compose (shift-t (- 'delta-t))
				 ((->C 'alpha 'omega 'omega0) 'delta-t))
			(->H-state 't0 'q0 'p0)))
(matrix-by-rows (list 0 0) (list 0 0))
so that is symplectic also


|#

#|

canonical tests

(named-lambda (real-test-canonical? C H Hprime)
  (- (compose (phase-space-derivative H) C)
     (* (D C) (phase-space-derivative Hprime))))

(print-expression
 ((real-test-canonical? ((->C 'alpha 'omega 'omega0) 'delta-t)
			(Hdo 'alpha 'omega 'omega0)
			(Hdo 'alpha 'omega 'omega0))
  (->H-state 't0 'q0 'p0)))
(vector 0 0 0)

ok


(print-expression
 ((real-test-canonical? (compose (shift-t (- 'delta-t))
				 ((->C 'alpha 'omega 'omega0) 'delta-t))
			(compose (Hdo 'alpha 'omega 'omega0)
				 (shift-t 'delta-t))
			(Hdo 'alpha 'omega 'omega0))
  (->H-state 't0 'q0 'p0)))
(vector 0 0 0)
ok


but 

(print-expression
 ((real-test-canonical? (compose (shift-t (- 'delta-t))
				 ((->C 'alpha 'omega 'omega0) 'delta-t))
			(Hdo 'alpha 'omega 'omega0)
			(compose (Hdo 'alpha 'omega 'omega0)
				 (shift-t (- 'delta-t))))
  (->H-state 't0 'q0 'p0)))
; mess

|#


;;;----------------------------------------------------------------
;;; now try to do this with the z = C(t, z') convention

;;; presumably this means
;;; (t, q(t), p(t)) = (Cpp D)(t+D, q(t+D), p(t+D))
;;; (t+D, q(t+D), p(t+D)) = advancer(D)(t, q(t), p(t))

(define ((advancer->Cpp advancer) delta-t) (advancer (- delta-t)))

;;; and the version that does not affect time
;;; (t+D, q(t+D), p(t+D)) = advancer(D)(t, q(t), p(t))
;;; (t, q(t), p(t)) = (Cppp D)(t, q(t+D), p(t+D)) takes an inconsistent state!

(define ((advancer->Cppp advancer) delta-t)
  (compose (advancer (- delta-t)) (shift-t delta-t)))

#|
symplectic 

(print-expression
 (symplectic-transform? 
  ((advancer->Cpp (make-advancer 'alpha 'omega 'omega0))
   'delta-t)
  (->H-state 't0 'q0 'p0)))
;(matrix-by-rows (list 0 0) (list 0 0))

(print-expression
 (symplectic-transform? 
  ((advancer->Cppp (make-advancer 'alpha 'omega 'omega0))
   'delta-t)
  (->H-state 't0 'q0 'p0)))
;(matrix-by-rows (list 0 0) (list 0 0))

|#

#|

(print-expression
 ((real-test-canonical? 
   ((advancer->Cpp (make-advancer 'alpha 'omega 'omega0)) 
    'delta-t)
   (Hdo 'alpha 'omega 'omega0)
   (Hdo 'alpha 'omega 'omega0))
  (->H-state 't0 'q0 'p0)))
;(vector 0 0 0)

(print-expression
 ((real-test-canonical? 
   ((advancer->Cppp (make-advancer 'alpha 'omega 'omega0)) 
    'delta-t)
   (Hdo 'alpha 'omega 'omega0)
   (compose (Hdo 'alpha 'omega 'omega0) (shift-t 'delta-t)))
  (->H-state 't0 'q0 'p0)))
;(vector 0 0 0)


|#


#|

(print-expression
 ((- (compose (Hdo 'alpha 'omega 'omega0)
	      ((advancer->Cpp (make-advancer 'alpha 'omega 'omega0)) 
	       'delta-t))
     (Hdo 'alpha 'omega 'omega0))
  (->H-state 't0 'q0 'p0)))
; mess presumably H' != H circ C

|#