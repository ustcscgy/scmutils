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


(define (((M-of-q->omega-of-t M-of-q) q) t)
  (define M-on-path (compose M-of-q q))
  (define (w-cross t)
    (* ((derivative M-on-path) t)
       (m:transpose (M-on-path t))))
  (antisymmetric->3vector-components (w-cross t)))

#|
(show-expression
 (((M-of-q->omega-of-t Euler->M)
   (vector (literal-function 'theta)
	   (literal-function 'phi)
	   (literal-function 'psi)))
  't))
(vector
 (+ (* (cos (phi t)) ((D theta) t))
    (* (sin (phi t)) ((D psi) t) (sin (theta t))))
 (+ (* -1 (cos (phi t)) ((D psi) t) (sin (theta t)))
    (* ((D theta) t) (sin (phi t))))
 (+ ((D phi) t) (* (cos (theta t)) ((D psi) t))))

|#

(define (((M-of-q->omega-body-of-t M-of-q) q) t)
  (* (m:transpose (M-of-q (q t)))
     (((M-of-q->omega-of-t M-of-q) q) t)))

#|
(show-expression
 (((M-of-q->omega-body-of-t Euler->M)
   (vector (literal-function 'theta)
	   (literal-function 'phi)
	   (literal-function 'psi)))
  't))
(vector
 (+ (* ((D phi) t) (sin (psi t)) (sin (theta t)))
    (* (cos (psi t)) ((D theta) t)))
 (+ (* ((D phi) t) (sin (theta t)) (cos (psi t)))
    (* -1 (sin (psi t)) ((D theta) t)))
 (+ (* (cos (theta t)) ((D phi) t)) ((D psi) t)))

|#

(define ((abstract-to-state-function f) state)
  (let ((t (state->t state))
        (q (state->q state))
        (qdot (state->qdot state)))
    (let ((osc-q (osculating-path t q qdot)))
      ((f osc-q) t))))

(define (M->omega M-of-q)
  (abstract-to-state-function 
    (M-of-q->omega-of-t M-of-q)))

(define (M->omega-body M-of-q)
  (abstract-to-state-function 
    (M-of-q->omega-body-of-t M-of-q)))

#|
(show-expression
 ((M->omega-body Euler->M)
  (->state 't 
	   (vector 'theta 'phi 'psi)
	   (vector 'thetadot 'phidot 'psidot))))
(vector (+ (* phidot (sin psi) (sin theta)) (* thetadot (cos psi)))
        (+ (* phidot (cos psi) (sin theta)) (* -1 thetadot (sin psi)))
        (+ (* phidot (cos theta)) psidot))


|#

#|

(define ((T-rigid-body A B C) state)
  (let ((omega-body ((M->omega-body Euler->M) state))
	(I (m:make-diagonal (vector A B C))))
    (* 1/2 omega-body (* I omega-body))))

|#


#|

(print-expression
 ((M-of-state->omega-body
   (literal-function 'q->M *exactly-one* *matrix* *vector*))
  (->state 't 
	   (vector 'theta 'phi 'psi)
	   (vector 'thetadot 'phidot 'psidot))))


(print-expression
 ((M-of-state->omega-body
   (matrix-by-rows
    (list (literal-function 'q->M11) (literal-function 'q->M12))
    (list (literal-function 'q->M21) (literal-function 'q->M22))))
  (->state 't 
	   'q
	   'qdot)))

|#