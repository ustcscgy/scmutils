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

;;;; Driven pendulum example

(define ((T-pend m l g ys) state)
  (let ((t (state->t state))
	(theta (state->q state))
	(thetadot (state->qdot state)))
    (let ((ysdot (D ys)))
      (* 1/2 m
	 (+ (* (square (* l thetadot)))
	    (* (square (ysdot t)))
	    (* 2 (ysdot t) l (sin theta) thetadot))))))

(define ((V-pend m l g ys) state)
  (let ((t (state->t state))
	(theta (state->q state))
	(thetadot (state->qdot state)))
    (* m g (- (ys t) (* l (cos theta))))))

(define L-pend (- T-pend V-pend))


(define (periodically-driven-pendulum m l g a omega)
  (L-pend m l g (lambda (t) (* a (cos (* omega t))))))

(define (H-pend-sysder m l g a omega)
  (phase-space-derivative
   (Lagrangian->Hamiltonian
    (periodically-driven-pendulum m l g a omega))))

(define H-pend-compiled-sysder
  (compile-sysder 1 H-pend-sysder))

(define (driven-pendulum-map mass l g a omega)
  (let ((map-period (/ 2pi omega)))
    (lambda (theta ptheta continue fail)
      (let ((ns
	     (ode-advancer
	      (H-pend-compiled-sysder mass l g a omega)
	      (->state			;state to start from
	       0.0			;t0 = phase/drive-freq.
	       theta
	       ptheta) 
	      map-period			
	      1.0e-10)))
	(continue 
	 ((principal-value pi) (vector-ref ns 1))
	 (vector-ref ns 2))))))
