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

; a particle is a pair of a mass and a function f: (t, q)->x
; ((V masses) time positions) -> potential-energy

(define ((make-Lagrangian V particles) state)
  (pp particles)
  (let ((t (state->t state))
	(q (state->q state))
	(qd (state->qdot state)))
    (let ((path (osculating-path t q qd)))
      (let ((x-paths 
	     (map (lambda (particle)
		    (let ((f (cdr particle)))
		      (lambda (t) (f t (path t)))))
		  particles))
	    (ms (map car particles)))
	(let ((xs (map (lambda (x) (x t)) x-paths))
	      (vs (map (lambda (x) ((derivative x) t)) x-paths)))
	  (let ((kinetic-energy 
		 (apply + 
			(map (lambda (m v) (* 1/2 m (square v)))
			     ms vs)))
		(potential-energy ((V ms) t xs)))
	    (- kinetic-energy
	       potential-energy)))))))

#|

;;;----------------------------------------------------------------

;;; useful in general

(define (height x) (vector-ref x 2))

(define (((V-uniform-gravity g) masses) time positions)
  (apply + (map (lambda (m x) (* m g (height x))) masses positions)))

;;;----------------------------------------------------------------
;;; driven pendulum

(define (L-dpend-auto m g l ys)
  (make-Lagrangian 
   (V-uniform-gravity g)
   (list 
    (cons m  
	  (lambda (t q) 
	    (vector (* l (sin q))
		    0
		    (- (ys t) (* l (cos q)))))))))
			

(print-expression 
 ((L-dpend-auto 'm 'g 'l (literal-function 'ys))
  (->state 't0 'theta 'thetadot)))

(+ (* g l m (cos theta))
   (* -1 g m (ys t0))
   (* 1/2 (expt l 2) m (expt thetadot 2))
   (* l m thetadot ((D ys) t0) (sin theta))
   (* 1/2 m (expt ((D ys) t0) 2)))
;Unspecified return value

;;; CORRECT

(print-expression 
 (((lagrange-equations
    (L-dpend-auto 'm 'g 'l (literal-function 'ys)))
   (literal-function 'theta))
  't))

(+ (* g l m (sin (theta t)))
   (* (expt l 2) m (((expt D 2) theta) t))
   (* l m (((expt D 2) ys) t) (sin (theta t))))
;Unspecified return value

;;; ALSO CORRECT

;;;----------------------------------------------------------------
;;; double pendulum

(define (L-double-pend m1 l1 m2 l2 g)
  (let ((p1 (lambda (t q)
	      (let ((t1 (vector-ref q 0)))
		(vector (* l1 (sin t1))
			0
			(- (* l1 (cos t1))))))))
    (let ((p2 (lambda (t q)
		(let ((t2 (vector-ref q 1)))
		  (+ (p1 t q)
		     (vector (* l2 (sin t2))
			     0
			     (- (* l2 (cos t2)))))))))
      (pp (list p1 p2))
      (make-lagrangian
       (V-uniform-gravity g)
       (list (cons m1 p1)
	     (cons m2 p2))))))


|#