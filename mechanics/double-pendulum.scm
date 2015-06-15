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

(define ((Lf2 m1 m2 g) local)
  (let ((q (coordinate local))
        (v (velocity local)))
    (let ((x1 (ref q 0))
	  (y1 (ref q 1))
	  (x2 (ref q 2))
	  (y2 (ref q 3))
	  (vx1 (ref v 0))
	  (vy1 (ref v 1))
	  (vx2 (ref v 2))
	  (vy2 (ref v 3)))
      (- (+ (* 1/2 m1 (+ (square vx1) (square vy1)))
	    (* 1/2 m2 (+ (square vx2) (square vy2))))
	 (+ (* m1 g y1) (* m2 g y2))))))


(define ((double-pend-coordinates l1 l2) t q)
  (let ((theta1 (ref q 0))
	(theta2 (ref q 1)))
    (let ((x1 (* l1 (sin theta1)))
	  (y1 (- l1 (* l1 (cos theta1)))))
      (let ((x2 (+ x1 (* l2 (sin theta2))))
	    (y2 (- y1 (* l2 (cos theta2)))))
	(coordinate-tuple x1 y1 x2 y2)))))


(define (L-double-pend m1 m2 l1 l2 g)
  (compose (Lf2 m1 m2 g)
           (F->C (double-pend-coordinates l1 l2))))


#|


(show-expression
 ((LE (L-double-pend 'm1 'm2 'l1 'l2 'g))
  (->local 't 
	   (coordinate-tuple 'theta_1 'theta_2)
	   (velocity-tuple 'thetadot_1 'thetadot_2)
	   (acceleration-tuple 'thetadotdot_1 'thetadotdot_2))))

(row
 (+ (* l1 l2 m2 (expt thetadot_2 2) (cos theta_2) (sin theta_1))
    (* -1 l1 l2 m2 (expt thetadot_2 2) (sin theta_2) (cos theta_1))
    (* l1 l2 m2 thetadotdot_2 (cos theta_2) (cos theta_1))
    (* l1 l2 m2 thetadotdot_2 (sin theta_2) (sin theta_1))
    (* g l1 m1 (sin theta_1))
    (* g l1 m2 (sin theta_1))
    (* (expt l1 2) m1 thetadotdot_1)
    (* (expt l1 2) m2 thetadotdot_1))
 (+ (* -1 l1 l2 m2 (expt thetadot_1 2) (cos theta_2) (sin theta_1))
    (* l1 l2 m2 (expt thetadot_1 2) (sin theta_2) (cos theta_1))
    (* l1 l2 m2 thetadotdot_1 (cos theta_2) (cos theta_1))
    (* l1 l2 m2 thetadotdot_1 (sin theta_2) (sin theta_1))
    (* g l2 m2 (sin theta_2))
    (* (expt l2 2) m2 thetadotdot_2)))


(print-expression
 ((Lagrangian->acceleration
   (L-double-pend 'm1 'm2 'l1 'l2 'g))
  (->local 't 
	   (coordinate-tuple 'theta_1 'theta_2)
	   (velocity-tuple 'thetadot_1 'thetadot_2))))

|#

(define ((double-pend-coordinates-x-theta l1 l2) t q)
  (let ((x1 (ref q 0))
	(theta (ref q 1)))
    (let ((y1 (- l1 (sqrt (- (square l1) (square x1))))))
      (let ((x2 (+ x1 (* l2 (sin theta))))
	    (y2 (- y1 (* l2 (cos theta)))))
	(coordinate-tuple x1 y1 x2 y2)))))


(define (L-double-pend-x-theta m1 m2 l1 l2 g)
  (compose (Lf2 m1 m2 g)
           (F->C (double-pend-coordinates-x-theta l1 l2))))


#|

(print-expression
 ((Lagrangian->acceleration
   (L-double-pend-x-theta 'm1 'm2 'l1 'l2 'g))
  (->local 't 
	   (coordinate-tuple 'x 'theta)
	   (velocity-tuple 'xdot 'thetadot))))

(column
 (/
  (+ (* (expt l1 4) l2 m2 (expt thetadot 2) (sin theta))
     (* -2 (expt l1 2) l2 m2 (expt thetadot 2) (expt x 2) (sin theta))
     (* l2 m2 (expt thetadot 2) (expt x 4) (sin theta))
     (* -1 l2 m2 (expt thetadot 2) x (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3) (cos theta))
     (* g (expt l1 4) m2 (cos theta) (sin theta))
     (* -2 g (expt l1 2) m2 (expt x 2) (cos theta) (sin theta))
     (* g m2 (expt x 4) (cos theta) (sin theta))
     (* -1 g m2 x (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3) (expt (cos theta) 2))
     (* -1 (expt l1 2) m2 x (expt xdot 2) (expt (cos theta) 2))
     (* m2 (expt x 2) (expt xdot 2) (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) (cos theta) (sin theta))
     (* m2 (expt xdot 2) (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3) (cos theta) (sin theta))
     (* -1 g m1 x (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3))
     (* -1 (expt l1 2) m1 x (expt xdot 2)))
  (+ (* (expt l1 4) m2 (expt (sin theta) 2))
     (* 3 (expt l1 2) m2 (expt x 2) (expt (cos theta) 2))
     (* -2 m2 (expt x 4) (expt (cos theta) 2))
     (* -2 m2 x (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3) (cos theta) (sin theta))
     (* (expt l1 4) m1)
     (* -1 (expt l1 2) m1 (expt x 2))
     (* -2 (expt l1 2) m2 (expt x 2))
     (* m2 (expt x 4))))
 (/
  (+ (* -1 (expt l1 4) l2 m2 (expt thetadot 2) (cos theta) (sin theta))
     (* 3 (expt l1 2) l2 m2 (expt thetadot 2) (expt x 2) (cos theta) (sin theta))
     (* -2 l2 m2 (expt thetadot 2) (expt x 4) (cos theta) (sin theta))
     (* 2 l2 m2 (expt thetadot 2) x (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3) (expt (cos theta) 2))
     (* -1 l2 m2 (expt thetadot 2) x (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3))
     (* -1 g (expt l1 4) m1 (sin theta))
     (* -1 g (expt l1 4) m2 (sin theta))
     (* 2 g (expt l1 2) m1 (expt x 2) (sin theta))
     (* 2 g (expt l1 2) m2 (expt x 2) (sin theta))
     (* -1 g m1 (expt x 4) (sin theta))
     (* g m1 x (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3) (cos theta))
     (* -1 g m2 (expt x 4) (sin theta))
     (* g m2 x (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3) (cos theta))
     (* (expt l1 2) m1 x (expt xdot 2) (cos theta))
     (* (expt l1 2) m2 x (expt xdot 2) (cos theta))
     (* -1 m1 (expt x 2) (expt xdot 2) (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) (sin theta))
     (* -1 m1 (expt xdot 2) (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3) (sin theta))
     (* -1 m2 (expt x 2) (expt xdot 2) (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) (sin theta))
     (* -1 m2 (expt xdot 2) (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3) (sin theta)))
  (+ (* (expt l1 4) l2 m2 (expt (sin theta) 2))
     (* 3 (expt l1 2) l2 m2 (expt x 2) (expt (cos theta) 2))
     (* -2 l2 m2 (expt x 4) (expt (cos theta) 2))
     (* -2 l2 m2 x (expt (sqrt (+ (expt l1 2) (* -1 (expt x 2)))) 3) (cos theta) (sin theta))
     (* (expt l1 4) l2 m1)
     (* -1 (expt l1 2) l2 m1 (expt x 2))
     (* -2 (expt l1 2) l2 m2 (expt x 2))
     (* l2 m2 (expt x 4)))))

limit for large l1 is

(show-expression 
 '(column
   (/
    (+ (*  l2 m2 (expt thetadot 2) (sin theta))
       (* g  m2 (cos theta) (sin theta)))
    (+ (*  m2 (expt (sin theta) 2))
       (*  m1)))
   (/
    (+ (* -1  l2 m2 (expt thetadot 2) (cos theta) (sin theta))
       (* -1 g  m1 (sin theta))
       (* -1 g  m2 (sin theta)))
    (+ (*  l2 m2 (expt (sin theta) 2))
       (*  l2 m1)))))

|#


#|  for sliding pendulum

;;; from systems.scm
(define ((L-sliding-pend m1 m2 b g) state)
  (let ((q (coordinate state))
	(qdot (velocity state)))
    (let* ((x (ref q 0))
	   (xdot (ref qdot 0))
	   (theta (ref q 1))
	   (thetadot (ref qdot 1))
	   (rel-pend-vel
	    (* b thetadot (vector (cos theta) (sin theta))))
	   (pend-vel (+ rel-pend-vel (vector xdot 0)))
	   (Tpend (* 1/2 m2 (square pend-vel)))
	   (Tsupport (* 1/2 m1 (square xdot)))
	   (V (- (* m2 g b (cos theta)))))
      (+ Tpend Tsupport (- V)))))

(show-expression
 ((Lagrangian->acceleration
   (L-sliding-pend 'm1 'm2 'l2 'g))
  (->local 't 
	   (coordinate-tuple 'x 'theta)
	   (velocity-tuple 'xdot 'thetadot))))

gives same answer

|#