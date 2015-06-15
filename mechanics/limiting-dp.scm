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

;;; find driven pendulum by taking small mass limit

; first with the pendulum constraint built in but with a dynamical pivot
(define ((L mp m l g k) local)
  (let ((q (coordinate local))
	(v (velocity local)))
    (let ((xp (ref q 0))
	  (yp (ref q 1))
	  (theta (ref q 2))
	  (vxp (ref v 0))
	  (vyp (ref v 1))
	  (thetadot (ref v 2)))
      (let ((x (+ xp (* l (sin theta))))
	    (y (- yp (* l (cos theta))))
	    (vx (+ vxp (* l thetadot (cos theta))))
	    (vy (+ vyp (* l thetadot (sin theta)))))
	(- (+ (* 1/2 (* mp (+ (square vxp) (square vyp))))
	      (* 1/2 (* m (+ (square (+ vx vxp)) (square (+ vy vyp))))))
	   (+ (* m g (+ y yp))
	      (* 1/2 k (+ (square xp) (square yp)))))))))
      
(show-expression
 (((Lagrange-equations (L 'm_p 'm 'l 'g 'k))
   (coordinate-tuple (literal-function 'x_p)
		     (literal-function 'y_p)
		     (literal-function 'theta)))
  't))
 

(show-expression
 ((LE (L 'm_p 'm 'l 'g 'k))
  (->local 't 
	   (coordinate-tuple 'x 'y 'theta)
	   (velocity-tuple 'v_x 'v_y 'thetadot)
	   (acceleration-tuple 'a_x 'a_y 'thetadotdot))))
(row
 (+ (* -2 l m (expt thetadot 2) (sin theta))
    (* 2 l m thetadotdot (cos theta))
    (* 4 a_x m)
    (* a_x m_p)
    (* k x))
 (+ (* 2 l m (expt thetadot 2) (cos theta))
    (* 2 l m thetadotdot (sin theta))
    (* 4 a_y m)
    (* a_y m_p)
    (* 2 g m)
    (* k y))
 (+ (* 2 a_x l m (cos theta))
    (* 2 a_y l m (sin theta))
    (* g l m (sin theta))
    (* (expt l 2) m thetadotdot)))

; ok in the limit m_p >> m we get the usual lagrange equations

(define ((L mp m l g k) local)
  (let ((q (coordinate local))
	(v (velocity local)))
    (let ((xp (ref q 0))
	  (yp (ref q 1))
	  (x (ref q 2))
	  (y (ref q 3))
	  (lambda (ref q 4))
	  (vxp (ref v 0))
	  (vyp (ref v 1))
	  (vx (ref v 2))
	  (vy (ref v 3)))
      (+ (- (+ (* 1/2 (* mp (+ (square vxp) (square vyp))))
	       (* 1/2 (* m (+ (square (+ vx vxp)) (square (+ vy vyp))))))
	    (+ (* m g (+ y yp))
	       (* 1/2 k (+ (square xp) (square yp)))))
	 (* lambda 
	    (- (+ (square x) (square y)) (square l)))))))


(show-expression
 ((LE (L 'm_p 'm 'l 'g 'k))
  (->local 't 
	   (coordinate-tuple 'xi 'eta 'x 'y 'lambda)
	   (velocity-tuple 'v_xi 'v_eta 'v_x 'v_y 'lambdadot)
	   (acceleration-tuple 'a_xi 'a_eta 'a_x 'a_y 'lambdadotdot))))