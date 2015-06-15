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

#|
;;; sliding pendulum 

;;; from lag, below is a better one...
(define ((L-sliding-pend-0 m1 m2 b g) state)
  (let ((q (coordinate state))
	(qdot (velocity state)))
    (let ((x (ref q 0))
	  (theta (ref q 1))
	  (xdot (ref qdot 0))
	  (thetadot (ref qdot 1)))
      (let* ((rel-pend-vel
	      (* b thetadot (velocity-tuple (cos theta) (sin theta))))
	     (pend-vel (+ rel-pend-vel (velocity-tuple xdot 0)))
	     (Tpend (* 1/2 m2 (square pend-vel)))
	     (Tsupport (* 1/2 m1 (square xdot)))
	     (V (- (* m2 g b (cos theta)))))
	(+ Tpend Tsupport (- V))))))

(show-expression
 (((Lagrange-equations (L-sliding-pend-0 'm_1 'm_2 'b 'g))
   (coordinate-tuple (literal-function 'x)
		     (literal-function 'theta)))
  't))
(down
 (+ (* -1 b m_2 (sin (theta t)) (expt ((D theta) t) 2))
    (* b m_2 (((expt D 2) theta) t) (cos (theta t)))
    (* m_1 (((expt D 2) x) t))
    (* m_2 (((expt D 2) x) t)))
 (+ (* (expt b 2) m_2 (((expt D 2) theta) t))
    (* b g m_2 (sin (theta t)))
    (* b m_2 (((expt D 2) x) t) (cos (theta t)))))
|#

(define ((F-sliding-pend l) t q)
  (let ((x (ref q 0))
	(theta (ref q 1)))
    (let ((pos1 (coordinate-tuple x 0))
	  (pos2 (coordinate-tuple 
		 (+ x (* l (sin theta))) 
		 (- 0 (* l (cos theta))))))
      (coordinate-tuple pos1 pos2))))

(define ((L-uniform-g m1 m2 g) local)
  (let ((q (coordinate local))
	(v (velocity local)))
    (let ((pos1 (ref q 0))
	  (pos2 (ref q 1))
	  (vel1 (ref v 0))
	  (vel2 (ref v 1)))
      (let ((T (+ (* 1/2 m1 (square vel1))
		  (* 1/2 m2 (square vel2))))
	    (V (+ (* m1 g (ref pos1 1))
		  (* m2 g (ref pos2 1)))))
	(- T V)))))

(define (L-sliding-pend m1 m2 l g)
  (compose (L-uniform-g m1 m2 g)
	   (F->C (F-sliding-pend l))))
  
#|

(pe
 ((- (L-sliding-pend 'm_1 'm_2 'l 'g) 
     (L-sliding-pend-0 'm_1 'm_2 'l 'g))
  (->local 't
	   (coordinate-tuple 'x 'theta)
	   (velocity-tuple 'xdot 'thetadot))))
0

(show-expression
 (((Lagrange-equations 
    (L-sliding-pend 'm_1 'm_2 'l 'g))
   (coordinate-tuple (literal-function 'x)
		     (literal-function 'theta)))
  't))
"\\boxit{ $$\\left[ \\matrix{ \\displaystyle{  - l m_{2} \\sin\\left(
\\theta\\left( t \\right) \\right) \\left( D\\theta\\left( t \\right)
\\right)^{2} + l m_{2} D^{2}\\theta\\left( t \\right) \\cos\\left(
\\theta\\left( t \\right) \\right) + m_{1} D^{2}x\\left( t \\right) +
m_{2} D^{2}x\\left( t \\right)} \\cr \\cr \\displaystyle{ g l m_{2}
\\sin\\left( \\theta\\left( t \\right) \\right) + l^{2} m_{2}
D^{2}\\theta\\left( t \\right) + l m_{2} D^{2}x\\left( t \\right)
\\cos\\left( \\theta\\left( t \\right) \\right)}} \\right]$$}"


(show-expression
 ((Euler-Lagrange-operator 
   (L-sliding-pend 'm_1 'm_2 'b 'g))
  (->local 't
	   (coordinate-tuple 'x 'theta)
	   (velocity-tuple 'xdot 'thetadot)
	   (acceleration-tuple 'xdotdot 'thetadotdot))))
"\\boxit{ $$\\left[ \\matrix{ \\displaystyle{  - b m_{2}
\\dot{\\theta}^{2} \\sin\\left( \\theta \\right) + b m_{2}
\\ddot{\\theta} \\cos\\left( \\theta \\right) + m_{1} \\ddot{x} +
m_{2} \\ddot{x}} \\cr \\cr \\displaystyle{ b^{2} m_{2} \\ddot{\\theta}
+ b g m_{2} \\sin\\left( \\theta \\right) + b m_{2} \\ddot{x}
\\cos\\left( \\theta \\right)}} \\right]$$}"

(show-expression
 ((Lagrangian->acceleration 
  (L-sliding-pend 'm_1 'm_2 'b 'g))
   (->local 't
	   (coordinate-tuple 'x 'theta)
	   (velocity-tuple 'xdot 'thetadot))))
;Value 203: 
"\\boxit{ $$\\left( \\matrix{ \\displaystyle{ {{b m_{2} \\dot{\\theta}^{2} \\sin\\left( \\theta \\right)}\\over {m_{2} \\left( \\sin\\left( \\theta \\right) \\right)^{2} + m_{1}}} + {{g m_{2} \\sin\\left( \\theta \\right) \\cos\\left( \\theta \\right)}\\over {m_{2} \\left( \\sin\\left( \\theta \\right) \\right)^{2} + m_{1}}}} \\cr \\cr \\displaystyle{  - {{m_{2} \\dot{\\theta}^{2} \\sin\\left( \\theta \\right) \\cos\\left( \\theta \\right)}\\over {m_{2} \\left( \\sin\\left( \\theta \\right) \\right)^{2} + m_{1}}} - {{g m_{1} \\sin\\left( \\theta \\right)}\\over {b m_{2} \\left( \\sin\\left( \\theta \\right) \\right)^{2} + b m_{1}}} - {{g m_{2} \\sin\\left( \\theta \\right)}\\over {b m_{2} \\left( \\sin\\left( \\theta \\right) \\right)^{2} + b m_{1}}}}} \\right)$$}"

(show-expression
 ((Lagrangian->state-derivative
  (L-sliding-pend 'm_1 'm_2 'b 'g))
   (->local 't
	   (coordinate-tuple 'x 'theta)
	   (velocity-tuple 'xdot 'thetadot))))
;Value 204: "\\boxit{ $$\\left( \\matrix{ \\displaystyle{ 1} \\cr \\cr \\displaystyle{ \\dot{x}} \\cr \\cr \\displaystyle{ \\dot{\\theta}} \\cr \\cr \\displaystyle{ {{b m_{2} \\dot{\\theta}^{2} \\sin\\left( \\theta \\right)}\\over {m_{2} \\left( \\sin\\left( \\theta \\right) \\right)^{2} + m_{1}}} + {{g m_{2} \\sin\\left( \\theta \\right) \\cos\\left( \\theta \\right)}\\over {m_{2} \\left( \\sin\\left( \\theta \\right) \\right)^{2} + m_{1}}}} \\cr \\cr \\displaystyle{  - {{m_{2} \\dot{\\theta}^{2} \\sin\\left( \\theta \\right) \\cos\\left( \\theta \\right)}\\over {m_{2} \\left( \\sin\\left( \\theta \\right) \\right)^{2} + m_{1}}} - {{g m_{1} \\sin\\left( \\theta \\right)}\\over {b m_{2} \\left( \\sin\\left( \\theta \\right) \\right)^{2} + b m_{1}}} - {{g m_{2} \\sin\\left( \\theta \\right)}\\over {b m_{2} \\left( \\sin\\left( \\theta \\right) \\right)^{2} + b m_{1}}}}} \\right)$$}"


|#

#|


(define (sliding-pend-sysder m_1 m_2 l g)
  (Lagrangian->state-derivative
   (L-sliding-pend m_1 m_2 l g)))

(define sliding-pend-sysder-compiled
  (compile-parametric 5 sliding-pend-sysder))

(define plot-win (frame 0. 100. -pi pi))

(define (do-it sysder istate0 final-t dt tol)
  (let loop ((istate istate0))
    (if (< (istate->t istate) final-t)
      (let ((ns
             (ode-advancer sysder istate dt tol)))
        (let ((local (istate->local ns)))
          (let ((theta 
                 ((principal-value pi) 
                  (ref (coordinate local) 1))))
            (plot-point plot-win
                        (time local)
                        theta)))
        (loop ns)))))

(graphics-clear plot-win)
(do-it (sliding-pend-sysder-compiled 
	0.1				; m_1 
	1.0				; m_2    
	1.0				; l
	9.8)				; g
       (->istate 0.0         
                 (coordinate-tuple 1. 1.)
		 (velocity-tuple 0. 1.))
       100.0                 
       0.01                  
       1.0e-15)      

;;; alternatively
(define (do-it sysder monitor istate0 final-t dt tol)
  (let loop ((istate istate0))
    (if (< (istate->t istate) final-t)
      (let ((ns (ode-advancer sysder istate dt tol)))
	(monitor (istate->local ns))
        (loop ns)))))

(define ((monitor-theta win) local)
  (let ((theta (ref (coordinate local) 1)))
    (plot-point win
		(time local)
		((principal-value pi) theta))
    local))

(define ((monitor-x win) local)
  (let ((x (ref (coordinate local) 0)))
    (plot-point win
		(time local)
		x)
    local))

(define win (frame 0. 10. -pi pi))

(do-it (sliding-pend-sysder-compiled 
	.1				; m_1 
	10.0				; m_2    
	1.0				; l
	9.8)				; g
       (compose (monitor-theta win)
		(monitor-x win))
       (->istate 0.0         
                 (coordinate-tuple 0. pi/2)
		 (velocity-tuple 0. 0.))
       10.0                 
       0.001                  
       1.0e-15)


(define ((animate-sp win m1 m2 l g dt) istate0)
  (let loop ((istate istate0))
    (let ((local (istate->local istate)))
      (let ((t (time local))
	    (q (coordinate local)))
	(let ((pos ((F-sliding-pend l) t q)))
	  (let ((x1 (ref (ref pos 0) 0))
		(y1 (ref (ref pos 0) 1))
		(x2 (ref (ref pos 1) 0))
		(y2 (ref (ref pos 1) 1)))
	    (graphics-clear win)
	    (plot-line win x1 y1 x2 y2)
					;(graphics-flush win)
	    (loop (ode-advancer (sliding-pend-sysder-compiled m1 m2 l g)
				istate 
				dt
				1.e-13))))))))

(define win (frame -2. 2. -2. 2.))
((animate-sp win .1 1. 1. 9.8 .001)
 (->istate 0. 0. 1. -.1 0.))

|#



