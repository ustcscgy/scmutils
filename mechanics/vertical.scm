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


;;; map out region of stability of vertical equilibrium for driven pendulum

;;; from linear-stability.scm
(define (fixed-point-eigen T xe ye eps cont)
  (let ((M00 ((richardson-derivative 
	       (lambda (dx)
		 (T (+ xe dx) ye 
		    (lambda (x y) ((principal-value pi) (- x xe)))
		    'failure))
	       eps) 0.0))
	(M01 ((richardson-derivative 
	       (lambda (dx)
		 (T xe (+ ye dx) 
		    (lambda (x y) ((principal-value pi) (- x xe)))
		    'failure))
	       eps) 0.0))
	(M10 ((richardson-derivative 
	       (lambda (dx)
		 (T (+ xe dx) ye (lambda (x y) y) 'failure))
	       eps) 0.0))
	(M11 ((richardson-derivative 
	       (lambda (dx)
		 (T xe (+ ye dx) (lambda (x y) y) 'failure))
	       eps) 0.0)))
    (let ((trace (+ M00 M11))
	  (determinant (- (* M00 M11) (* M01 M10))))
      (quadratic 
       1. (- trace) determinant 
       (lambda (root1 root2)
	 (cont root1 M01 (- root1 M00)
	       root2 M01 (- root2 M00)))))))

;;; (dp-map ml^2 mlg mlAw^2 omega)

(define vertical-stable?
  (let ((ml^2 1.)
	(mlg 9.8))
    (lambda (A/l w/w0)
      (let ((w0 (sqrt (/ mlg ml^2)))
	    (mlAw^2 (* ml^2 A/l (square w/w0) (/ mlg ml^2) )))
	(fixed-point-eigen (dp-map ml^2 mlg mlAw^2 (* w0 w/w0)) 
			   pi 0.0
			   1.e-2
			   (lambda (root1 x1 y1 root2 x2 y2)
			     (not (real? root1))))))))
			   

#|
(graphics-close win)
(graphics-clear win)
(define win (frame 0 3 -3 0 600 600))
(plot-line win 0 (log10 (sqrt 2)) 3 (- (log10 (sqrt 2)) 3))

(let loopw ((log-w/w0 2.9))
  (if (> log-w/w0 0.)
      (let loopA ((log-A/l -2.9))
	(if (< log-A/l 0)
	    (begin
	      (if (vertical-stable? (expt 10. log-A/l) (expt 10. log-w/w0))
		  (plot-point win log-w/w0 log-A/l))
	      (loopA (+ log-A/l 0.1)))
	    (loopw (- log-w/w0 0.1))))))

(graphics-clear win)

(let loopw ((log-w/w0 2.9))
  (if (> log-w/w0 0.)
      (let loopA ((log-A/l -2.9))
	(if (< log-A/l 0)
	    (begin
	      (if (vertical-stable? (expt 10. log-A/l) (expt 10. log-w/w0))
		  (graphics-operation win 'fill-circle log-w/w0 log-A/l .01))
	      (loopA (+ log-A/l 0.1)))
	    (loopw (- log-w/w0 0.1))))))
|#


#|

(vertical-stable? (expt 10. -0.2) 10)
;Value: ()
(vertical-stable? (expt 10. -0.2) 100)
;Value: ()
(expt 10. -0.2)
;Value: .6309573444801932

|#


#|

the perturbation theory 

H = H0 + H1 + H2

{H0, W} + H1 = 0

H' = H0 + H2 + 1/2 {H1, W} + {H2, W}

|#

(define ((H0 I alpha beta omega) state)
  (let ((p (state->p state)))
    (let ((ptime (ref p 0))
	  (ptheta (ref p 1)))
      (+ ptime (/ (square ptheta) (* 2 I))))))
      
(define ((H1 I alpha beta omega) state)
  (let ((q (state->q state)))
    (let ((time (ref q 0))
	  (theta (ref q 1)))
      (* beta (+ (cos (- theta (* omega time)))
		 (cos (+ theta (* omega time))))))))

(define ((H2 I alpha beta omega) state)
  (let ((q (state->q state)))
    (let ((time (ref q 0))
	  (theta (ref q 1)))
      (* -1 alpha (cos theta)))))

(define ((W1 I alpha beta omega) state)
  (let ((q (state->q state))
	(p (state->p state)))
    (let ((time (ref q 0))
	  (theta (ref q 1))
	  (ptime (ref p 0))
	  (ptheta (ref p 1)))
      (let ((omega0 (/ ptheta I)))
	(* -1 
	   (+ (/ (* alpha (sin theta)) omega0)
	      (* -1 beta 
		 (+ (/ (sin (- theta (* omega time)))
		       (- omega0 omega))
		    (/ (sin (+ theta (* omega time)))
		       (+ omega0 omega))))))))))

#|

(print-expression
 ((+ (Poisson-bracket (H0 'I 'alpha 'beta 'omega)
		      (W1 'I 'alpha 'beta 'omega))
     (+ (H1 'I 'alpha 'beta 'omega)
	(H2 'I 'alpha 'beta 'omega)))
  (->H-state 'tau
	     (vector 'time 'theta)
	     (vector 'ptime 'ptheta))))
0
;Unspecified return value

|#

;;; exclude the zero frequency term

(define ((W2 I alpha beta omega) state)
  (let ((q (state->q state))
	(p (state->p state)))
    (let ((time (ref q 0))
	  (theta (ref q 1))
	  (ptime (ref p 0))
	  (ptheta (ref p 1)))
      (let ((omega0 (/ ptheta I)))
	(* -1 
	   (* -1 beta 
	      (+ (/ (sin (- theta (* omega time)))
		    (- omega0 omega))
		 (/ (sin (+ theta (* omega time)))
		    (+ omega0 omega)))))))))

#|

(print-expression
 ((+ (Poisson-bracket (H0 'I 'alpha 'beta 'omega)
		      (W2 'I 'alpha 'beta 'omega))
     (H1 'I 'alpha 'beta 'omega))
  (->H-state 'tau
	     (vector 'time 'theta)
	     (vector 'ptime 'ptheta))))
0
;Unspecified return value

(define (L_W2 I alpha beta omega)
  (Lie-derivative (W2 I alpha beta omega)))

(print-expression
 ((+ ((L_W2 'I 'alpha 'beta 'omega)
      (H0 'I 'alpha 'beta 'omega))
     (H1 'I 'alpha 'beta 'omega))
  (->H-state 'tau
	     (vector 'time 'theta)
	     (vector 'ptime 'ptheta))))
0
;Unspecified return value

;;; check L_W^2 H0 = - L_W H1

(print-expression
 ((+ ((square (L_W2 'I 'alpha 'beta 'omega))
      (H0 'I 'alpha 'beta 'omega))
     ((L_W2 'I 'alpha 'beta 'omega)
      (H1 'I 'alpha 'beta 'omega)))
  (->H-state 'tau
	     (vector 'time 'theta)
	     (vector 'ptime 'ptheta))))
0
;Unspecified return value

|#

#|

now compute second order terms

(print-expression
 ((+ (* 1/2 ((L_W2 'I 'alpha 'beta 'omega)
	     (H1 'I 'alpha 'beta 'omega)))
     ((L_W2 'I 'alpha 'beta 'omega)
      (H2 'I 'alpha 'beta 'omega)))
  (->H-state 'tau
	     (vector 'time 'theta)
	     (vector 'ptime 'ptheta))))

a big mess -- need averaging program

|#

#|

just draw contours of H

(define ((Hpp I alpha beta omega) state)
  (let ((q (state->q state))
	(p (state->p state)))
    (+ (* 1/2 (/ (square p) I))
       (* -1 alpha (cos q))
       (* (/ (square beta) (* 4 I))
	  (- 1 (cos (* 2 q)))
	  (+ (/ 1 (square (- (/ p I) omega)))
	     (/ 1 (square (+ (/ p I) omega))))))))

(define (der I alpha beta omega)
  (phase-space-derivative (Hpp I alpha beta omega)))

(define compiled-der 
  (compile-sysder 1 der))

(define ((Hmap I alpha beta omega dt) x y continue fail)
  (let ((state (->H-state 0. x y)))
    (let ((nstate (ode-advancer (compiled-der I alpha beta omega) state dt 1.e-12)))
      (continue ((principal-value pi)
		 (state->q nstate))
		(state->p nstate)))))


alpha = mlg 
beta = mlAw^2/2 
A/l = .1   w = 31.304951684997057
alpha = 9.8
beta = (* .1 (square 31.304951684997057) 1/2) = 49

(define win2 (frame -pi pi -10 10 500 500))
(graphics-clear win2)
(graphics-close win2)
(let ((m 1)
      (l 1)
      (g 9.8)
      (A 0.03)
      (omega (* 100 (sqrt 9.8))))
  (let ((I (* m l l))
	(alpha (* m l g))
	(beta (* 1/2 m l A (square omega))))
    (explore-map win2 (Hmap I alpha beta omega .001))))

(define win2 (frame -pi pi -100 100 500 500))
(graphics-clear win2)
(graphics-close win2)
(let ((m 1)
      (l 1)
      (g 9.8)
      (A 0.2)
      (omega (* 100 (sqrt 9.8))))
  (let ((I (* m l l))
	(alpha (* m l g))
	(beta (* 1/2 m l A (square omega))))
    (explore-map win2 (Hmap I alpha beta omega .001))))

|#

