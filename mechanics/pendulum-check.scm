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

(define ((check-action alpha beta) state)
  (let ((E ((Hpendulum alpha beta) state)))
    (write-line (if (< E beta) 'librating 'circulating))
    (if (< E beta)
	;; oscillating
	(let ((theta-max (acos (/ (- E) beta))))
	  (* 4  (/ 1 2pi)
	     (definite-integral 
	       (lambda (theta) (sqrt (* 2 alpha (+ E (* beta (cos theta))))))
	       0.
	       theta-max
	       1.e-13)))
	;; circulating
	(* 2 (/ 1 2pi)
	   (definite-integral 
	     (lambda (theta) (sqrt (* 2 alpha (+ E (* beta (cos theta))))))
	     0.
	     pi
	     1.e-13)))))

|#


#|

((check-action 2. 9.8) (->H-state 0 4. 4.))
circulating
;Value: 5.951793337936512

((pendulum-action 2. 9.8) (->H-state 0 4. 4.))
circulating
;Value: 5.951793337937046

((check-action 2. 9.8) (->H-state 0 1. 3.))
librating
;Value: 3.203964986352005

((pendulum-action 2. 9.8) (->H-state 0 1. 3.))
librating
;Value: 3.203964986292995

|#

#|

; thetadot = p/alpha
(define ((check-frequency-oscillation alpha beta) theta0)
  (let ((theta0 (abs theta0)))
    (let loop ((state (->H-state 0. theta0 0.)))
      (if (< (state->q state) 0.0)
	  (let refine ((state state))
	    (if (< (abs ((principal-value pi) (state->q state))) 1.e-13)
		(/ 2pi (* 4 (state->t state)))
		(refine (ode-advancer 
			 (pendulum-sysder-compiled alpha beta)
			 state
			 (- (/ (state->q state) (/ (state->p state) alpha)))
			 1.e-13))))
	  (loop (ode-advancer (pendulum-sysder-compiled alpha beta) state .01 1.e-13))))))

((check-frequency-oscillation 2. 9.8) 1.)
;Value: 2.0758916552604294
(pendulum-frequency 2. 9.8 ((Hpendulum 2. 9.8) (->H-state 0. 1. 0.)))
;Value: 2.0758916552604245

(sqrt (/ 9.8 2.))
;Value: 2.2135943621178655
(pendulum-frequency 2. 9.8 ((Hpendulum 2. 9.8) (->H-state 0. 0. 0.)))
;Value: 2.2135943621178655

(define ((check-frequency-circulation alpha beta) ptheta0)
  (let ((ptheta0 (abs ptheta0)))
    (let loop ((state (->H-state 0. -pi ptheta0)))
      (if (> (state->q state) 0)
	  (let refine ((state state))
	    (if (< (abs ((principal-value pi) (state->q state))) 1.e-13)
		(/ 2pi (* 2 (state->t state)))
		(refine (ode-advancer 
			 (pendulum-sysder-compiled alpha beta)
			 state
			 (- (/ (state->q state) (/ (state->p state) alpha)))
			 1.e-13))))
	  (loop (ode-advancer (pendulum-sysder-compiled alpha beta) state .01 1.e-13))))))

((check-frequency-circulation 2. 9.8) 1.)
;Value: 1.9539382280411355

(pendulum-frequency 2. 9.8 ((Hpendulum 2. 9.8) (->H-state 0. -pi 1.)))
;Value: 1.9539382280411273

|#

#|

(let ((alpha 2.) (beta 9.8))
  (let ((state0 (->state 0. 0. 2.)))
    (let ((E ((Hpendulum alpha beta) state0))
	  (time 2.1))
      (let ((theta ((pendulum-oscillating-angle alpha beta E) time))
	    (ptheta ((pendulum-oscillating-angular-momentum alpha beta E) time)))
	(write-line (list time theta ptheta))
	(write-line (((check-solution alpha beta 1.e-13) state0) time))))))

|#


#|

(let ((alpha 2.) (beta 9.8))
  (let ((state0 (->state 0. 0. 15.)))
    (let ((E ((Hpendulum alpha beta) state0))
	  (time 1.1))
      (let ((theta ((pendulum-circulating-angle alpha beta E) time))
	    (ptheta ((pendulum-circulating-angular-momentum alpha beta E) time)))
	(write-line (list time theta ptheta))
	(write-line (((check-solution alpha beta 1.e-13) state0) time))))))

|#


#|

(let ((alpha 2.) (beta 9.8) (time 20.))
  ((Hpendulum alpha beta) 
   (->H-state time 
	      ((pendulum-separatrix-angle alpha beta) time)
	      ((pendulum-separatrix-angular-momentum alpha beta) time))))
|#
