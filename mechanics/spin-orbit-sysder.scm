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


;;; spin-orbit

;;; C, n = 1

(define ((spinorbit-sysder eps e) s)
  (let ((t (ref s 0))
	(theta (ref s 1))
	(thetadot (ref s 2))
	(f (ref s 3)))
    (let ((a/R (/ (+ 1 (* e (cos f))) (- 1 (square e)))))
      (vector 1.
	      thetadot
	      (* -1/2 (square eps) (cube a/R) 
		 (sin (* 2 (- theta f))))
	      (* (sqrt (- 1 (square e))) (square a/R))))))

(define spinorbit-sysder-compiled
  (compile-parametric 4 spinorbit-sysder))

(define (spinorbit-map eps e)
  (let ((n 5) (pv (principal-value pi)))
    (let ((dt (/ 2pi n)))
      (lambda (theta thetadot cont fail)
	(let loop ((s (->local 0. theta thetadot 0.)) (n n))
	  (if (fix:= n 0)
	      (let ((theta (pv (ref s 1)))
		    (thetadot (velocity s)))
		(cont theta thetadot))
	      (loop (ode-advancer (spinorbit-sysder-compiled eps e)
				  s
				  dt
				  1.e-13) 
		    (fix:- n 1))))))))

#|

(define win (frame (* -pi 5/4) (* pi 5/4) -0.25 2.25 500 500))

(define (do-it eps e)
  (plot-line win -pi 0. pi 0.)
  (plot-line win -pi 2. pi 2.)
  (plot-line win pi 0 pi 2.)
  (plot-line win -pi 0 -pi 2.)

  (plot-line win -3.5 (- 1. eps) -3.5 (+ 1. eps))
  (plot-line win -3.5 1. -3.4 1.)
  (plot-line win -3.5 (+ 1. eps) -3.4 (+ 1. eps))
  (plot-line win -3.5 (- 1. eps) -3.4 (- 1. eps))

  (plot-line win -3.5 (- 1.5 (* eps (sqrt (* 7/2 e)))) -3.5 (+ 1.5 (* eps (sqrt (* 7/2 e)))))
  (plot-line win -3.5 1.5 -3.4 1.5)
  (plot-line win -3.5 (+ 1.5 (* eps (sqrt (* 7/2 e)))) -3.4 (+ 1.5 (* eps (sqrt (* 7/2 e)))))
  (plot-line win -3.5 (- 1.5 (* eps (sqrt (* 7/2 e)))) -3.4 (- 1.5 (* eps (sqrt (* 7/2 e)))))

  (plot-line win -3.5 (- 0.5 (* eps (sqrt (* 1/2 e)))) -3.5 (+ 0.5 (* eps (sqrt (* 1/2 e)))))
  (plot-line win -3.5 0.5 -3.4 0.5)
  (plot-line win -3.5 (+ 0.5 (* eps (sqrt (* 1/2 e)))) -3.4 (+ 0.5 (* eps (sqrt (* 1/2 e)))))
  (plot-line win -3.5 (- 0.5 (* eps (sqrt (* 1/2 e)))) -3.4 (- 0.5 (* eps (sqrt (* 1/2 e)))))

  (explore-map win (spinorbit-map eps e) 1000)
  )

(do-it 0.15 0.1) ; show resonance widths
(do-it 0.27 0.1) ; show connected chaos
(do-it 0.25 0.1) ; 

predicted overlap 
(/ 1 (+ 2 (sqrt (* 14 .1))))
;Value: .3141477089923372




|#
