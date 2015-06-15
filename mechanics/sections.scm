#| -*-Scheme-*-

$Id: copyright.scm,v 1.4 2005/12/13 06:41:00 cph Exp $

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

(declare (usual-integrations))

;;; Now to do the Poincare section.  
;;;  Map explorer:
;;;   Left button starts a trajectory.
;;;   Middle button continues a trajectory.
;;;   Right button interrogates coordinates.

(define (explore-map window poincare-map #!optional n)
  (define (iterate-map i x y)
    (if (fix:> i 0) 
	(poincare-map 
	 x y
	 (lambda (nx ny)
	   (plot-point window x y)
	   (iterate-map (fix:- i 1) nx ny))
	 (lambda ()  
	   (newline)
	   (display "Illegal point: ")
	   (write (list  x y))
	   (button-loop x y)))
	(button-loop x y)))
  (define (button-loop ox oy)
    (pointer-coordinates 
     window 
     (lambda (x y button)
       (case button
	 ((0)
	  (write-line (list x y))
	  (display " started.")
	  (iterate-map n x y))
	 ((1) 
	  (write-line (list ox oy))
	  (display " continued.")
	  (iterate-map n ox oy))
	 ((2)
	  (write-line (list x y))
	  (display " hit.")
	  (button-loop ox oy))))))
  (if (default-object? n) (set! n 1000))
  (newline)
  (display "Left button starts a trajectory.")
  (newline)
  (display "Middle button continues a trajectory.")
  (newline)
  (display "Right button interrogates coordinates.")
  (button-loop 9. 9.))

(define (pointer-coordinates window continue)
  (beep)
  (get-pointer-coordinates window continue))
  
;;; For iterating maps of the shape used by explore-map.

(define ((iterated-map map n) x y continue fail)
  (if (fix:< n 0) (error "iterated-map: cannot invert map"))
  (let loop ((x x) (y y) (i n))
    (if (fix:= i 0) 
	(continue x y)
	(map x y
	     (lambda (nx ny)
	       (loop nx ny (fix:- i 1)))
	     fail))))


(define (display-map window poincare-map x y n)
  (plot-point window x y)
  (if (fix:> n 0)
      (poincare-map 
       x y
       (lambda (nx ny)
	 (display-map window poincare-map nx ny (fix:- n 1)))
       (lambda ()  
	   (newline)
	   (display "Illegal point: ")
	   (write (list x y))))))


;;; For example, first a simple problem.

(define ((standard-map K) x y continue fail)
  (let ((yp (flo:pv (flo:+ y (flo:* K (flo:sin x))))))
    (continue (flo:pv (flo:+ x yp)) yp)))

(define ((standard-map-inverse K) x y continue fail)
  (let ((xp (flo:pv (flo:- x y))))
    (continue xp (flo:pv (flo:- y (flo:* K (flo:sin xp)))))))

;;; This is the 0-2pi principal value:

(define (flo:pv x)
  (flo:- x (flo:* 2pi (flo:floor (flo:/ x 2pi)))))

#|
(define win (frame 0.0 2pi 0.0 2pi))
(explore-map win (standard-map 1.0) 1000)
(graphics-close win)
|#

#| 
(define ((T-pend m l g ys) local)
  (let ((t (time local))
        (theta (coordinate local))
        (thetadot (velocity local)))
    (let ((ysdot (D ys)))
      (* 1/2 m
         (+ (square (* l thetadot))
            (square (ysdot t))
            (* 2 (ysdot t) l (sin theta) thetadot))))))

(define ((V-pend m l g ys) local)
  (let ((t (time local))
        (theta (coordinate local)))
    (* m g (- (ys t) (* l (cos theta))))))

(define L-pend (- T-pend V-pend))

(define ((periodic-drive amplitude frequency phase) t)
  (* amplitude (cos (+ (* frequency t) phase))))

(define (L-periodically-driven-pendulum m l g a omega)
  (let ((ys (periodic-drive a omega 0)))
    (L-pend m l g ys)))

(define (H-pend-sysder m l g a omega)
  (Hamiltonian->state-derivative
    (Lagrangian->Hamiltonian
      (L-periodically-driven-pendulum m l g a omega))))

;;; The driven pendulum map is very interesting.
(set! *ode-integration-method* 'bulirsch-stoer)

(define (driven-pendulum-map mass l g a omega)
  (let ((map-period (/ 2pi omega)))
    (lambda (theta p-theta continue fail)
      (let ((ns
	     ((state-advancer H-pend-sysder
			      mass l g a omega)
	      (up			;state to start from
	       0.0			;t0 = phase/drive-freq.
	       theta
	       p-theta)
	      map-period			
	      1.0e-10)))
	(continue 
	 ((principal-value pi) (vector-ref ns 1))
	 (vector-ref ns 2))))))

;;; driven at twice the small-oscillation resonance.
;;; driven-pend-fig6.ps
(define win (frame -pi pi -10.0 10.0))

(let* ((m 1.0)				;m=1kg
       (l 1.0)				;l=1m
       (g 9.8)				;g=9.8m/s^2
       (small-amplitude-frequency (sqrt (/ g l)))
       (drive-amplitude 0.1)		;a=1/10 m
       (drive-frequency
	(* 2.0 small-amplitude-frequency)))
  (explore-map
   win
   (driven-pendulum-map m l g drive-amplitude drive-frequency)
   1000))

(graphics-close win)

;;; driven off resonance, at 4.2 times the natural frequency.
;;; driven-pend-fig7.ps

(define win (frame -pi pi -20 20))

(let* ((m 1.0)				;m=1kg
       (l 1.0)				;l=1m
       (g 9.8)				;g=9.8m/s^2
       (small-amplitude-frequency (sqrt (/ g l)))
       (drive-amplitude 0.05)		;a=1/20 m
       (drive-frequency
	(* 4.2 small-amplitude-frequency)))
  (explore-map
   win
   (driven-pendulum-map m l g drive-amplitude drive-frequency)
   1000))

(graphics-close win)

;;; Driven at high frequency, with larger amplitude.
;;; This shows the stable inverted pendulum island.
;;; driven-pend-fig8.ps

(define win (frame -pi pi -20 20))

(let* ((m 1.0)				;m=1kg
       (l 1.0)				;l=1m
       (g 9.8)				;g=9.8m/s^2
       (small-amplitude-frequency (sqrt (/ g l)))
       (drive-amplitude 0.2)		;a=1/20 m
       (drive-frequency
	(* 10.1 small-amplitude-frequency)))
  (explore-map
   win
   (driven-pendulum-map m l g drive-amplitude drive-frequency)
   1000))

(graphics-close win)
|#

#|
;;; An alternative, using construction and selection rather than
;;; continuations Here a map returns a new vector of x,y or an object
;;; describing the reason why it cannot.

(define (explore-map window poincare-map #!optional n)
  (define (iterate-map i x y)
    (if (fix:> i 0) 
	(let ((nxy (poincare-map x y)))
	  (if (vector? nxy)
	      (begin (plot-point window x y)
		     (iterate-map (fix:- i 1)
				  (vector-ref nxy 0)
				  (vector-ref nxy 1)))
	      (begin (newline)
		     (display "Illegal point: ")
		     (write (list  x y))
		     (button-loop x y))))
	(button-loop x y)))
  (define (button-loop ox oy)
    (pointer-coordinates 
     window 
     (lambda (x y button)
       (case button
	 ((0)
	  (write-line (list x y))
	  (display " started.")
	  (iterate-map n x y))
	 ((1) 
	  (write-line (list ox oy))
	  (display " continued.")
	  (iterate-map n ox oy))
	 ((2)
	  (write-line (list x y))
	  (display " hit.")
	  (button-loop ox oy))))))
  (if (default-object? n) (set! n 1000))
  (newline)
  (display "Left button starts a trajectory.")
  (newline)
  (display "Middle button continues a trajectory.")
  (newline)
  (display "Right button interrogates coordinates.")
  (button-loop 9. 9.))

;;; In this form the map is slightly different, but everything 
;;; else is the same.

(define ((standard-map K) x y)
  (let ((yp (flo:pv (flo:+ y (flo:* K (flo:sin x))))))
    (vector (flo:pv (flo:+ x yp)) yp)))

(define (driven-pendulum-map mass l g a omega)
  (let ((map-period (/ 2pi omega)))
    (lambda (theta p-theta)
      (let ((ns
	     ((state-advancer H-pend-sysder
			      mass l g a omega)
	      (up			;state to start from
	       0.0			;t0 = phase/drive-freq.
	       theta
	       p-theta)
	      map-period			
	      1.0e-10)))
	(vector
	 ((principal-value pi) (vector-ref ns 1))
	 (vector-ref ns 2))))))

|#
