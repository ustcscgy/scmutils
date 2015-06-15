#| -*-Scheme-*-

Copyright (C) 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994,
    1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005,
    2006, 2007, 2008, 2009, 2010 Massachusetts Institute of Technology

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


(define (S^n-stereographic offset)
  (lambda (manifold)
    (let ((offset (offset (manifold 'dimension))))
      (define (me m)
	(case m
	  ((check-coordinates)
	   (lambda (coords)
	     (let ((coords (- coords offset)))
	       (and (up? coords)
		    (fix:= (s:dimension coords) (manifold 'dimension))
		    (let ((remaining-coords (cdr (up-structure->list
		    coords))))
		      (every (lambda (coord)
			       (or (not (number? coord))
				   (not (= coord 0))))
			     remaining-coords))))))
	  ((check-point)
	   (lambda (point) (my-manifold-point? point manifold)))
	  ((coords->point)
	   (lambda (coords)
	     (if (not ((me 'check-coordinates) coords))
		 (error "Bad coordinates: S^n" coords))
	     (let ((offset-coords (- coords offset))
		   (n (manifold 'dimension)))
	       (let* ((delta (dot-product offset-coords offset-coords))
		      (x0 (/ (- 1 delta) (+ 1 delta)))
		      (pt (s:generate
			   (fix:+ n 1)
			   'up
			   (lambda (i)
			     (if (fix:= i 1) x0
				 (/ (* 2 (ref offset-coords
					 (if (fix:= i 0) 0
					     (- i 1))))
				    (+ 1 delta)))))))
		 (make-manifold-point pt manifold me coords)))))

	  ((point->coords)
	   (lambda (point)
	     (if (not ((me 'check-point) point))
		 (error "Bad point: S^n" point))
	     (let ((pt (manifold-point-representation point)))
	       (s:generate
		(manifold 'dimension)
		'up
		(lambda (i)
		  (/ (ref pt (if (= i 0) 0
				 (+ i 1)))
		     (+ 1 (ref pt 1))))))))
	  (else
	   (error "S^n: Bad message" m))))
      me)))



;; (define (S^n-stereographic pole)
;;   (lambda (manifold)
;;     (let ((offset (offset (manifold 'dimension))))
;;       (define (me m)
;; 	(case m
;; 	  ((check-coordinates)
;; 	   (lambda (coords)
;; 	     (and (up? coords)
;; 		  (fix:= (s:dimension coords) (manifold 'dimension))
;; 		  (let ((remaining-coords (cdr (up-structure->list
;; coords))))
;; 		    (every (lambda (coord)
;; 			     (or (not (number? coord))
;; 				 (not (= coord 0))))
;; 			   remaining-coords)))))
;; 	  ((check-point)
;; 	   (lambda (point) (my-manifold-point? point manifold)))
;; 	  ((coords->point)
;; 	   (lambda (coords)
;; 	     (if (not ((me 'check-coordinates) coords))
;; 		 (error "Bad coordinates: S^n" coords))
;; 	     (let ((n (manifold 'dimension)))
;; 	       (let* ((delta (dot-product coords coords))
;; 		      (x0 (/ (- 1 delta) (+ 1 delta)))
;; 		      (pt (s:generate
;; 			   (fix:+ n 1)
;; 			   'up
;; 			   (lambda (i)
;; 			     (if (fix:= i 1) x0
;; 				 (/ (* 2 (ref coords
;; 					 (if (fix:= i 0) 0
;; 					     (- i 1))))
;; 				    (+ 1 delta)))))))
;; 		 (make-manifold-point pt manifold me coords)))))
;; 	  ((point->coords)
;; 	   (lambda (point)
;; 	     (if (not ((me 'check-point) point))
;; 		 (error "Bad point: S^n" point))
;; 	     (let ((pt (manifold-point-representation point)))
;; 	       (s:generate
;; 		(manifold 'dimension)
;; 		'up
;; 		(lambda (i)
;; 		  (/ (ref pt (if (= i 0) 0
;; 				 (+ i 1)))
;; 		     (+ 1 (ref pt 1))))))))
;; 	  (else
;; 	   (error "S^n: Bad message" m))))
;;       me)))

(attach-coordinate-system 'stereographic 'north-pole S^n-type
  (S^n-stereographic
   (lambda (n) (s:generate n 'up (lambda (i) 0)))))

(attach-coordinate-system 'stereographic 'tilted S^n-type
  (S^n-stereographic
   (lambda (n) (s:generate n 'up (lambda (i) 0)))))


#| The Circle

(define S1 (make-manifold S^n-type 1))
(define circular (coordinate-system-at 'spherical 'north-pole S1))
(define slope (coordinate-system-at 'stereographic 'north-pole S1))

(define-coordinates (up theta) circular)
(define-coordinates (up s) slope)

(define f (compose cos theta))

(define mtheta  ((circular '->point) (up 'theta)))
(define ms  ((slope '->point) (up 's)))

(pe ((slope '->coords) mtheta))
(pe ((circular '->coords) mtheta))
(pe ((circular '->coords) ms))
(se ((slope '->coords) ms))

(pe ((d/dtheta f) mtheta))
(pe ((d/dtheta f) ms))
(se ((d/ds theta) ms))

(se ((d/dtheta theta) mtheta))
(se ((d/dtheta s) ms))
(se ((d/ds theta) mtheta))

(se (down (up  ((d/ds theta) mtheta)
 	       ((d/dtheta s) mtheta))
 	  (up ((d/ds theta) ms)
 	      (simplify ((d/dtheta s) ms)))))

|#

#| The Riemann Sphere

(define S2 (make-manifold S^n-type 2))
(define latlong (coordinate-system-at 'spherical 'north-pole S2))
(define riemann (coordinate-system-at 'stereographic 'north-pole S2))

(define-coordinates (up theta1 theta2) latlong)
(define-coordinates (up x y) riemann)



(define mlatlong  ((latlong '->point) (up 'theta1 'theta2)))
(define mriemann  ((riemann '->point) (up 'x 'y)))

(pe ((latlong '->coords) mlatlong))
(se ((riemann '->coords) mlatlong))
(pe ((riemann '->coords) mriemann))
(se ((latlong '->coords) mriemann))

(define f (compose cos (* theta1 theta2)))

(pe ((d/dtheta1 f) mlatlong))
(pe ((d/dtheta2 f) mlatlong))

(se ((d/dtheta1 x) mlatlong))
(pe ((d/dx theta1) mlatlong))
(pe ((d/dy theta1) mlatlong))

(pe ((d/dx theta1) mlatlong))


(se (down (up  ((d/ds theta) mtheta)
 	       ((d/dtheta s) mtheta))
 	  (up ((d/ds theta) ms)
 	      (simplify ((d/dtheta s) ms)))))

|#
