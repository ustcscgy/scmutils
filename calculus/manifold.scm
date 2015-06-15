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

;;;; Manifolds are declared

(define (define-manifold manifold-name
	  manifold-type manifold-dimension
	  point-type point-dimension)
  (let ((point-prototype
	 (s:generate point-dimension 'up
		     (lambda (i)
		       (string->symbol
			(string-append
			 (symbol->string manifold-name)
			 "-"
			 (number->string i))))))
	(point-chains
	 (s:generate point-dimension 'up list))
	(coordinate-systems '()))
    (define (the-manifold m)
      (case m
	((name manifold-name) manifold-name)

	((type manifold-type) manifold-type)
	((dimension manifold-dimension) manifold-dimension)

	((point-type) point-type)
	((point-dimension) point-dimension)

	((typical-point) (typical-object point-prototype))
	((point-chains) point-chains)

	((new-coordinates)
	 (lambda (new)
	   (set! coordinate-systems
		 (cons new coordinate-systems))))

	(else
	 (let ((cs
		(find-matching-item coordinate-systems
		  (lambda (coordinate-system)
		    (eq? m (coordinate-system 'coordinate-system-name))))))
	   (if cs
	       cs
	       (error "Unknown message: manifold" manifold-name m))))))
    (if (environment-bound? generic-environment manifold-name)
	(write-line `(clobbering ,manifold-name)))
    (environment-define generic-environment manifold-name the-manifold)
    manifold-name))

(define (define-coordinate-system manifold coordinate-system-name
	  coordinate-prototype
	  coordinates->point point->coordinates)
  (let* ((access-chains
	  (s:map-chain (lambda (element chain) chain)
		       coordinate-prototype))
	 (dual-chains (flip-indices access-chains)))
    (let ((n (s:dimension coordinate-prototype))
	  (coordinate-functions #f)
	  (coordinate-basis-vector-fields #f)
	  (coordinate-basis-1form-fields #f))
      (if (not (fix:= n (manifold 'dimension)))
	  (error "Coordinate system does not have dimension of manifold"
		 coordinate-system-name (manifold 'manifold-name)))
      (define (the-coordinate-system m)
	(case m
	  ((manifold) manifold)
	  ((name coordinate-system-name) coordinate-system-name)
	  ((dimension coordinate-system-dimension)  n)
	  ((->point) coordinates->point)
	  ((->coords) point->coordinates)
	  ((typical-coords) (typical-object coordinate-prototype))
	  ((access-chains) access-chains)
	  ((dual-chains) dual-chains)
	  ((coordinate-functions)
	   (if (not coordinate-functions)
	       (set! coordinate-functions
		     (s:map/r (lambda (coordinate-name access-chain)
				(list coordinate-name
				      (compose (apply component access-chain)
					       point->coordinates)))
			      coordinate-prototype
			      access-chains)))
	   coordinate-functions)
	  ((coordinate-basis-vector-fields)
	   (if (not coordinate-basis-vector-fields)
	       (set! coordinate-basis-vector-fields
		     (s:map/r (lambda (coordinate-name access-chain)
				(let ((oname
				       (string->symbol
					(string-append "d/d"
						       (symbol->string coordinate-name)))))
				  (list oname
					(apply coordinate-basis-vector-field
					       the-coordinate-system
					       oname
					       access-chain))))
			      coordinate-prototype
			      access-chains)))
	   coordinate-basis-vector-fields)

	  ((coordinate-basis-1form-fields)
	   (if (not coordinate-basis-1form-fields)
	       (set! coordinate-basis-1form-fields
		     (s:map/r (lambda (coordinate-name access-chain)
				(let ((oname
				       (string->symbol
					(string-append "d"
						       (symbol->string coordinate-name)))))
				  (list oname
					(apply coordinate-basis-1form-field
					       the-coordinate-system
					       oname
					       access-chain))))
			      coordinate-prototype
			      access-chains)))
	   coordinate-basis-1form-fields)
	  ((coordinate-basis)
	   (make-basis (the-coordinate-system 'coordinate-basis-vector-fields)
		       (the-coordinate-system 'coordinate-basis-1form-fields)))
	  (else
	   (error "Unknown message: coordinate-system"
		  coordinate-system-name (manifold 'manifold-name)))))
      ((manifold 'new-coordinates) the-coordinate-system)))
  (list coordinate-system-name (manifold 'name)))






(define (install-coordinates coordinate-system)
  (define (install-symbol name value)
    (if (environment-bound? generic-environment name)
	(write-line `(clobbering ,name)))
    (environment-define generic-environment name value))
  (define (install-symbols s)
    (s:foreach (lambda (symval)
		 (install-symbol (car symval) (cadr symval)))
	       s))
  (install-symbols (coordinate-system 'coordinate-functions))
  (install-symbols (coordinate-system 'coordinate-basis-vector-fields))
  (install-symbols (coordinate-system 'coordinate-basis-1form-fields))
  (list (coordinate-system 'name) ((coordinate-system 'manifold) 'name)))

#|
(define-manifold 'Euclidean-plane Real 2 Real 2)

(define-coordinate-system Euclidean-plane 'rectangular (up 'x 'y)
  (lambda (coords)
    (if (and (up? coords) (fix:= (s:dimension coords) 2))
	coords
	(error "Bad coordinates: real-tuple" coords)))
  (lambda (point)
    (if (and (up? point) (fix:= (s:dimension point) 2))
	point
	(error "Bad point: real-tuple" point)))
  )

(define-coordinate-system Euclidean-plane 'polar (up 'r 'theta)
  (lambda (coords)
    (if (and (up? coords) (fix:= (s:dimension coords) 2))
	(let ((r (ref coords 0))
	      (theta (ref coords 1)))
	  (s:generate (s:dimension coords) 'up
		      (lambda (i)
			(cond ((= i 0) (* r (cos theta)))
			      ((= i 1) (* r (sin theta)))))))
	(error "Bad coordinates: polar" coords)))
  (lambda (point)
    (if (and (up? point) (fix:= (s:dimension point) 2))
	(let ((x (ref point 0))
	      (y (ref point 1)))
	  (s:generate (s:dimension point) 'up
		      (lambda (i)
			(cond ((= i 0) (sqrt (+ (square x) (square y))))
			      ((= i 1) (atan y x))
			      (else (ref point i))))))
	(error "Bad point: polar" point)))
  )

(install-coordinates (Euclidean-plane 'rectangular))

(install-coordinates (Euclidean-plane 'polar))


(define m (((Euclidean-plane 'rectangular) '->point) (up 'x0 'y0)))

(define circular (- (* x d/dy) (* y d/dx)))

(pec ((circular (+ (* 2 x) (* 3 y))) m))
#| Result:
(+ (* 3 x0) (* -2 y0))
|#

(pec ((circular theta) m))
#| Result:
1
|#

(pec ((dr circular) m))
#| Result:
0
|#
|#

#|
(define-manifold 'S2 real 2 real 3)

(define-coordinate-system S2 'colatitude-longitude (up 'theta 'phi)
  (lambda (coords)			;coordinates->point
    (if (and (up? coords) (fix:= (s:dimension coords) 2))
	(let ((theta (ref coords 0))
	      (phi   (ref coords 1)))
	  (up (* (sin theta) (cos phi)) ;x
	      (* (sin theta) (sin phi)) ;y
	      (* (cos theta))))		;z
	(error "Bad coordinates: S2" coords)))
  (lambda (point)			;point->coordinates 
    (if (and (up? point) (fix:= (s:dimension point) 3))
	(let ((x (ref point 0))
	      (y (ref point 1))
	      (z (ref point 2)))
	  ;;(assert (= 1
	  ;;           (+ (square x) (square y) (square z))))
	  (up (acos z)			;theta
	      (atan y x)))		;phi
	(error "Bad point: S2" point)))
  )
|#


