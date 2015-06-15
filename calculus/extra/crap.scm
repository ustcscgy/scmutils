#| -*-Scheme-*-

Copyright (C) 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994,
    1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005,
    2006, 2007, 2008, 2009, 2010, 2011 Massachusetts Institute of
    Technology

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

(define-manifold 'S2M Real 2 (up Real Real Real))

(define-coordinate-system S2M 'colatitude-longitude (up 'theta 'phi)
  (lambda (coords)			;coordinates->point
    (if (and (up? coords) (fix:= (s:dimension coords) 2))
	(let ((theta (ref coords 0))
	      (phi   (ref coords 1)))
	  (up (* (sin theta) (cos phi)) ;x
	      (* (sin theta) (sin phi)) ;y
	      (cos theta)))		;z
	(error "Bad coordinates: S2M" coords)))
  (lambda (point)			;point->coordinates 
    (if (and (up? point) (fix:= (s:dimension point) 3))
	(let ((x (ref point 0))
	      (y (ref point 1))
	      (z (ref point 2)))
	  ;;(assert (= 1
	  ;;           (+ (square x) (square y) (square z))))
	  (up (acos z)			;theta
	      (atan y x)))		;phi
	(error "Bad point: S2M" point))))


(define spherical-coordinates (S2M 'colatitude-longitude))


(define-manifold 'spacetime Real 4 (up Real Real Real Real))

(define-coordinate-system spacetime 'rectangular (up 't 'x 'y 'z)
  (lambda (coords)
    (if (and (up? coords) (fix:= (s:dimension coords) 4))
	coords
	(error "Bad coordinates: spacetime-rectangular" coords)))
  (lambda (point)
    (if (and (up? point) (fix:= (s:dimension point) 4))
	point
	(error "Bad point: spacetime-rectangular" point))))

(define-coordinate-system spacetime 'spherical (up 't 'r 'theta 'phi)
  (lambda (coords)
    (if (and (up? coords) (fix:= (s:dimension coords) 4))
	(let ((t (ref coords 0))
	      (r (ref coords 1))
	      (theta (ref coords 2))
	      (phi   (ref coords 3)))
	  (up t
	      (* r (sin theta) (cos phi)) 
	      (* r (sin theta) (sin phi)) 
	      (* r (cos theta))))
  	(error "Bad coordinates: spacetime-spherical" coords)))
  (lambda (point)
    (if (and (up? point) (fix:= (s:dimension point) 4))
	(let ((t (ref point 0))
	      (x (ref point 1))
	      (y (ref point 2))
	      (z (ref point 3)))
	  (let ((r (sqrt (+ (square x) (square y) (square z)))))
	    (up t
		r
		(acos (/ z r))
		(atan y x))))
	(error "Bad point: spacetime-spherical" point))))

(define spacetime-rectangular (spacetime 'rectangular))
(define spacetime-spherical (spacetime 'spherical))
