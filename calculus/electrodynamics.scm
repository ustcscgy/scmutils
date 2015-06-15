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

;;; Electrodynamics...

(define SR R4-rect)
(define-coordinates (up ct x y z) SR)
(define SR-basis (coordinate-system->basis SR))
(define an-event ((point SR) (up 'ct0 'x0 'y0 'z0)))

(define (g-Lorentz u v)
  (+ (* (dx u) (dx v))
     (* (dy u) (dy v))
     (* (dz u) (dz v))
     (* -1 (dct u) (dct v))))

(define L4-basis
  (coordinate-system->basis SR))


(define SR-star
  (Hodge-star g-Lorentz L4-basis))

(((SR-star
   (* (literal-manifold-function 'Bx SR)
      (wedge dy dz)))
  d/dct
  d/dx)
 an-event)
#| Result:
(Bx (up ct0 x0 y0 z0))
|#

;;; Fields E, B.  From MTW p.108

(define (Faraday Ex Ey Ez Bx By Bz)
  (+ (* Ex (wedge dx dct))
     (* Ey (wedge dy dct))
     (* Ez (wedge dz dct))
     (* Bx (wedge dy dz))
     (* By (wedge dz dx))
     (* Bz (wedge dx dy))))

(define (Maxwell Ex Ey Ez Bx By Bz)
  (+ (* -1 Bx (wedge dx dct))
     (* -1 By (wedge dy dct))
     (* -1 Bz (wedge dz dct))
     (* Ex (wedge dy dz))
     (* Ey (wedge dz dx))
     (* Ez (wedge dx dy))))

(((- (SR-star (Faraday 'Ex 'Ey 'Ez 'Bx 'By 'Bz))
     (Maxwell 'Ex 'Ey 'Ez 'Bx 'By 'Bz))
  (literal-vector-field 'u SR)
  (literal-vector-field 'v SR))
 an-event)
#| 0 |#


;;; **F + F = 0

(((+ ((compose SR-star SR-star) (Faraday 'Ex 'Ey 'Ez 'Bx 'By 'Bz))
     (Faraday 'Ex 'Ey 'Ez 'Bx 'By 'Bz))
  (literal-vector-field 'u SR)
  (literal-vector-field 'v SR))
 an-event)
#| 0 |#

;;; Defining the 4-current density J.

;;; Charge density is a manifold function.  Current density is a
;;; vector field having only spatial components.

(define (J charge-density Jx Jy Jz)
  (- (* (/ 1 :c) (+ (* Jx dx) (* Jy dy) (* Jz dz)))
     (* charge-density dct)))

(define rho (literal-manifold-function 'rho SR))

(define 4-current
  (J rho
     (literal-manifold-function 'Ix SR)
     (literal-manifold-function 'Iy SR)
     (literal-manifold-function 'Iz SR)))

(((d (SR-star 4-current))
  (literal-vector-field 'a SR)
  (literal-vector-field 'b SR)
  (literal-vector-field 'c SR)
  (literal-vector-field 'd SR))
 an-event)
#| Result:
;;; The charge conservation equations are too ugly to include.
|#

(((SR-star 4-current) d/dx d/dy d/dz) an-event)
#| 
(rho (up ct0 x0 y0 z0))
|#

(((SR-star 4-current) d/dct d/dy d/dz) an-event)
#| 
(/ (* -1 (Ix (up t0 x0 y0 z0))) :c)
|#

(((SR-star 4-current) d/dct d/dz d/dx) an-event)
#| 
(/ (* -1 (Iy (up t0 x0 y0 z0))) :c)
|#

(((SR-star 4-current) d/dct d/dx d/dy) an-event)
#| 
(/ (* -1 (Iz (up t0 x0 y0 z0))) :c)
|#

;;; Maxwell's equations in the form language are:
;;; dF=0, d(*F)=4pi *J
  
(define F
  (Faraday (literal-manifold-function 'Ex SR)
	   (literal-manifold-function 'Ey SR)
	   (literal-manifold-function 'Ez SR)
	   (literal-manifold-function 'Bx SR)
	   (literal-manifold-function 'By SR)
	   (literal-manifold-function 'Bz SR)))


;;; div B = 0
(((d F) d/dx d/dy d/dz) an-event)
#| Result: 
(+ (((partial 1) Bx) (up ct0 x0 y0 z0))
   (((partial 2) By) (up ct0 x0 y0 z0))
   (((partial 3) Bz) (up xt0 x0 y0 z0)))
|#

;;; curl E = -1/c dB/dt

(((d F) d/dct d/dy d/dz) an-event)
#|
(+ (((partial 0) Bx) (up ct0 x0 y0 z0))
   (((partial 2) Ez) (up ct0 x0 y0 z0))
   (* -1 (((partial 3) Ey) (up ct0 x0 y0 z0))))
|#

(((d F) d/dct d/dz d/dx) an-event)
#|
(+ (((partial 0) By) (up ct0 x0 y0 z0))
   (((partial 3) Ex) (up ct0 x0 y0 z0))
   (* -1 (((partial 1) Ez) (up ct0 x0 y0 z0))))
|#

(((d F) d/dct d/dx d/dy) an-event)
#|
(+ (((partial 0) Bz) (up ct0 x0 y0 z0))
   (((partial 1) Ey) (up ct0 x0 y0 z0))
   (* -1 (((partial 2) Ex) (up ct0 x0 y0 z0))))
|#

;;; div E = 4pi rho

(((- (d (SR-star F)) (* 4 :pi (SR-star 4-current)))
  d/dx d/dy d/dz)
 an-event)
#|
(+ (* -4 :pi (rho (up ct0 x0 y0 z0)))
   (((partial 1) Ex) (up ct0 x0 y0 z0))
   (((partial 2) Ey) (up ct0 x0 y0 z0))
   (((partial 3) Ez) (up ct0 x0 y0 z0)))
|#


;;; curl B = 1/c dE/dt + 4pi I

(((- (d (SR-star F)) (* 4 :pi (SR-star 4-current)))
  d/dct d/dy d/dz)
 an-event)
#|
(+ (((partial 0) Ex) (up ct0 x0 y0 z0))
   (* -1 (((partial 2) Bz) (up ct0 x0 y0 z0)))
   (((partial 3) By) (up ct0 x0 y0 z0))
   (/ (* 4 :pi (Ix (up ct0 x0 y0 z0))) :c))
|#

(((- (d (SR-star F)) (* 4 :pi (SR-star 4-current)))
  d/dct d/dz d/dx)
 an-event)
#|
(+ (((partial 0) Ey) (up ct0 x0 y0 z0))
   (* -1 (((partial 3) Bx) (up ct0 x0 y0 z0)))
   (((partial 1) Bz) (up ct0 x0 y0 z0))
   (/ (* 4 :pi (Iy (up ct0 x0 y0 z0))) :c))
|#

(((- (d (SR-star F)) (* 4 :pi (SR-star 4-current)))
  d/dct d/dx d/dy)
 an-event)
#|
(+ (((partial 0) Ez) (up ct0 x0 y0 z0))
   (* -1 (((partial 1) By) (up ct0 x0 y0 z0)))
   (((partial 2) Bx) (up ct0 x0 y0 z0))
   (/ (* 4 :pi (Iz (up ct0 x0 y0 z0))) :c))
|#

