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

;;; From covariant-derivative.scm

(define M (make-manifold S^2 2 3))
(define spherical
  (coordinate-system-at 'spherical 'north-pole M))
(define-coordinates (up theta phi) spherical)
(define-coordinates t the-real-line)
(define spherical-basis (coordinate-system->basis spherical))

(define G-S2-1
  (make-Christoffel
   (let ((zero  (lambda (point) 0))) 
     (down (down (up zero zero)
                 (up zero (/ 1 (tan theta))))
           (down (up zero (/ 1 (tan theta)))
                 (up (- (* (sin theta) (cos theta))) zero))))
   spherical-basis))


(define gamma:N->M
  (compose (spherical '->point)
           (up (lambda (t) :pi/2)
               (lambda (t) t))
           (the-real-line '->coords)))

(define basis-over-gamma
  (basis->basis-over-map gamma:N->M spherical-basis))


(define ((s u) t)
  (up t (u t)))

(define (g gamma cartan)
  (let ((omega
	 ((Cartan->forms (Cartan->Cartan-over-map cartan gamma))
	  ((differential gamma) d/dt))))
    (define (the-state-derivative state)
      (let ((t ((the-real-line '->point) (ref state 0)))
            (u (ref state 1)))
	(up 1 (*  -1 (omega t) u))))
    the-state-derivative))


(define (Lg g)
  (make-operator
   (lambda (h)
     (lambda (state)
       (* ((D h) state)
	  (g state))))
   `(Lg ,g)))


(series:for-each print-expression
  (((exp (Lg (g gamma:N->M (Christoffel->Cartan G-S2-1))))
    (component 1))
   (up 0 (up 0 1)))
  4)
(up 0 1)
(up -6.123031769111886e-17 0.)
(up 0. -1.8745759022776718e-33)
(up 3.8260292677525877e-50 0.)
;Value: ...


(series:for-each print-expression
  (((exp (Lg (g gamma:N->M (Christoffel->Cartan G-S2-1))))
    (component 1))
   (up 0 (up 1 0)))
  4)
(up 1 0)
(up 0. 6.123031769111886e-17)
(up -1.8745759022776718e-33 0.)
(up 0. -3.8260292677525877e-50)
;Value: ...

(series:for-each print-expression
  (((exp (Lg (g gamma:N->M (Christoffel->Cartan G-S2-1))))
    (component 1))
   (up 0 (up 1 1)))
  4)
(up 1 1)
(up -6.123031769111886e-17 6.123031769111886e-17)
(up -1.8745759022776718e-33 -1.8745759022776718e-33)
(up 3.8260292677525877e-50 -3.8260292677525877e-50)
;Value: ...



(define gamma1:N->M
  (compose (spherical '->point)
           (up (lambda (t) (+ ':pi/2 t))
               (lambda (t) t))
           (the-real-line '->coords)))
;Value: gamma1:N->M



(series:for-each (compose pp
			  (lambda (x) (eval x user-generic-environment))
			  simplify
			  expression)
  (((exp (Lg (g gamma1:N->M (Christoffel->Cartan G-S2-1))))
    (component 1))
   (up 0 (up 0 1)))
  6)
#(0 1)
#(-6.123031769111886e-17 6.123031769111886e-17)
#(.5 -.5)
#(7.143537063963867e-17 2.041010589703962e-17)
#(-.2916666666666667 -.08333333333333333)
#(-2.3471621781595565e-17 -1.4797326775353724e-17)
;Value: ...

(define gamma2:N->M
  (compose (spherical '->point)
           (up (lambda (t) 1)
               (lambda (t) t))
           (the-real-line '->coords)))
;Value: gamma2:N->M

(series:for-each (compose pp
			  (lambda (x) (eval x user-generic-environment))
			  simplify
			  expression)
  (((exp (Lg (g gamma2:N->M (Christoffel->Cartan G-S2-1))))
    (component 1))
   (up 0 (up 0 1)))
  6)
#(0 1)
#(-.4546487134128409 0.)
#(0. -.14596329086321444)
#(.02212067413215491 0.)
#(0. 3.5508803799365564e-3)
#(-3.22880639244211e-4 0.)
;Value: ...

(series:for-each (compose pp
			  (lambda (x) (eval x user-generic-environment))
			  simplify
			  expression)
  (((exp (Lg (g gamma2:N->M (Christoffel->Cartan G-S2-1))))
    (component 1))
   (up 0 (up 1 0)))
  6)
#(1 0)
#(0. .6420926159343308)
#(-.14596329086321444 0.)
#(0. -.03124065042024832)
#(3.5508803799365564e-3 0.)
#(0. 4.5599881440467076e-4)
;Value: ...

(series:for-each (compose pp
			  (lambda (x) (eval x user-generic-environment))
			  simplify
			  expression)
  (((exp (Lg (g gamma2:N->M (Christoffel->Cartan G-S2-1))))
    (component 1))
   (up 0 (up 1 0)))
  6)

(define gamma2:N->M
  (compose (spherical '->point)
           (up (lambda (t) 1/100)
               (lambda (t) t))
           (the-real-line '->coords)))
;Value: gamma2:N->M

(series:for-each (compose pp
			  (lambda (x) (eval x user-generic-environment))
			  simplify
			  expression)
  (((exp (* :pi/2 (Lg (g gamma2:N->M (Christoffel->Cartan G-S2-1)))))
    (component 1))
   (up 0 (up 1 0)))
  6)
#(1 0)
#(0. 157.07439665682676)
#(-1.2335771841934366 0.)
#(0. -64.58779731227044)
#(.2536187782271013 0.)
#(0. 7.9674033141726985)
;Value: ...

(series:for-each (compose pp
			  (lambda (x) (eval x user-generic-environment))
			  simplify
			  expression)
  (((exp (* :pi/2 (Lg (g gamma2:N->M (Christoffel->Cartan G-S2-1)))))
    (component 1))
   (up 0 (up 1 0)))
  8)
#(1 0)
#(0. 157.07439665682676)
#(-1.2335771841934366 0.)
#(0. -64.58779731227044)
#(.2536187782271013 0.)
#(0. 7.9674033141726985)
#(-2.0857222553597822e-2 0.)
#(0. -.46801937836336255)
;Value: ...

(+ 157.07439665682676 -64.58779731227044 7.9674033141726985 -.46801937836336255)
;Value: 99.98598328036566

(+ 1 -1.2335771841934366 .2536187782271013 -2.0857222553597822e-2)
;Value: -8.15628519933144e-4

(define (transform coords)
  (let ((colat (ref coords 0))
	(long (ref coords 1)))
    (let ((x (* (sin colat) (cos long)))
	  (y (* (sin colat) (sin long)))
	  (z (cos colat)))
      (let ((vp ((rotate-x (/ :pi 4)) (up x y z))))
	(let ((colatp (acos (ref vp 2)))
	      (longp (atan (ref vp 1) (ref vp 0))))
	  (up colatp longp))))))

(define (transform-inverse coords)
  (let ((colat (ref coords 0))
	(long (ref coords 1)))
    (let ((x (* (sin colat) (cos long)))
	  (y (* (sin colat) (sin long)))
	  (z (cos colat)))
      (let ((vp ((rotate-x (- (/ :pi 4))) (up x y z))))
	(let ((colatp (acos (ref vp 2)))
	      (longp (atan (ref vp 1) (ref vp 0))))
	  (up colatp longp))))))

(define gamma3:N->M
  (let ()
    (define (coords t) (transform (up pi/2 t)))
    (compose (spherical '->point)
	     coords
	     (the-real-line '->coords))))

#|
(series:for-each (compose pp
			  (lambda (x) (eval x user-generic-environment))
			  simplify	
			  expression)
  (((exp (* 1/10 ':pi (Lg (g gamma3:N->M (Christoffel->Cartan G-S2-1)))))
    (component 1))
   (up 0 (up (/ (sqrt 2) 2) (/ (sqrt 2) 2))))
  6)
|#



(define ((gg) state)
  ((g gamma3:N->M (Christoffel->Cartan G-S2-1)) state))


((state-advancer gg)
 (up 0 (up (/ (sqrt 2) 2) (/ (sqrt 2) 2)))
 pi/2)
;Value: #(1.570796326794897 #(.9999999999998719 2.584433346486528e-13))

(* ((D transform) (up pi/2 0)) (up 1 0))
;Value: #(.7071067811865476 .7071067811865475)

(* ((D transform) (up pi/2 pi/2)) (up 1 0))
;Value: #(1. 8.659274570719353e-17)


((state-advancer gg)
 (up 0 (up (/ (sqrt 2) 2) (/ (sqrt 2) 2)))
 1)
;Value: #(1. #(.8797941764813766 .5914444826417807))

(* ((D transform) (up pi/2 1)) (up 1 0))
;Value: #(.879794176481595 .5914444826416471)

