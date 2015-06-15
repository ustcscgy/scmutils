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

;;;; Test sequence for code from /usr/local/scmutils/dg/calculus.txt

;;; Caution, this test defines X, which is also used as a 
;;; part of the function type constructor in other code.

(define R2 (make-manifold R^n 2))

(define U (patch 'origin R2))

(define R2-rect (coordinate-system-at 'rectangular 'origin R2))

(define R2-rect (coordinate-system 'rectangular U))

(define R2-polar (coordinate-system 'polar/cylindrical U))

(define R2-rect-chi (R2-rect '->coords))
(define R2-rect-chi-inverse (R2-rect '->point))
(define R2-polar-chi (R2-polar '->coords))
(define R2-polar-chi-inverse (R2-polar '->point))


(print-expression 
  ((compose R2-polar-chi R2-rect-chi-inverse) 
   (up 'x0 'y0)))
#;(up (sqrt (+ (expt x0 2) (expt y0 2))) (atan y0 x0))

(print-expression 
  ((compose R2-rect-chi R2-polar-chi-inverse) 
   (up 'r0 'theta0)))
#;(up (* r0 (cos theta0)) (* r0 (sin theta0)))


(print-expression 
  ((D (compose R2-rect-chi R2-polar-chi-inverse))
   (up 'r0 'theta0)))
#;(down (up (cos theta0) (sin theta0))
        (up (* -1 r0 (sin theta0)) (* r0 (cos theta0))))

(define f (literal-manifold-function 'f-rect R2-rect))

(define R2->R (-> (UP Real Real) Real))

(define f
  (compose (literal-function 'f-rect R2->R)
           R2-rect-chi))

(define R2-point (R2-rect-chi-inverse (up 'x0 'y0)))

(define polar-point 
  (R2-polar-chi-inverse 
   (up (sqrt (+ (square 'x0) (square 'y0))) 
       (atan 'y0 'x0))))

(print-expression (f R2-point))
#;(f-rect (up x0 y0))


(print-expression (f polar-point))
#;(f-rect (up x0 y0))

(define-coordinates (up x y) R2-rect)
(define-coordinates (up r theta) R2-polar)

(print-expression (x (R2-rect-chi-inverse (up 'x0 'y0))))
#;x0

(print-expression (x (R2-polar-chi-inverse (up 'r0 'theta0))))
#;(* r0 (cos theta0))

(print-expression (r (R2-polar-chi-inverse (up 'r0 'theta0))))
#;r0

(print-expression (r (R2-rect-chi-inverse (up 'x0 'y0))))
#;(sqrt (+ (expt x0 2) (expt y0 2)))

(print-expression (theta (R2-rect-chi-inverse (up 'x0 'y0))))
#;(atan y0 x0)

(define h (+ (* x (square r)) (cube y)))

(print-expression (h R2-point))
#;(+ (expt x0 3) (* x0 (expt y0 2)) (expt y0 3))

(print-expression (h (R2-polar-chi-inverse (up 'r0 'theta0))))
#;(+ (* (expt r0 3) (expt (sin theta0) 3)) (* (expt r0 3) (cos theta0)))


(print-expression
 ((D (compose h R2-polar-chi-inverse))
  (up 'r0 'theta0)))
#;
(down
 (+ (* 3 (expt r0 2) (expt (sin theta0) 3)) (* 3 (expt r0 2) (cos theta0)))
 (+ (* 3 (expt r0 3) (cos theta0) (expt (sin theta0) 2))
    (* -1 (expt r0 3) (sin theta0))))


(define w (literal-vector-field 'v R2-rect))

(define v
  (components->vector-field
   (up (literal-function 'v^0 R2->R)
       (literal-function 'v^1 R2->R))
   R2-rect))

(print-expression
  ((v (literal-manifold-function 'f-rect R2-rect)) R2-point))
#;(+ (* (v^0 (up x0 y0)) (((partial 0) f-rect) (up x0 y0)))
     (* (v^1 (up x0 y0)) (((partial 1) f-rect) (up x0 y0))))

(print-expression
 (((coordinatize v R2-rect) (literal-function 'f-rect R2->R))
  (up 'x0 'y0)))
#;(+ (* (v^0 (up x0 y0)) (((partial 0) f-rect) (up x0 y0)))
     (* (v^1 (up x0 y0)) (((partial 1) f-rect) (up x0 y0))))

(define-coordinates (up x y) R2-rect)
(define-coordinates (up r theta) R2-polar)

(print-expression ((d/dx (square r)) R2-point))
#;(* 2 x0)

(print-expression 
  (((+ d/dx (* 2 d/dy)) (+ (square r) (* 3 x))) R2-point))
#;(+ 3 (* 2 x0) (* 4 y0))


(define J (- (* x d/dy) (* y d/dx)))

(series:for-each print-expression
                 (((exp (* 'a J)) R2-rect-chi)
                  ((R2-rect '->point) (up 1 0)))
                 6)
#|
(up 1 0)
(up 0 a)
(up (* -1/2 (expt a 2)) 0)
(up 0 (* -1/6 (expt a 3)))
(up (* 1/24 (expt a 4)) 0)
(up 0 (* 1/120 (expt a 5)))
;Value: ...
|#

(print-expression
 ((((evolution 6) 'a J) R2-rect-chi)
  ((R2-rect '->point) (up 1 0))))
#;(up (+ 1 (* -1/720 (expt a 6)) (* 1/24 (expt a 4)) (* -1/2 (expt a 2)))
      (+ (* 1/120 (expt a 5)) (* -1/6 (expt a 3)) a))


(define omega
  (components->1form-field
   (down (literal-function 'a\_0 R2->R)
         (literal-function 'a\_1 R2->R))
   R2-rect))

(pe ((omega (down d/dx d/dy)) R2-point))
#;(down (a_0 (up x0 y0)) (a_1 (up x0 y0)))

(eq? 'a\_0 'a_0)
;Value: #t


(print-expression ((dx d/dy) R2-point))
#;0

(print-expression ((dx d/dx) R2-point))
#;1

(print-expression ((dx J) R2-point))
#;(* -1 y0)

(print-expression ((dy J) R2-point))
#;x0

(print-expression ((dr J) R2-point))
#;0

(print-expression ((dtheta J) R2-point))
#;1

(define f (literal-manifold-function 'f R2-rect))

(print-expression (((- J d/dtheta) f) R2-point))
#;0

(print-expression
  ((omega (literal-vector-field 'v R2-rect)) R2-point))
#;(+ (* (v^0 (up x0 y0)) (a_0 (up x0 y0)))
     (* (v^1 (up x0 y0)) (a_1 (up x0 y0))))

(define e0 
  (+ (* (literal-manifold-function 'e0x R2-rect) d/dx)
     (* (literal-manifold-function 'e0y R2-rect) d/dy)))

(define e1 
  (+ (* (literal-manifold-function 'e1x R2-rect) d/dx)
     (* (literal-manifold-function 'e1y R2-rect) d/dy)))

(define e-vector-basis (down e0 e1))

(define e-dual-basis 
  (vector-basis->dual e-vector-basis R2-polar))

(print-expression ((e-dual-basis e-vector-basis) R2-point))
#;(up (down 1 0) (down 0 1))

(define v
  (* (up (literal-manifold-function 'bx R2-rect)
         (literal-manifold-function 'by R2-rect))
     e-vector-basis))

(print-expression ((e-dual-basis v) R2-point))
#;(up (bx (up x0 y0)) (by (up x0 y0)))

(let* ((polar-basis (coordinate-system->basis R2-polar))
       (polar-vector-basis (basis->vector-basis polar-basis))
       (polar-dual-basis (basis->1form-basis polar-basis))
       (f (literal-manifold-function 'f R2-rect)))
  (print-expression
    ((- ((commutator e0 e1) f)
        (* (- (e0 (polar-dual-basis e1))
              (e1 (polar-dual-basis e0)))
           (polar-vector-basis f)))
     R2-point)))
#;0


(define R3 (make-manifold R^n 3))
(define R3-rect (coordinate-system-at 'rectangular 'origin R3))
(define-coordinates (up x y z) R3-rect)
(define R3-point ((R3-rect '->point) (up 'x0 'y0 'z0)))
(define R3-cyl (coordinate-system-at 'polar/cylindrical 'origin R3))
(define-coordinates (up r theta zeta) R3-cyl)
(define R3->R (-> (UP Real Real Real) Real))
(define g (literal-manifold-function 'g R3-rect))

(define Jz (- (* x d/dy) (* y d/dx)))
(define Jx (- (* y d/dz) (* z d/dy)))
(define Jy (- (* z d/dx) (* x d/dz)))

(print-expression (((+ (commutator Jx Jy) Jz) g) R3-point))
#;0

(print-expression (((+ (commutator Jy Jz) Jx) g) R3-point))
#;0

(print-expression (((+ (commutator Jz Jx) Jy) g) R3-point))
#;0

(define v (+ (* 'a d/dx) (* 'b d/dy)))
(define w (+ (* 'c d/dx) (* 'd d/dy)))

(print-expression (((wedge dx dy) v w) R3-point))
#;(+ (* a d) (* -1 b c))

(define u (+ (* 'a d/dx) (* 'b d/dy) (* 'c d/dz)))
(define v (+ (* 'd d/dx) (* 'e d/dy) (* 'f d/dz)))
(define w (+ (* 'g d/dx) (* 'h d/dy) (* 'i d/dz)))

(print-expression
  (((wedge dx dy dz) u v w) R3-point))
#;(+ (* a e i) (* -1 a f h) (* -1 b d i) (* b f g) (* c d h) (* -1 c e g))

(print-expression
  (- (((wedge dx dy dz) u v w) R3-point)
     (determinant
      (matrix-by-rows (list 'a 'b 'c)
                      (list 'd 'e 'f)
                      (list 'g 'h 'i)))))
#;0

(define a (literal-manifold-function 'alpha R3-rect))
(define b (literal-manifold-function 'beta R3-rect))
(define c (literal-manifold-function 'gamma R3-rect))

(define theta (+ (* a dx) (* b dy) (* c dz)))

(define X (literal-vector-field 'X R3-rect))
(define Y (literal-vector-field 'Y R3-rect))

(print-expression
 (((- (d theta)
      (+ (wedge (d a) dx)
         (wedge (d b) dy)
         (wedge (d c) dz)))
   X Y)
  R3-point))
#;0

(define omega
  (+ (* a (wedge dy dz))
     (* b (wedge dz dx))
     (* c (wedge dx dy))))

(define Z (literal-vector-field 'Z R3-rect))

(print-expression
 (((- (d omega)
      (+ (wedge (d a) dy dz)
         (wedge (d b) dz dx)
         (wedge (d c) dx dy)))
   X Y Z)
  R3-point))
#;0

(print-expression
 (((d (d theta)) X Y Z) R3-point))
#;0

(define v (literal-vector-field 'v R2-rect))
(define w (literal-vector-field 'w R2-rect))

(define alpha (literal-function 'alpha R2->R))
(define beta (literal-function 'beta R2->R))

(define R2-rect-basis (coordinate-system->basis R2-rect))

(print-expression
 (let ((dx (ref (basis->1form-basis R2-rect-basis) 0))
       (dy (ref (basis->1form-basis R2-rect-basis) 1)))
   (((- (d (+ (* (compose alpha R2-rect-chi) dx)
              (* (compose beta R2-rect-chi) dy)))
        (* (compose (- ((partial 0) beta)
                       ((partial 1) alpha))
                    R2-rect-chi)
           (wedge dx dy)))
     v w)
    R2-point)))
#;0

(define domega
  (* (+ (d/dx a) (d/dy b) (d/dz c))
     (wedge dx dy dz)))

(print-expression
 (((- (d omega) domega) X Y Z) R3-point))
#;0

(define S2 (make-manifold S^2 2 3))
(define S2-spherical
  (coordinate-system-at 'spherical 'north-pole S2))
(define-coordinates (up theta phi) S2-spherical)
(define S2-basis (coordinate-system->basis S2-spherical))

(define mu 
  (compose (S2-spherical '->point)
           (up (literal-function 'alpha)
               (literal-function 'beta))
           (the-real-line '->coords)))

(define S2-basis-over-mu 
  (basis->basis-over-map mu S2-basis))

(define h
  (compose (literal-function 'h R2->R)
           (S2-spherical '->coords)))


(print-expression
  (((basis->vector-basis S2-basis-over-mu) h)
   ((the-real-line '->point) 't0)))
#;(down (((partial 0) h) (up (alpha t0) (beta t0)))
        (((partial 1) h) (up (alpha t0) (beta t0))))

(print-expression
  (((basis->1form-basis S2-basis-over-mu)
    (basis->vector-basis S2-basis-over-mu))
   ((the-real-line '->point) 't0)))
#;(up (down 1 0) (down 0 1))

(define-coordinates t the-real-line)

(print-expression
  (((basis->1form-basis S2-basis-over-mu)
    ((differential mu) d/dt))
   ((the-real-line '->point) 't0)))
#;(up ((D alpha) t0) ((D beta) t0))

(define-coordinates (up x y z) R3-rect)
(define-coordinates (up u v) R2-rect)
(define X2 (literal-vector-field 'X R2-rect))
(define Y2 (literal-vector-field 'Y R2-rect))

(define mu
  (compose (R3-rect '->point)
           (up (literal-function 'mux R2->R)
               (literal-function 'muy R2->R)
               (literal-function 'muz R2->R))
           (R2-rect '->coords)))

(define f (literal-manifold-function 'f R3-rect))

(print-expression
 (((- ((pullback mu) (d f)) (d ((pullback mu) f)))
   X2)
  R2-point))
#;0

(define theta (literal-1form-field 'theta R3-rect))

(print-expression
 (((- ((pullback mu) (d theta)) (d ((pullback mu) theta)))
   X2 Y2)
  R2-point))
#;0

;;;; Lie Derivative

(define Jz (- (* x d/dy) (* y d/dx)))


(series:for-each print-expression
  ((((exp (* 'a (Lie-derivative Jz))) d/dy) 
    (literal-manifold-function 'f R3-rect))
   ((R3-rect '->point) (up 1 0 0)))
  3)
#|					;Takes a very long time
(((partial 1) f) (up 1 0 0))
(* a (((partial 0) f) (up 1 0 0)))
(* -1/2 (expt a 2) (((partial 1) f) (up 1 0 0)))
;Value: ...
|#
 
(define V (literal-vector-field 'V R3-rect))

(print-expression
  (((- ((Lie-derivative V) (d theta))
       (d ((Lie-derivative V) theta)))
    X Y)
   R3-point))
#;0

(print-expression
  (((- ((Lie-derivative V) (d omega))
       (d ((Lie-derivative V) omega)))
    X Y Z)
   R3-point))
#;0

(print-expression
 ((((- (commutator (Lie-derivative X) (Lie-derivative Y))
       (Lie-derivative (commutator X Y)))
    theta)
   Z)
  R3-point))
#;0

(print-expression
 ((((- (commutator (Lie-derivative X) (Lie-derivative Y))
       (Lie-derivative (commutator X Y)))
    omega)
   Z V)
  R3-point))
#;0

;;;; Interior Product

(define X (literal-vector-field 'X R3-rect))
(define Y (literal-vector-field 'Y R3-rect))
(define Z (literal-vector-field 'Z R3-rect))

(define alpha (literal-manifold-function 'alpha R3-rect))
(define beta (literal-manifold-function 'beta R3-rect))
(define gamma (literal-manifold-function 'gamma R3-rect))

(define omega
  (+ (* alpha (wedge dx dy))
     (* beta (wedge dy dz))
     (* gamma (wedge dz dx))))

(define ((L1 X) omega)
  (+ ((interior-product X) (d omega))
     (d ((interior-product X) omega))))

(print-expression
  ((- (((Lie-derivative X) omega) Y Z)
      (((L1 X) omega) Y Z))
   R3-point))
#;0

;;;; Covariant Derivative of Vector Fields

(define R2-rect-basis (coordinate-system->basis R2-rect))
(define R2-polar-basis (coordinate-system->basis R2-polar))
(define-coordinates (up x y) R2-rect)
(define-coordinates (up r theta) R2-polar)

(define R2-rect-Christoffel
  (make-Christoffel
   (let ((zero (lambda (m) 0)))
     (down (down (up zero zero)
                 (up zero zero))
           (down (up zero zero)
                 (up zero zero))))
   R2-rect-basis))

(define R2-polar-Christoffel
  (make-Christoffel
   (let ((zero (lambda (m) 0)))
     (down (down (up zero zero)
                 (up zero (/ 1 r)))
           (down (up zero (/ 1 r))
                 (up (* -1 r) zero))))
   R2-polar-basis))


(define R2-rect-Cartan
  (Christoffel->Cartan R2-rect-Christoffel))

(define R2-polar-Cartan
  (Christoffel->Cartan R2-polar-Christoffel))


(define f
  (compose (literal-function 'f-rect R2->R) R2-rect-chi))

(print-expression
 (((((covariant-derivative R2-rect-Cartan) d/dx)
    J)
   f)
  R2-point))
#;(((partial 1) f-rect) (up x0 y0))

(print-expression ((d/dy f) R2-point))
#;(((partial 1) f-rect) (up x0 y0))


(print-expression
 (((((covariant-derivative R2-polar-Cartan) d/dx)
    J)
   f)
  R2-point))
#;(((partial 1) f-rect) (up x0 y0))


(define v (literal-vector-field 'v R2-rect))
(define w (literal-vector-field 'w R2-rect))

(print-expression
 (((((- (covariant-derivative R2-rect-Cartan)
        (covariant-derivative R2-polar-Cartan))
     v)
    w)
   f)
  R2-point))
#;0


;;;; Parallel Transport


(define M (make-manifold R^n 2))
(define spherical (coordinate-system-at 'rectangular 'origin M))
(define-coordinates (up theta phi) spherical)
(define spherical-basis (coordinate-system->basis spherical))

(define G-S2-1
  (make-Christoffel
   (let ((zero  (lambda (point) 0))) 
     (down (down (up zero zero)
                 (up zero (/ 1 (tan theta))))
           (down (up zero (/ 1 (tan theta)))
                 (up (- (* (sin theta) (cos theta))) zero))))
   spherical-basis))


(define v (literal-vector-field 'v spherical))

(define f (literal-manifold-function 'f spherical))

(define omega (Christoffel->Cartan G-S2-1))

(define ((((F v) t) f) m)
  (define ((shift t) m)
    ((spherical '->point)
     (up (+ (theta m) t)
	 (phi m))))
  (define (sin-theta m)
    (sin (theta m)))
  (let ((vphi ((dphi v) ((shift t) m)))
	(vtheta ((dtheta v) ((shift t) m))))
    (((+ (* vtheta d/dtheta)
	 (* vphi
	    d/dphi 
	    (/ (sin-theta ((shift t) m))
	       (sin-theta m))))
      f)
     m)))


(pec ((((D (F v)) 0) f) ((spherical '->point) (up 'theta0 'phi0))))

#| Result:
(+
 (* (((partial 0) v^0) (up theta0 phi0)) (((partial 0) f) (up theta0 phi0)))
 (* (((partial 0) v^1) (up theta0 phi0)) (((partial 1) f) (up theta0 phi0)))
 (/ (* (v^1 (up theta0 phi0)) (cos theta0) (((partial 1) f) (up theta0 phi0)))
    (sin theta0)))
|#

(pec (- ((((D (F v)) 0) f) ((spherical '->point) (up 'theta0 'phi0)))
       (((((covariant-derivative omega) d/dtheta) v)
	 f)
	((spherical '->point) (up 'theta0 'phi0)))))

#| Result:
0
|#


(define M (make-manifold R^n 2))
(define spherical (coordinate-system-at 'rectangular 'origin M))
(define-coordinates (up theta phi) spherical)
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
           (up (literal-function 'alpha)
               (literal-function 'beta))
           (the-real-line '->coords)))

(define basis-over-gamma
  (basis->basis-over-map gamma:N->M spherical-basis))

(define w
  (basis-components->vector-field
   (up (compose (literal-function 'w0)
                (the-real-line '->coords))
       (compose (literal-function 'w1)
                (the-real-line '->coords)))
   (basis->vector-basis basis-over-gamma)))

(define sphere-Cartan (Christoffel->Cartan G-S2-1))

(print-expression
 (s:map/r 
  (lambda (omega)
    ((omega
      (((covariant-derivative-over-map sphere-Cartan gamma:N->M) 
        d/dt) 
       w))
     ((the-real-line '->point) 'tau)))
  (basis->1form-basis basis-over-gamma)))
#;(up
   (+ (* -1 (w1 tau) ((D beta) tau) (cos (alpha tau)) (sin (alpha tau))) ((D w0) tau))
   (+ (/ (* (w1 tau) (cos (alpha tau)) ((D alpha) tau)) (sin (alpha tau)))
      (/ (* ((D beta) tau) (cos (alpha tau)) (w0 tau)) (sin (alpha tau)))
      ((D w1) tau)))


;;;; Geodesic Motion

(define gamma:N->M
  (compose (spherical '->point)
           (up (literal-function 'alpha)
               (literal-function 'beta))
           (the-real-line '->coords)))

(print-expression
 (((((covariant-derivative-over-map  sphere-Cartan gamma:N->M) d/dt)
    ((differential gamma:N->M) d/dt))
   (spherical '->coords))
  ((the-real-line '->point) 't)))
#;
(up
 (+ (* -1 (expt ((D beta) t) 2) (cos (alpha t)) (sin (alpha t)))
    (((expt D 2) alpha) t))
 (+ (/ (* 2 ((D beta) t) (cos (alpha t)) ((D alpha) t)) (sin (alpha t)))
    (((expt D 2) beta) t)))


(define (Lfree s)
  (* 1/2 (square (velocity s))))

(define (sphere->R3 s)
  (let ((q (coordinate s)))
    (let ((theta (ref q 0)) (phi (ref q 1)))
      (up (* (sin theta) (cos phi))
          (* (sin theta) (sin phi))
          (cos theta)))))

(define Lsphere
  (compose Lfree (F->C sphere->R3)))


(print-expression
  (((Lagrange-equations Lsphere)
    (up (literal-function 'alpha)
        (literal-function 'beta)))
   't))
#;
(down
 (+ (* -1 (expt ((D beta) t) 2) (cos (alpha t)) (sin (alpha t)))
    (((expt D 2) alpha) t))
 (+ (* 2 ((D beta) t) (cos (alpha t)) (sin (alpha t)) ((D alpha) t))
    (* (expt (sin (alpha t)) 2) (((expt D 2) beta) t))))


;;; Curvature


(define sphere-Cartan (Christoffel->Cartan G-S2-1))

(print-expression
 (((Riemann (covariant-derivative sphere-Cartan))
   dphi d/dtheta d/dphi d/dtheta)
  ((spherical '->point) (up 'theta0 'phi0))))
#;1

(define R4 (make-manifold R^n 4))
(define R4-rect (coordinate-system-at 'rectangular 'origin R4))
(define states R4-rect)
(define-coordinates (up Theta Phi w0 w1) states)

(define (G v)
  (let ((alphadot (dTheta v)) (betadot (dPhi v)))
    (+ v
       (* (compose (* sin cos) Theta) betadot w1 d/dw0)
       (* -1
          (compose (/ cos sin) Theta)
          (+ (* w0 betadot) (* w1 alphadot))
          d/dw1))))


(define Gu (G d/dTheta))
(define Gv (G d/dPhi))

(define (initial-state initial-coords w)
  (let ((Theta0 (ref initial-coords 0))
        (Phi0 (ref initial-coords 1)))
    (let ((m ((spherical '->point) (up Theta0 Phi0))))
      ((states '->point)
       (up Theta0 Phi0 ((dtheta w) m) ((dphi w) m))))))

(print-expression
 ((dw0 (commutator Gu Gv))
  (initial-state (up 'Theta0 'Phi0) d/dphi)))
#;(* -1 (expt (sin Theta0) 2))

(print-expression
 ((dw1 (commutator Gu Gv))
  (initial-state (up 'Theta0 'Phi0) d/dtheta)))
#;1

;;; Torsion

(let ((nabla (covariant-derivative sphere-Cartan)))
  (for-each
   (lambda (x)
     (for-each
      (lambda (y)
	(print-expression
	 ((((torsion nabla) x y) 
	   (literal-manifold-function 'f spherical))
	  ((spherical '->point) (up 'theta0 'phi0)))))
      (list  d/dtheta d/dphi)))
   (list  d/dtheta d/dphi)))
#|
0
0
0
0
|#

;;;; Metrics

(define ((g-sphere R) u v)
  (* (square R)
     (+ (* (dtheta u) (dtheta v))
        (* (compose (square sin) theta)
           (dphi u)
           (dphi v)))))

(print-expression
 ((Christoffel->symbols
   (metric->Christoffel-1 (g-sphere 'R) spherical-basis))
  ((spherical '->point) (up 'theta0 'phi0))))
#;
(down
 (down (down 0 0)
       (down 0 (* (expt R 2) (sin theta0) (cos theta0))))
 (down (down 0 (* (expt R 2) (sin theta0) (cos theta0)))
       (down (* -1 (expt R 2) (sin theta0) (cos theta0)) 0)))



(print-expression
 ((Christoffel->symbols
   (metric->Christoffel-2 (g-sphere 'R) spherical-basis))
  ((spherical '->point) (up 'theta0 'phi0))))
#;
(down
 (down (up 0 0)
       (up 0 (/ (cos theta0) (sin theta0))))
 (down (up 0 (/ (cos theta0) (sin theta0)))
       (up (* -1 (sin theta0) (cos theta0)) 0)))


;;;; Electrodynamics

(define SR R4-rect)
(define-coordinates (up t x y z) SR)
(define an-event ((SR '->point) (up 't0 'x0 'y0 'z0)))
(define c 'c)         ; We like units.

(define (g-Lorentz u v)
  (+ (* (dx u) (dx v))
     (* (dy u) (dy v))
     (* (dz u) (dz v))
     (* -1 (square c) (dt u) (dt v))))

(define SR-vector-basis
  (down (* (/ 1 c) d/dt) d/dx d/dy d/dz))

(define SR-1form-basis
  (up (* c dt) dx dy dz))

(define SR-basis
  (make-basis SR-vector-basis
              SR-1form-basis))

(print-expression
 ((SR-1form-basis SR-vector-basis)
  an-event))
#;(up (down 1 0 0 0) (down 0 1 0 0) (down 0 0 1 0) (down 0 0 0 1))

(define (Faraday Ex Ey Ez Bx By Bz)
  (+ (* Ex c (wedge dx dt))
     (* Ey c (wedge dy dt))
     (* Ez c (wedge dz dt))
     (* Bx (wedge dy dz))
     (* By (wedge dz dx))
     (* Bz (wedge dx dy))))

(define (Maxwell Ex Ey Ez Bx By Bz)
  (+ (* -1 Bx c (wedge dx dt))
     (* -1 By c (wedge dy dt))
     (* -1 Bz c (wedge dz dt))
     (* Ex (wedge dy dz))
     (* Ey (wedge dz dx))
     (* Ez (wedge dx dy))))

(define SR-star
  (Hodge-star g-Lorentz SR-basis))


(print-expression
 (((- (SR-star (Faraday 'Ex 'Ey 'Ez 'Bx 'By 'Bz))
      (Maxwell 'Ex 'Ey 'Ez 'Bx 'By 'Bz))
   (literal-vector-field 'u SR)
   (literal-vector-field 'v SR))
  an-event))
#;0

(define (J charge-density Ix Iy Iz)
  (- (* (/ 1 c) (+ (* Ix dx) (* Iy dy) (* Iz dz)))
     (* charge-density c dt)))

(define F
  (Faraday (literal-manifold-function 'Ex SR)
           (literal-manifold-function 'Ey SR)
           (literal-manifold-function 'Ez SR)
           (literal-manifold-function 'Bx SR)
           (literal-manifold-function 'By SR)
           (literal-manifold-function 'Bz SR)))

(define 4-current
  (J (literal-manifold-function 'rho SR)
     (literal-manifold-function 'Ix SR)
     (literal-manifold-function 'Iy SR)
     (literal-manifold-function 'Iz SR)))

(print-expression
 (((d F) d/dx d/dy d/dz) an-event))
#;
(+ (((partial 1) Bx) (up t0 x0 y0 z0))
   (((partial 2) By) (up t0 x0 y0 z0))
   (((partial 3) Bz) (up t0 x0 y0 z0)))

(print-expression
 (((d F) (* (/ 1 c) d/dt) d/dy d/dz) an-event))
#;
(+ (((partial 2) Ez) (up t0 x0 y0 z0))
   (* -1 (((partial 3) Ey) (up t0 x0 y0 z0)))
   (/ (((partial 0) Bx) (up t0 x0 y0 z0)) c))


(print-expression
 (((d F) (* (/ 1 c) d/dt) d/dz d/dx) an-event))
#;
(+ (((partial 3) Ex) (up t0 x0 y0 z0))
   (* -1 (((partial 1) Ez) (up t0 x0 y0 z0)))
   (/ (((partial 0) By) (up t0 x0 y0 z0)) c))

(print-expression
 (((d F) (* (/ 1 c) d/dt) d/dx d/dy) an-event))
#;
(+ (((partial 1) Ey) (up t0 x0 y0 z0))
   (* -1 (((partial 2) Ex) (up t0 x0 y0 z0)))
   (/ (((partial 0) Bz) (up t0 x0 y0 z0)) c))


(print-expression
 (((- (d (SR-star F)) (* '4pi (SR-star 4-current)))
   d/dx d/dy d/dz)
  an-event))
#;
(+ (* -1 4pi (rho (up t0 x0 y0 z0)))
   (((partial 1) Ex) (up t0 x0 y0 z0))
   (((partial 2) Ey) (up t0 x0 y0 z0))
   (((partial 3) Ez) (up t0 x0 y0 z0)))

(print-expression
 (((- (d (SR-star F)) (* '4pi (SR-star 4-current)))
   (* (/ 1 c) d/dt) d/dy d/dz)
  an-event))
#;
(+ (/ (* 4pi (Ix (up t0 x0 y0 z0))) c)
   (* -1 (((partial 2) Bz) (up t0 x0 y0 z0)))
   (((partial 3) By) (up t0 x0 y0 z0))
   (/ (((partial 0) Ex) (up t0 x0 y0 z0)) c))

(print-expression
 (((- (d (SR-star F)) (* '4pi (SR-star 4-current)))
   (* (/ 1 c) d/dt) d/dz d/dx)
  an-event))
#;
(+ (/ (* 4pi (Iy (up t0 x0 y0 z0))) c)
   (* -1 (((partial 3) Bx) (up t0 x0 y0 z0)))
   (((partial 1) Bz) (up t0 x0 y0 z0))
   (/ (((partial 0) Ey) (up t0 x0 y0 z0)) c))

(print-expression
 (((- (d (SR-star F)) (* '4pi (SR-star 4-current)))
   (* (/ 1 c) d/dt) d/dx d/dy)
  an-event))
#;
(+ (/ (* 4pi (Iz (up t0 x0 y0 z0))) c)
   (* -1 (((partial 1) By) (up t0 x0 y0 z0)))
   (((partial 2) Bx) (up t0 x0 y0 z0))
   (/ (((partial 0) Ez) (up t0 x0 y0 z0)) c))


;;; Testing Covariant Derivative on forms.

(define (Gijk i j k)
  (literal-manifold-function
   (string->symbol
    (string-append "G^"
		   (number->string i)
		   "_"
		   (number->string j)
		   (number->string k)))
   R2-rect))
				    
(define G
  (down (down (up (Gijk 0 0 0)
		  (Gijk 1 0 0))
	      (up (Gijk 0 1 0)
		  (Gijk 1 1 0)))
	(down (up (Gijk 0 0 1)
		  (Gijk 1 0 1))
	      (up (Gijk 0 1 1)
		  (Gijk 1 1 1)))))

(define X (literal-vector-field 'X R2-rect))

(define Y (literal-vector-field 'Y R2-rect))

(define omega (literal-1form-field 'omega R2-rect))

(define q_0 (up 'q_x 'q_y))

(define m_0
  ((R2-rect '->point) q_0))

(define F (literal-manifold-function 'F R2-rect))

(let* ((CF (Christoffel->Cartan
	    (make-Christoffel G
			      (coordinate-system->basis R2-rect))))
       (D_x ((covariant-derivative CF) X)))
  (pe (- (+ (((D_x omega) Y) m_0)
	    ((omega (D_x Y)) m_0))
	 ((D_x (omega Y)) m_0))))
#;0


(define tau (literal-1form-field 'tau R2-rect))

(define Z (literal-vector-field 'Z R2-rect))

(let* ((CF (Christoffel->Cartan
	    (make-Christoffel G
			      (coordinate-system->basis R2-rect))))
       (D_x ((covariant-derivative CF) X)))
  (pe (- (((D_x (wedge omega tau)) Y Z) m_0)
	 (+ (((wedge omega (D_x tau)) Y Z) m_0)
	    (((wedge (D_x omega) tau) Y Z) m_0)))))
#;0

(let* ((CF (Christoffel->Cartan
	    (make-Christoffel G
			      (coordinate-system->basis R2-polar))))
       (D_x ((covariant-derivative CF) X)))
  (pe (- (((D_x (wedge omega tau)) Y Z) m_0)
	 (+ (((wedge omega (D_x tau)) Y Z) m_0)
	    (((wedge (D_x omega) tau) Y Z) m_0)))))
#;0

;;;; Answers below, for checking and timing
#|
;;; On zohar (ThinkPad T42)

cat /proc/cpuinfo
processor	: 0
vendor_id	: GenuineIntel
cpu family	: 6
model		: 13
model name	: Intel(R) Pentium(R) M processor 1.80GHz
stepping	: 6
cpu MHz		: 600.000
cache size	: 2048 KB
fdiv_bug	: no
hlt_bug		: no
f00f_bug	: no
coma_bug	: no
fpu		: yes
fpu_exception	: yes
cpuid level	: 2
wp		: yes
flags		: fpu vme de pse tsc msr mce cx8 sep mtrr pge mca cmov pat clflush dts acpi mmx fxsr sse sse2 ss tm pbe up est tm2
bogomips	: 1196.95
clflush size	: 64

Image saved on Monday October 27, 2008 at 10:07:48 PM
  Release 7.7.90.+                 || Microcode 15.1
  Runtime 15.7                     || SF 4.41
  LIAR/i386 4.118                  || Edwin 3.116
  ScmUtils Mechanics . Summer 2008 || SOS 1.8
  IMAIL 1.21

(gc-flip)
;Value: 11906842

(show-time
 (lambda ()
   (load "/usr/local/scmutils/src/calculus/test")))
;Loading "/usr/local/scmutils/src/calculus/test.scm"... 
(up (sqrt (+ (expt x0 2) (expt y0 2))) (atan y0 x0))
(up (* r0 (cos theta0)) (* r0 (sin theta0)))
(down (up (cos theta0) (sin theta0))
      (up (* -1 r0 (sin theta0)) (* r0 (cos theta0))))
(f-rect (up x0 y0))
(f-rect (up x0 y0))
x0
(* r0 (cos theta0))
r0
(sqrt (+ (expt x0 2) (expt y0 2)))
(atan y0 x0)
(+ (expt x0 3) (* x0 (expt y0 2)) (expt y0 3))
(+ (* (expt r0 3) (expt (sin theta0) 3)) (* (expt r0 3) (cos theta0)))
(down
 (+ (* 3 (expt r0 2) (expt (sin theta0) 3)) (* 3 (expt r0 2) (cos theta0)))
 (+ (* 3 (expt r0 3) (cos theta0) (expt (sin theta0) 2))
    (* -1 (expt r0 3) (sin theta0))))
(+ (* (v^0 (up x0 y0)) (((partial 0) f-rect) (up x0 y0)))
   (* (v^1 (up x0 y0)) (((partial 1) f-rect) (up x0 y0))))
(+ (* (v^0 (up x0 y0)) (((partial 0) f-rect) (up x0 y0)))
   (* (v^1 (up x0 y0)) (((partial 1) f-rect) (up x0 y0))))
(* 2 x0)
(+ 3 (* 2 x0) (* 4 y0))
(up 1 0)
(up 0 a)
(up (* -1/2 (expt a 2)) 0)
(up 0 (* -1/6 (expt a 3)))
(up (* 1/24 (expt a 4)) 0)
(up 0 (* 1/120 (expt a 5)))
(up (+ 1 (* -1/720 (expt a 6)) (* 1/24 (expt a 4)) (* -1/2 (expt a 2)))
    (+ (* 1/120 (expt a 5)) (* -1/6 (expt a 3)) a))
(down (a_0 (up x0 y0)) (a_1 (up x0 y0)))
0
1
(* -1 y0)
x0
0
1
0
(+ (* (v^0 (up x0 y0)) (a_0 (up x0 y0))) (* (v^1 (up x0 y0)) (a_1 (up x0 y0))))
(up (down 1 0) (down 0 1))
(up (bx (up x0 y0)) (by (up x0 y0)))
0
0
0
0
(+ (* a d) (* -1 b c))
(+ (* a e i) (* -1 a f h) (* -1 b d i) (* b f g) (* c d h) (* -1 c e g))
0
0
0
0
0
0
(down (((partial 0) h) (up (alpha t0) (beta t0)))
      (((partial 1) h) (up (alpha t0) (beta t0))))
(up (down 1 0) (down 0 1))
(up ((D alpha) t0) ((D beta) t0))
0
0
(((partial 1) f) (up 1 0 0))
(* a (((partial 0) f) (up 1 0 0)))
(* -1/2 (expt a 2) (((partial 1) f) (up 1 0 0)))
0
0
0
0
(((partial 1) f-rect) (up x0 y0))
(((partial 1) f-rect) (up x0 y0))
(((partial 1) f-rect) (up x0 y0))
0
(+ (/ (* (cos theta0) (v^1 (up theta0 phi0)) (((partial 1) f) (up theta0 phi0))) (sin theta0))
   (* (((partial 1) f) (up theta0 phi0)) (((partial 0) v^1) (up theta0 phi0)))
   (* (((partial 0) f) (up theta0 phi0)) (((partial 0) v^0) (up theta0 phi0))))
0
(up
 (+ (* -1 (w1 tau) ((D beta) tau) (cos (alpha tau)) (sin (alpha tau))) ((D w0) tau))
 (+ (/ (* (w1 tau) (cos (alpha tau)) ((D alpha) tau)) (sin (alpha tau)))
    (/ (* ((D beta) tau) (cos (alpha tau)) (w0 tau)) (sin (alpha tau)))
    ((D w1) tau)))
(up
 (+ (* -1 (expt ((D beta) t) 2) (cos (alpha t)) (sin (alpha t)))
    (((expt D 2) alpha) t))
 (+ (/ (* 2 ((D beta) t) (cos (alpha t)) ((D alpha) t)) (sin (alpha t)))
    (((expt D 2) beta) t)))
(down
 (+ (* -1 (expt ((D beta) t) 2) (cos (alpha t)) (sin (alpha t)))
    (((expt D 2) alpha) t))
 (+ (* 2 ((D beta) t) (cos (alpha t)) (sin (alpha t)) ((D alpha) t))
    (* (expt (sin (alpha t)) 2) (((expt D 2) beta) t))))
1
(* -1 (expt (sin Theta0) 2))
1
0
0
0
0
(down
 (down (down 0 0) (down 0 (* (expt R 2) (cos theta0) (sin theta0))))
 (down (down 0 (* (expt R 2) (cos theta0) (sin theta0)))
       (down (* -1 (expt R 2) (cos theta0) (sin theta0)) 0)))
(down
 (down (up 0 0) (up 0 (/ (cos theta0) (sin theta0))))
 (down (up 0 (/ (cos theta0) (sin theta0)))
       (up (* -1 (cos theta0) (sin theta0)) 0)))
(up (down 1 0 0 0) (down 0 1 0 0) (down 0 0 1 0) (down 0 0 0 1))
0
(+ (((partial 1) Bx) (up t0 x0 y0 z0))
   (((partial 2) By) (up t0 x0 y0 z0))
   (((partial 3) Bz) (up t0 x0 y0 z0)))
(+ (((partial 2) Ez) (up t0 x0 y0 z0))
   (* -1 (((partial 3) Ey) (up t0 x0 y0 z0)))
   (/ (((partial 0) Bx) (up t0 x0 y0 z0)) c))
(+ (((partial 3) Ex) (up t0 x0 y0 z0))
   (* -1 (((partial 1) Ez) (up t0 x0 y0 z0)))
   (/ (((partial 0) By) (up t0 x0 y0 z0)) c))
(+ (((partial 1) Ey) (up t0 x0 y0 z0))
   (* -1 (((partial 2) Ex) (up t0 x0 y0 z0)))
   (/ (((partial 0) Bz) (up t0 x0 y0 z0)) c))
(+ (* -1 4pi (rho (up t0 x0 y0 z0)))
   (((partial 1) Ex) (up t0 x0 y0 z0))
   (((partial 2) Ey) (up t0 x0 y0 z0))
   (((partial 3) Ez) (up t0 x0 y0 z0)))
(+ (/ (* 4pi (Ix (up t0 x0 y0 z0))) c)
   (* -1 (((partial 2) Bz) (up t0 x0 y0 z0)))
   (((partial 3) By) (up t0 x0 y0 z0))
   (/ (((partial 0) Ex) (up t0 x0 y0 z0)) c))
(+ (/ (* 4pi (Iy (up t0 x0 y0 z0))) c)
   (* -1 (((partial 3) Bx) (up t0 x0 y0 z0)))
   (((partial 1) Bz) (up t0 x0 y0 z0))
   (/ (((partial 0) Ey) (up t0 x0 y0 z0)) c))
(+ (/ (* 4pi (Iz (up t0 x0 y0 z0))) c)
   (* -1 (((partial 1) By) (up t0 x0 y0 z0)))
   (((partial 2) Bx) (up t0 x0 y0 z0))
   (/ (((partial 0) Ez) (up t0 x0 y0 z0)) c))
0
0
0
;... done
;;; Binah: 9 December 2009, unmemoized derivative
;; process time: 17950 (16740 RUN + 1210 GC); real time: 17945
;;; Binah: 9 December 2009, memoized derivative
;; process time: 36930 (35370 RUN + 1560 GC); real time: 36928
;;; Earlier
;process time: 67340 (58580 RUN + 8760 GC); real time: 67583
;Unspecified return value


(show-memoizer-statistics)
(44 8 memoized-simplify)
(0 0 complex-rules)
(0 0 exp-expand)
(0 0 exp-contract)
(0 0 exp->sincos)
(0 0 sincos->exp2)
(0 0 sincos->exp1)
(42 28 sincos-random)
(38 24 flush-obvious-ones)
(39 23 split-high-degree-sines)
(38 24 split-high-degree-cosines)
(0 0 cos^2->sin^2)
(38 24 sin^2->cos^2)
(0 0 contract-expt-trig)
(0 0 contract-multiangle)
(0 0 expand-multiangle)
(10 32 angular-parity)
(8 23 sincos->trig)
(25 129 trig->sincos)
(12 45 canonicalize-partials)
(0 0 log-expand)
(0 0 log-contract)
(0 0 logexp->specfun)
(0 0 specfun->logexp)
(4 13 sqrt-contract)
(24 22 sqrt-expand)
(0 0 special-trig)
(2 6 triginv)
(12 12 simsqrt)
(142 137 miscsimp)
(142 137 magsimp)
(0 0 logexp)
(0 0 diff:cosh)
(0 0 diff:sinh)
(21 114 diff:atan2)
(0 0 diff:atan1)
(0 0 diff:acos)
(0 0 diff:asin)
(23 33 diff:cos)
(35 69 diff:sin)
(0 0 diff:log)
(0 0 diff:exp)
(0 0 diff:expt)
(0 0 diff:power)
(49 114 diff:sqrt)
(0 0 diff:invert)
(1 262 diff:negate)
(45 37 diff:/)
(313 16974 diff:*)
(243 798 diff:-)
(980 7843 diff:+)
;Value: done
|#