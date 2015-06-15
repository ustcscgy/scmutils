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

;;;; nabla_X V  =  covariant derivative of V wrt X
;;;   V is a vector field, X is a vector field

;;; More complete covariant derivative procedure

(define (covariant-derivative Cartan)
  (assert (Cartan? Cartan))
  (define (nabla X)
    (define (nabla_X V)
      (cond ((function? V) (X V))
	    ((vector-field? V)
	     (((covariant-derivative-vector Cartan) X) V))
	    ((form-field? V)
	     (((covariant-derivative-form Cartan) X) V))
	    (else
	     (error "Bad input -- covariant-derivative"))))
    (make-operator nabla_X `(nabla ,(diffop-name X))))
  nabla)


(define (covariant-derivative-vector Cartan)
  (let ((basis (Cartan->basis Cartan))
        (Cartan-forms (Cartan->forms Cartan)))
    (let ((vector-basis (basis->vector-basis basis))
          (1form-basis (basis->1form-basis basis)))
      (lambda (V)
	(lambda (U)
	  (let ((u-components (1form-basis U)))
	    (let ((deriv-components
		   (+ (V u-components)
		      (* (Cartan-forms V) u-components))))
	      (define (the-derivative f)
		(* (vector-basis f) deriv-components))
	      (procedure->vector-field the-derivative
	       `((nabla ,(diffop-name V)) ,(diffop-name U))))))))))
	  

(define (((covariant-derivative-form Cartan) V) tau)
  (let ((k (get-rank tau))
	(nabla_V ((covariant-derivative-vector Cartan) V)))
    (procedure->nform-field
     (lambda vectors
       (assert (= k (length vectors)))
       (- (V (apply tau vectors))
	  (sigma (lambda (i)
		   (apply tau
			  (list-with-substituted-coord vectors i
				(nabla_V (list-ref vectors i)))))
		 0 (- k 1))))
     k
     `((nabla ,(diffop-name V)) ,(diffop-name tau)))))

;;; also nabla V (X), where nabla V is covariant differential
;;; nabla V(X)

(define (((covariant-differential Cartan) V) X)
  (((covariant-derivative Cartan) X) V))


(define (Cartan->Christoffel Cartan)
  (assert (Cartan? Cartan))
  (let ((basis (Cartan->basis Cartan))
	(Cartan-forms (Cartan->forms Cartan)))
    (make-Christoffel
     (s:map/r Cartan-forms
	      (basis->vector-basis basis))
     basis)))

(define (Christoffel->Cartan Christoffel)
  (assert (Christoffel? Christoffel))
  (let ((basis (Christoffel->basis Christoffel))
	(Christoffel-symbols
	 (Christoffel->symbols Christoffel)))
    (make-Cartan
     (* Christoffel-symbols (basis->1form-basis basis))
     basis)))

;;; Constructors and Selectors

(define (Cartan-transform cartan basis-prime)
  (let ((basis (Cartan->basis cartan)) ;; tuple of basis vectors
	(forms (Cartan->forms cartan))
	(prime-dual-basis (basis->1form-basis basis-prime))
	(prime-vector-basis (basis->vector-basis basis-prime)))
    (let ((vector-basis (basis->vector-basis basis))
	  (1form-basis (basis->1form-basis basis)))
      (let ((J-inv (s:map/r 1form-basis prime-vector-basis))
	    (J (s:map/r prime-dual-basis vector-basis)))
	(let ((omega-prime-forms
	       (procedure->1form-field
		(lambda (u)
		  (+ (* J (u J-inv))
		     (* J (* (forms u) J-inv)))))))
	  (make-Cartan omega-prime-forms basis-prime))))))
		  

(define (make-Cartan forms basis)
  (list '*Cartan* forms basis))

(define (Cartan? thing)
  (and (pair? thing)
       (eq? (car thing) '*Cartan*)))

(define (Cartan->forms thing) (cadr thing))

(define (Cartan->basis thing) (caddr thing))


(define (make-Christoffel symbols basis)
  (list '*Christoffel* symbols basis))

(define (Christoffel? thing)
  (and (pair? thing)
       (eq? (car thing) '*Christoffel*)))

(define (Christoffel->symbols thing) (cadr thing))

(define (Christoffel->basis thing) (caddr thing))

#|
;;; Fun with Christoffel symbols.

(install-coordinates R2-rect (up 'x 'y))

(define R2-rect-basis
  (coordinate-system->basis R2-rect))
(define R2-rect-point
  ((R2-rect '->point) (up 'x0 'y0)))

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
	      
(clear-arguments)
(suppress-arguments '((up x0 y0)))

(pec (G R2-rect-point)
    (compose arg-suppressor simplify))
#| Result:
(down (down (up G^0_00 G^1_00) (up G^0_10 G^1_10))
      (down (up G^0_01 G^1_01) (up G^0_11 G^1_11)))
|#

(define CG (make-Christoffel G R2-rect-basis))

(define CF (Christoffel->Cartan CG))

(pec (((Cartan->forms CF) (literal-vector-field 'X R2-rect))
      R2-rect-point))
#| Result:
(down
 (up
  (+ (* (G^0_00 (up x0 y0)) (X^0 (up x0 y0)))
     (* (G^0_01 (up x0 y0)) (X^1 (up x0 y0))))
  (+ (* (G^1_00 (up x0 y0)) (X^0 (up x0 y0)))
     (* (G^1_01 (up x0 y0)) (X^1 (up x0 y0)))))
 (up
  (+ (* (G^0_10 (up x0 y0)) (X^0 (up x0 y0)))
     (* (G^0_11 (up x0 y0)) (X^1 (up x0 y0))))
  (+ (* (G^1_10 (up x0 y0)) (X^0 (up x0 y0)))
     (* (G^1_11 (up x0 y0)) (X^1 (up x0 y0))))))
|#

(pec ((Christoffel->symbols
       (Cartan->Christoffel (Christoffel->Cartan CG)))
      R2-rect-point)
    (compose arg-suppressor simplify))
#| Result:
(down (down (up G^0_00 G^1_00) (up G^0_10 G^1_10))
      (down (up G^0_01 G^1_01) (up G^0_11 G^1_11)))
|#

;; Transformation of Cartan to polar leaves covariant derivative
;; invariant.

(pec (((((- (covariant-derivative CF)
	    (covariant-derivative
	     (Cartan-transform CF (R2-polar 'coordinate-basis))))
	 (literal-vector-field 'A R2-rect))
	(literal-vector-field 'B R2-polar))
       (literal-scalar-field 'f R2-polar))
      R2-rect-point))

#| Result:
0
|#

;; Example from the text: 

(define-coordinates (up x y) R2-rect)
(define-coordinates (up r theta) R2-polar)

(define v (literal-vector-field 'v R2-rect))
(define w (literal-vector-field 'w R2-rect))
(define f (literal-manifold-function 'f R2-rect))

(define R2-rect-basis (coordinate-system->basis R2-rect))
(define R2-polar-basis (coordinate-system->basis R2-polar))

(define R2-rect-Christoffel
  (make-Christoffel
   (let ((zero (lambda (m) 0)))
     (down (down (up zero zero)
                 (up zero zero))
           (down (up zero zero)
                 (up zero zero))))
   R2-rect-basis))

(define R2-rect-Cartan
  (Christoffel->Cartan R2-rect-Christoffel))


(define R2-polar-Christoffel
  (make-Christoffel
   (let ((zero (lambda (m) 0)))
     (down (down (up zero zero)
                 (up zero (/ 1 r)))
           (down (up zero (/ 1 r))
                 (up (* -1 r) zero))))
   R2-polar-basis))

(define R2-polar-Cartan
  (Christoffel->Cartan R2-polar-Christoffel))


(pec
 (((((- (covariant-derivative R2-rect-Cartan)
        (covariant-derivative R2-polar-Cartan))
     v)
    w)
   f)
  R2-rect-point))
#| Result:
0
|#

(pec
 (((((- (covariant-derivative R2-polar-Cartan)
        (covariant-derivative 
         (Cartan-transform R2-polar-Cartan R2-rect-basis)))
     v)
    w)
   f)
  R2-rect-point))
#| Result:
0
|#

(define X (literal-vector-field 'X R2-rect))

(define V (literal-vector-field 'V R2-rect))

(pec (((((covariant-derivative CF) X) V)
       (literal-manifold-function 'F R2-rect))
      R2-rect-point)
    (compose arg-suppressor simplify))
#| Result:
(+ (* G^0_00 V^0 ((partial 0) F) X^0)
   (* G^1_00 V^0 ((partial 1) F) X^0)
   (* G^0_10 ((partial 0) F) V^1 X^0)
   (* G^1_10 ((partial 1) F) V^1 X^0)
   (* G^0_01 V^0 ((partial 0) F) X^1)
   (* G^1_01 V^0 ((partial 1) F) X^1)
   (* G^0_11 ((partial 0) F) V^1 X^1)
   (* G^1_11 ((partial 1) F) V^1 X^1)
   (* ((partial 0) F) ((partial 0) V^0) X^0)
   (* ((partial 0) F) ((partial 1) V^0) X^1)
   (* ((partial 1) F) ((partial 0) V^1) X^0)
   (* ((partial 1) F) ((partial 1) V^1) X^1))
|#

|#

#|
(define-coordinates (up x y) R2-rect)
(define rect-basis (coordinate-system->basis R2-rect))

(define-coordinates (up r theta) R2-polar)
(define polar-basis (coordinate-system->basis R2-polar))

(define rect-chi (R2-rect '->coords))
(define rect-chi-inverse (R2-rect '->point))
(define polar-chi (R2-polar '->coords))
(define polar-chi-inverse (R2-polar '->point))
(define m2 (rect-chi-inverse (up 'x0 'y0)))

(define rect-Christoffel
  (make-Christoffel
   (let ((zero (lambda (m) 0)))
     (down (down (up zero zero)
		 (up zero zero))
	   (down (up zero zero)
		 (up zero zero))))
   rect-basis))

(define polar-Christoffel
  (make-Christoffel
   (let ((zero (lambda (m) 0)))
     (down (down (up zero zero)
		 (up zero (/ 1 r)))
	   (down (up zero (/ 1 r))
		 (up (* -1 r) zero))))
   polar-basis))

(define rect-Cartan
  (Christoffel->Cartan rect-Christoffel))

(define polar-Cartan
  (Christoffel->Cartan polar-Christoffel))

(define J (- (* x d/dy) (* y d/dx)))

(define f (literal-scalar-field 'f R2-rect))

(pec
 (((((covariant-derivative rect-Cartan) 
     d/dx)
    J)
   f)
  m2))
#| Result:
(((partial 1) f) (up x0 y0))
|#

(pec
 (((((covariant-derivative polar-Cartan) 
     d/dx)
    J)
   f)
  m2))
#| Result:
(((partial 1) f) (up x0 y0))
|#
|#

#|
;;; More generally, can show independence here

(define v (literal-vector-field 'v R2-rect))
(define w (literal-vector-field 'w R2-rect))

(pec
 (((((- (covariant-derivative rect-Cartan)
	(covariant-derivative polar-Cartan))
     v)
    w)
   f)
  m2))
#| Result:
0
|#

(define v (literal-vector-field 'v R2-polar))
(define w (literal-vector-field 'w R2-polar))

(pec
 (((((- (covariant-derivative rect-Cartan)
	(covariant-derivative polar-Cartan))
     v)
    w)
   f)
  m2))
#| Result:
0
|#

|#

;;; Over a map

(define (covariant-derivative-over-map Cartan map)
  (assert (Cartan? Cartan))
  (define (nabla X)
    (define (nabla_X V)
      (cond ((function? V) (X V))
	    ((vector-field? V)
	     (((covariant-derivative-over-map-vector Cartan map) X) V))
	    ((form-field? V)
	     (((covariant-derivative-over-map-form Cartan map) X) V))
	    (else
	     (error "Bad input -- covariant-derivative"))))
    (make-operator nabla_X `(nabla ,(diffop-name X))))
  nabla)

(define (((covariant-derivative-over-map-vector Cartan map) V) U-over-map)
  (let ((Cartan-over-map (Cartan->Cartan-over-map Cartan map)))
    (let ((basis (Cartan->basis Cartan-over-map))
	  (Cartan-forms (Cartan->forms Cartan-over-map)))
      (let ((vector-basis (basis->vector-basis basis))
	    (1form-basis (basis->1form-basis basis)))
	(let ((u-components (1form-basis U-over-map)))
	  (let ((deriv-components
		 (+ (V u-components)
		    (* (Cartan-forms ((differential map) V)) u-components))))
	    (define (the-derivative f)
	      (* (vector-basis f) deriv-components))
	    (procedure->vector-field the-derivative
	      `((nabla ,(diffop-name V)) ,(diffop-name U-over-map)))))))))

(define (((covariant-derivative-over-map-form Cartan map) V) tau)
  (let ((k (get-rank tau))
	(nabla_V ((covariant-derivative-over-map-vector Cartan map) V)))
    (procedure->nform-field
     (lambda vectors
       (assert (= k (length vectors)))
       (- (V (apply tau vectors))
	  (sigma (lambda (i)
		   (apply tau
			  (list-with-substituted-coord vectors i
				(nabla_V (list-ref vectors i)))))
		 0 (- k 1))))
     k
     `((nabla ,(diffop-name V)) ,(diffop-name tau)))))

(define (Cartan->Cartan-over-map Cartan map)
  (let ((basis (basis->basis-over-map map (Cartan->basis Cartan)))
	(Cartan-forms
	 (s:map/r (form-field->form-field-over-map map)
		  (Cartan->forms Cartan))))
    (make-Cartan Cartan-forms basis)))

#|
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

(pec
 (s:map/r 
  (lambda (omega)
    ((omega
      (((covariant-derivative-over-map sphere-Cartan gamma:N->M) 
        d/dt) 
       w))
     ((the-real-line '->point) 'tau)))
  (basis->1form-basis basis-over-gamma)))

#| Result:
(up
 (+ (* -1 (w1 tau) ((D beta) tau) (cos (alpha tau)) (sin (alpha tau)))
    ((D w0) tau))
 (/
  (+ (* (w1 tau) (cos (alpha tau)) ((D alpha) tau))
     (* ((D beta) tau) (cos (alpha tau)) (w0 tau))
     (* (sin (alpha tau)) ((D w1) tau)))
  (sin (alpha tau))))
|#

;;; Geodesic equation

(pec
 (s:map/r
  (lambda (omega)
    ((omega
      (((covariant-derivative-over-map sphere-Cartan gamma:N->M)
	d/dt)
       ((differential gamma:N->M) d/dt)))
     ((the-real-line '->point) 't)))
  (basis->1form-basis basis-over-gamma)))

#| Result:
(up
 (+ (* -1 (sin (alpha t)) (expt ((D beta) t) 2) (cos (alpha t)))
    (((expt D 2) alpha) t))
 (+ (/ (* 2 ((D beta) t) (cos (alpha t)) ((D alpha) t)) (sin (alpha t)))
    (((expt D 2) beta) t)))
|#

|#

#|
;;; Geodesic equation

(define-coordinates (up x y) R2-rect)

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


(define CG
  (make-Christoffel G (coordinate-system->basis R2-rect)))

(define gamma:N->M
  (compose (R2-rect '->point)
           (up (literal-function 'alpha)
               (literal-function 'beta))
           (the-real-line '->coords)))

(define basis-over-gamma
  (basis->basis-over-map gamma:N->M
			 (coordinate-system->basis R2-rect)))

(define u
  (basis-components->vector-field
   (up (compose (literal-function 'u0)
                (the-real-line '->coords))
       (compose (literal-function 'u1)
                (the-real-line '->coords)))
   (basis->vector-basis basis-over-gamma)))


(pec
 (s:map/r
  (lambda (omega)
    ((omega
      (((covariant-derivative-over-map (Christoffel->Cartan CG) gamma:N->M)
	d/dt)
       u))
     ((the-real-line '->point) 't)))
  (basis->1form-basis basis-over-gamma)))
#| Result:
(up
 (+ (* ((D beta) t) (u0 t) (G^0_01 (up (alpha t) (beta t))))
    (* ((D beta) t) (u1 t) (G^0_11 (up (alpha t) (beta t))))
    (* ((D alpha) t) (u0 t) (G^0_00 (up (alpha t) (beta t))))
    (* ((D alpha) t) (u1 t) (G^0_10 (up (alpha t) (beta t))))
    ((D u0) t))
 (+ (* ((D beta) t) (u0 t) (G^1_01 (up (alpha t) (beta t))))
    (* ((D beta) t) (u1 t) (G^1_11 (up (alpha t) (beta t))))
    (* ((D alpha) t) (u0 t) (G^1_00 (up (alpha t) (beta t))))
    (* ((D alpha) t) (u1 t) (G^1_10 (up (alpha t) (beta t))))
    ((D u1) t)))
|#

(pec
 (s:map/r
  (lambda (omega)
    ((omega
      (((covariant-derivative-over-map (Christoffel->Cartan CG) gamma:N->M)
	d/dt)
       ((differential gamma:N->M) d/dt)))
     ((the-real-line '->point) 't)))
  (basis->1form-basis basis-over-gamma)))

#| Result:
(up
 (+ (* (expt ((D beta) t) 2) (G^0_11 (up (alpha t) (beta t))))
    (* ((D beta) t) ((D alpha) t) (G^0_01 (up (alpha t) (beta t))))
    (* ((D beta) t) ((D alpha) t) (G^0_10 (up (alpha t) (beta t))))
    (* (expt ((D alpha) t) 2) (G^0_00 (up (alpha t) (beta t))))
    (((expt D 2) alpha) t))
 (+ (* (expt ((D beta) t) 2) (G^1_11 (up (alpha t) (beta t))))
    (* ((D beta) t) ((D alpha) t) (G^1_01 (up (alpha t) (beta t))))
    (* ((D beta) t) ((D alpha) t) (G^1_10 (up (alpha t) (beta t))))
    (* (expt ((D alpha) t) 2) (G^1_00 (up (alpha t) (beta t))))
    (((expt D 2) beta) t)))
|#

|#

#|
;;;; Geodesic Equations = Lagrange Equations 

;;; Here I restrict everything to the unit sphere.
;;; The coordinates on the unit sphere

(define-coordinates t R1-rect)
(define-coordinates (up theta phi) S2-spherical)

(define 2-sphere-basis (coordinate-system->basis S2-spherical))

;;; The Christoffel symbols (for r=1) (p.341 MTW) are:
 
(define G-S2-1
  (make-Christoffel
   (let ((zero  (lambda (point) 0))) 
     (down (down (up zero zero)
		 (up zero (/ 1 (tan theta))))
	   (down (up zero (/ 1 (tan theta)))
		 (up (- (* (sin theta) (cos theta))) zero))))
   2-sphere-basis))

(pec (let ((mu:N->M (compose (S2-spherical '->point)
			     (up (literal-function 'mu-theta)
				 (literal-function 'mu-phi))
			     (R1-rect '->coords)))
	   (Cartan (Christoffel->Cartan G-S2-1)))
       (s:map/r 
	(lambda (w)
	  ((w
	    (((covariant-derivative-over-map Cartan mu:N->M) d/dt)
	     ((differential mu:N->M) d/dt)))
	   ((R1-rect '->point) 'tau)))
	(basis->1form-basis
	 (basis->basis-over-map mu:N->M
				(Cartan->basis Cartan))))))
#| Result:
(up (+ (* -1
	  (expt ((D mu-phi) tau) 2)
	  (cos (mu-theta tau))
	  (sin (mu-theta tau)))
       (((expt D 2) mu-theta) tau))
    (+ (/ (* 2 ((D mu-phi) tau)
	     (cos (mu-theta tau))
	     ((D mu-theta) tau))
	  (sin (mu-theta tau)))
       (((expt D 2) mu-phi) tau)))
|#

;;; We can get the geodesic equations as ordinary Lagrange
;;; equations of a free particle constrained to the surface
;;; of the sphere:

(define ((Lfree m) s)
  (let ((t (time s))
	(q (coordinate s))
	(v (velocity s)))
    (* 1/2 m (square v))))

#|
;;; F is really the embedding map, from the coordinates on the sphere
;;; to the 3-space coordinates in the embedding manifold.

;;; This hides the assumption that the R3 manifold is the same one as
;;; the embedding manifold.

(define F
  (compose (R3-rect '->coords)
	   (S2-spherical '->point)
	   coordinate))

;;; Actually (9 June 2009--GJS) this no longer works, because R3-rect
;;; does not accept an S2-spherical point as in the same manifold.  

;;; Fixed by explicit transfer of a point -- see manifold.scm
|#


(define F
  (compose (R3-rect '->coords)
	   (transfer-point S2-spherical R3-rect)
	   (S2-spherical '->point)
	   coordinate))

(define Lsphere
  (compose (Lfree 1) (F->C F)))

(pec (((Lagrange-equations Lsphere)
       (up (literal-function 'theta)
	   (literal-function 'phi)))
      't))
#| Result:
(down
 (+ (((expt D 2) theta) t)
    (* -1 (cos (theta t)) (sin (theta t)) (expt ((D phi) t) 2)))
 (+ (* (expt (sin (theta t)) 2) (((expt D 2) phi) t))
    (* 2 (cos (theta t)) (sin (theta t)) ((D phi) t) ((D theta) t))))
|#


;;; Note these are DOWN while the geodesic equations are UP.  This is
;;; due to the fact that the geodesic equations are raised by the
;;; metric, which is diagonal, here R=1, and cancels an instance
;;; of(expt (sin theta) 2).

;;; Also see p.345 MTW for computing Christoffel symbols from Lagrange
;;; equations.
|#

#|
;;; Exercise on computation of Christoffel symbols.

(install-coordinates R3-rect (up 'x 'y 'z))
(define R3-rect-point ((R3-rect '->point) (up 'x0 'y0 'z0)))

(install-coordinates R3-cyl (up 'r 'theta 'zeta))
(define R3-cyl-point ((R3-cyl '->point) (up 'r0 'theta0 'z0)))

(define mpr (R3-rect '->coords))

(pec (((* d/dr d/dr) mpr) R3-rect-point))
#| Result:
(up 0 0 0)
|#
;;; So \Gamma^r_{rr} = 0, \Gamma^\theta_{rr} = 0

(pec (((* d/dtheta d/dr) mpr) R3-rect-point))
#| Result:
(up (/ (* -1 y0) (sqrt (+ (expt x0 2) (expt y0 2))))
    (/ x0 (sqrt (+ (expt x0 2) (expt y0 2))))
    0)
|#
;;; by hand = -sint d/dx + cost d/dy = 1/r d/dtheta
;;; Indeed.

(pec (((* d/dtheta d/dr) mpr) R3-cyl-point))
#| Result:
(up (* -1 (sin theta0)) (cos theta0) 0)
|#
;;; So \Gamma^r_{r\theta} = 0, \Gamma^\theta_{r\theta} = 1/r

(pec (((* d/dr d/dtheta) mpr) R3-rect-point))
#| Result:
(up (/ (* -1 y0) (sqrt (+ (expt x0 2) (expt y0 2))))
    (/ x0 (sqrt (+ (expt x0 2) (expt y0 2))))
    0)
|#
;;; by hand = -sint d/dx + cost d/dy = 1/r d/dtheta

(pec (((* d/dr d/dtheta) mpr) R3-cyl-point))
#| Result:
(up (* -1 (sin theta0)) (cos theta0) 0)
|#
;;; So \Gammar_{\theta r} = 0, \Gamma\theta_{\theta r} = 1/r

(pec (((* d/dtheta d/dtheta) mpr) R3-rect-point))
#| Result:
(up (* -1 x0) (* -1 y0) 0)
|#
;;; by hand = -r cost d/dx - r sint d/dy = -r d/dr

(pec (((* d/dtheta d/dtheta) mpr) R3-cyl-point))
#| Result:
(up (* -1 r0 (cos theta0)) (* -1 r0 (sin theta0)) 0)
|#
;;; So \Gammar_{\theta \theta} = -r, \Gamma\theta_{\theta \theta} = 0

;;; These are correct Christoffel symbols...
|#

#|
;;; Computation of Covariant derivatives by difference quotient.
;;; CD below is parallel in definition to the Lie Derivative.
;;; Does not seem to depend on a derivative of basis vectors, in fact
;;; the derivative of the basis vectors is multiplied by zero in the 
;;; product rule output.

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

(define q_0 (up 'q_x 'q_y))

(define m_0
  ((R2-rect '->point) q_0))

(define F (literal-manifold-function 'F R2-rect))


(define (((((CD CF chart) v) u) F) m)

  (define (Sigma state) (ref state 0))
  (define (U state) (ref state 1))
  (define (sigma-u sigma u) (up sigma u))

  (define chi (chart '->coords))
  (define chi^-1 (chart '->point))

  ;; ((gamma m) delta) is the point on gamma advanced by delta.

  (define ((gamma m) delta)
    (chi^-1 (+ (chi m) (* delta ((v chi) m)))))

  (let ((basis (Cartan->basis CF)))
    (let ((vector-basis (basis->vector-basis basis))
	  (1form-basis (basis->1form-basis basis)))
      (let ((u^i (1form-basis u)))
	(let ((initial-state
	       (sigma-u (chi m) (u^i m))))

	  (define (bs state)
	    (let ((sigma (Sigma state)))
	      (let ((m_0 (chi^-1 sigma)))
		(up ((v chi) m_0) 
		    (* -1
		       (((Cartan->forms CF) v) m_0)
		       (u^i m_0))))))

	  (define (vs fs)
	    (* (D fs) bs))

	  ;; First-order approximation to A

	  (define (Au delta)
	    (+ (U initial-state)
	       (* delta ((vs U) initial-state))))

	  (define (g delta)
	    (let ((advanced-m ((gamma m) delta)))
	      (* (- (u^i advanced-m) (Au delta))
		 ((vector-basis F) advanced-m))))

	  ((D g) 0))))))

;;; A bit simpler, but lacking in motivation?

(define (((((CD CF chart) v) u) F) m)

  (define (Sigma state) (ref state 0))
  (define (U state) (ref state 1))
  (define (sigma-u sigma u) (up sigma u))

  (define chi (chart '->coords))
  (define chi^-1 (chart '->point))

  ;; ((gamma m) delta) is the point on gamma advanced by delta.

  (define ((gamma m) delta)
    (chi^-1 (+ (chi m) (* delta ((v chi) m)))))

  (let ((basis (Cartan->basis CF)))
    (let ((vector-basis (basis->vector-basis basis))
	  (1form-basis (basis->1form-basis basis)))
      (let ((u^i (1form-basis u)))
	(let ((initial-state
	       (sigma-u (chi m) (u^i m))))

	  ;; First-order approximation to A

	  (define (Au delta)
	    (- (u^i m)
	       (* delta
		  (((Cartan->forms CF) v) m)
		  (u^i m))))

	  (define (g delta)
	    (let ((advanced-m ((gamma m) delta)))
	      (* (- (u^i advanced-m) (Au delta))
		 ((vector-basis F) advanced-m))))

	  ((D g) 0))))))

(let ((CF (Christoffel->Cartan
	   (make-Christoffel G
			     (coordinate-system->basis R2-rect)))))
  (pe (- (((((CD CF R2-rect) X) Y) F) m_0)
	 (((((covariant-derivative CF) X) Y) F) m_0))))
0

(let ((CF (Christoffel->Cartan
	   (make-Christoffel G
			     (coordinate-system->basis R2-polar)))))
  (pe (- (((((CD CF R2-rect) X) Y) F) m_0)
	 (((((covariant-derivative CF) X) Y) F) m_0))))
0

(let ((CF (Christoffel->Cartan
	   (make-Christoffel G
			     (coordinate-system->basis R2-rect)))))
  (pe (- (((((CD CF R2-polar) X) Y) F) m_0)
	 (((((covariant-derivative CF) X) Y) F) m_0))))
0

(let ((CF (Christoffel->Cartan
	   (make-Christoffel G
			     (coordinate-system->basis R2-polar)))))
  (pe (- (((((CD CF R2-polar) X) Y) F) m_0)
	 (((((covariant-derivative CF) X) Y) F) m_0))))
;Too slow...
|#

#|
;;; Testing on forms.

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
0


(define tau (literal-1form-field 'tau R2-rect))

(define Z (literal-vector-field 'Z R2-rect))

(let* ((CF (Christoffel->Cartan
	    (make-Christoffel G
			      (coordinate-system->basis R2-rect))))
       (D_x ((covariant-derivative CF) X)))
  (pe (- (((D_x (wedge omega tau)) Y Z) m_0)
	 (+ (((wedge omega (D_x tau)) Y Z) m_0)
	    (((wedge (D_x omega) tau) Y Z) m_0)))))
0

(let* ((CF (Christoffel->Cartan
	    (make-Christoffel G
			      (coordinate-system->basis R2-polar))))
       (D_x ((covariant-derivative CF) X)))
  (pe (- (((D_x (wedge omega tau)) Y Z) m_0)
	 (+ (((wedge omega (D_x tau)) Y Z) m_0)
	    (((wedge (D_x omega) tau) Y Z) m_0)))))
0
|#