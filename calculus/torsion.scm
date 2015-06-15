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

;;; Experimental definitions of torsion and Riemann...

#|
;;; Old version takes Cartan forms

;;; Hawking and Ellis page 34.

(define ((torsion-vector Cartan) X Y)
  (let ((nabla (covariant-derivative Cartan)))
    (+ ((nabla X) Y)
       (* -1 ((nabla Y) X))
       (* -1 (commutator X Y)))))


;;; Hawking and Ellis equation 2.18, page 35.

(define ((Riemann-curvature Cartan) u v)
  (let ((nabla (covariant-derivative Cartan)))
    (- (commutator (nabla u) (nabla v))
       (nabla (commutator u v)))))
|#

;;; This version takes nabla, allowing use over the map, etc.

(define ((torsion-vector nabla) X Y)
  (+ ((nabla X) Y)
     (* -1 ((nabla Y) X))
     (* -1 (commutator X Y))))


;;; Hawking and Ellis equation 2.18, page 35.

(define ((Riemann-curvature nabla) u v)
  (- (commutator (nabla u) (nabla v))
     (nabla (commutator u v))))


;;; Riemann tensor now takes the nabla too.

(define (Riemann nabla)
  (define (the-Riemann-tensor w x u v)
    (w (((Riemann-curvature nabla) u v) x)))
  (declare-argument-types! the-Riemann-tensor
			   (list 1form-field?
				 vector-field?
				 vector-field?
				 vector-field?))
  the-Riemann-tensor)

(define (torsion nabla)
  (define (the-torsion w x y)
    (w ((torsion-vector nabla) x y)))
  (declare-argument-types! the-torsion
			   (list 1form-field?
				 vector-field?
				 vector-field?))
  the-torsion)


;;; use the connection derived from Lagrange equations on a sphere
;;; compute torsion for the non-symmetrized connection
;;; compute curvature to see if different from the symmetrized connection

;;; find the formula for the connection from the Lagrangian

(define ((Lfree m) s)
  (* 1/2 m (square (velocity s))))


;;; constraint to be on the sphere

(define ((F R) s)
  (let ((q (coordinate s)))
    (let ((theta (ref q 0))
	  (phi (ref q 1)))
      (up (* R (sin theta) (cos phi))
	  (* R (sin theta) (sin phi))
	  (* R (cos theta))))))


;;; Thus the Lagrangian is:

(define (Lsphere m R)
  (compose (Lfree m) (F->C (F R))))

#|

(pe (((Lagrange-equations (Lsphere 'm 'R))
      (up (literal-function 'theta)
	  (literal-function 'phi)))
     't))
(down
 (+ (* -1 (expt R 2) m (cos (theta t)) (sin (theta t)) (expt ((D phi) t) 2))
    (* (expt R 2) m (((expt D 2) theta) t)))
 (+
  (* 2 (expt R 2) m ((D theta) t) (cos (theta t)) (sin (theta t)) ((D phi) t))
  (* (expt R 2) m (((expt D 2) phi) t) (expt (sin (theta t)) 2))))

solving for the highest order terms...

(up (+ (((expt D 2) theta) t)
       (* -1  (cos (theta t)) (sin (theta t)) (expt ((D phi) t) 2)))
    (+ (((expt D 2) phi) t)
       (/ (* 2 (cos (theta t)) ((D theta) t) ((D phi) t)) (sin (theta t)))))
|#

#|

(define-coordinates t the-real-line)

(define M R2-rect)

(define-coordinates (up theta phi) M)

(define M-basis (coordinate-system->basis M))


;;; Unsymmetrized Christoffel symbols. 

(define Gamma-sphere
  (make-Christoffel
   (let ((zero (lambda (point) 0)))
     (down (down (up zero zero)
		 (up zero (/ 2 (tan theta))))
	   (down (up zero zero)
		 (up (- (* (sin theta) (cos theta))) zero))))
   M-basis))

(define sphere-Cartan (Christoffel->Cartan Gamma-sphere))


(define v (literal-vector-field 'v M))
(define w (literal-vector-field 'w M))
(define f (literal-manifold-function 'f M))

(set! *divide-out-terms* #f)


(pec ((((torsion-vector (covariant-derivative sphere-Cartan)) v w) f)
      ((M '->point) (up 'theta 'phi))))
#| Result:
(/ (+ (* 2
	 (cos theta)
	 (((partial 1) f) (up theta phi))
	 (w^1 (up theta phi))
	 (v^0 (up theta phi)))
      (* -2
	 (cos theta)
	 (((partial 1) f) (up theta phi))
	 (w^0 (up theta phi))
	 (v^1 (up theta phi))))
   (sin theta))

;;; Hand simplified
(* 2
   (/ (cos theta) (sin theta))
   (((partial 1) f) (up theta phi))
   (- (* (w^1 (up theta phi)) (v^0 (up theta phi)))
      (* (w^0 (up theta phi)) (v^1 (up theta phi)))))
|#
;;; So the torsion is certainly not zero.



;;; An arbitrary path on the sphere.

(define gamma
  (compose (M '->point)
	   (up (literal-function 'theta)
	       (literal-function 'phi))
	   (the-real-line '->coords)))

(define basis-over-gamma
  (basis->basis-over-map gamma M-basis))

;;; geodesic equations

(pe
 (((((covariant-derivative-over-map sphere-Cartan gamma) d/dt)
    ((differential gamma) d/dt))
   (M '->coords))
  ((the-real-line '->point) 't)))
(up (+ (((expt D 2) theta) t)
       (* -1 (sin (theta t)) (cos (theta t)) (expt ((D phi) t) 2)))
    (+ (((expt D 2) phi) t)
       (/ (* 2 ((D theta) t) (cos (theta t)) ((D phi) t))
	  (sin (theta t)))))

agrees with the Lagrange equations.
|#

#|
torsion for this connection

(define a-function
  (compose (literal-function 'f (-> (UP Real Real) Real))
	   (M '->coords)))

(for-each
 (lambda (x)
   (for-each
    (lambda (y)
      (pe
       ((((torsion-vector (covariant-derivative sphere-Cartan)) x y) a-function)
	((M '->point) (up 'theta 'phi)))))
    (list d/dtheta d/dphi)))
 (list d/dtheta d/dphi))
0
(/ (* 2 (((partial 1) f) (up theta phi)) (cos theta)) (sin theta))
(/ (* -2 (((partial 1) f) (up theta phi)) (cos theta)) (sin theta))
0

nonzero torsion
|#

#|

now compute curvature

(pec (((Riemann (covariant-derivative sphere-Cartan)) dphi d/dtheta d/dphi d/dtheta)
      ((M '->point) (up 'theta 'phi))))
#| Result:
0
|#

;;; But the symmetrized connection is

(define G-S2-1
  (make-Christoffel
   (let ((zero  (lambda (point) 0))) 
     (down (down (up zero zero)
                 (up zero (/ 1 (tan theta))))
           (down (up zero (/ 1 (tan theta)))
                 (up (- (* (sin theta) (cos theta))) zero))))
   M-basis))

(define symmetrized-Cartan (Christoffel->Cartan G-S2-1))

(pec (((Riemann (covariant-derivative symmetrized-Cartan)) dphi d/dtheta d/dphi d/dtheta)
      ((M '->point) (up 'theta 'phi))))
#| Result:
1
|#

Note that the curvature is different for the 
symmetrized (torsion-free) connection.



(for-each
 (lambda (alpha)
   (for-each
    (lambda (beta)
      (for-each
       (lambda (gamma)
	 (for-each
	  (lambda (delta)
	    (newline)
	    (pe `(,alpha ,beta ,gamma ,delta))
	    (pe (((Riemann (covariant-derivative sphere-Cartan)) 
		  alpha beta gamma delta)
		 ((M '->point) (up 'theta 'phi)))))
	  (list d/dtheta d/dphi)))
       (list d/dtheta d/dphi)))
    (list d/dtheta d/dphi)))
 (list dtheta dphi))
;;; nonzero elements are...
(dtheta d/dphi d/dtheta d/dphi)
1

(dtheta d/dphi d/dphi d/dtheta)
-1

;;; NOT SAME AS SYMMETRIC CASE ....

|#

#|

check with usual connection...

(for-each
 (lambda (alpha)
   (for-each
    (lambda (beta)
      (for-each
       (lambda (gamma)
	 (for-each
	  (lambda (delta)
	    (newline)
	    (pe `(,alpha ,beta ,gamma ,delta))
	    (pe (((Riemann (covariant-derivative symmetrized-Cartan)) 
		  alpha beta gamma delta)
		 ((M '->point) (up 'theta 'phi)))))
	  (list d/dtheta d/dphi)))
       (list d/dtheta d/dphi)))
    (list d/dtheta d/dphi)))
 (list dtheta dphi))
;;; nonzero components

(dtheta d/dphi d/dtheta d/dphi)
(expt (sin theta) 2)

(dtheta d/dphi d/dphi d/dtheta)
(* -1 (expt (sin theta) 2))

(dphi d/dtheta d/dtheta d/dphi)
-1

(dphi d/dtheta d/dphi d/dtheta)
1

|#

#|

modified from ricci.scm

(define ((Ricci nabla basis) u v)
  (let ((R (Riemann-curvature nabla)))
    (contract
     (lambda (d/dx dx) (dx ((R d/dx v) u)))
     basis)))

(for-each
 (lambda (alpha)
   (for-each
    (lambda (beta)
      (pe `(,alpha ,beta))
      (pe (((Ricci (covariant-derivative symmetrized-Cartan)
		   M-basis) 
		  alpha beta)
	   ((M '->point) (up 'theta 'phi)))))
    (list d/dtheta d/dphi)))
 (list d/dtheta d/dphi))

(d/dtheta d/dtheta)
1
(d/dtheta d/dphi)
0
(d/dphi d/dtheta)
0
(d/dphi d/dphi)
(expt (sin theta) 2)

(for-each
 (lambda (alpha)
   (for-each
    (lambda (beta)
      (pe `(,alpha ,beta))
      (pe (((Ricci (covariant-derivative sphere-Cartan)
		   M-basis) 
		  alpha beta)
	   ((M '->point) (up 'theta 'phi)))))
    (list d/dtheta d/dphi)))
 (list d/dtheta d/dphi))
(d/dtheta d/dtheta)
0
(d/dtheta d/dphi)
0
(d/dphi d/dtheta)
0
(d/dphi d/dphi)
1

;;; Ricci curvatures are different

;;; to get the scalar curvature we need to raise one index 
;;; and then contract =>  so need a metric
|#

;;; investigate relation between geodesic deviation
;;; and variational evolution equations

#|

(define (F-sphere-evolution state)
  (let ((theta (ref state 0))
	(phi (ref state 1))
	(theta-dot (ref state 2))
	(phi-dot (ref state 3)))
    (up theta-dot
	phi-dot
	(* (cos theta) (sin theta) (expt phi-dot 2))
	(* -1 (/ (* 2 (cos theta) theta-dot phi-dot) (sin theta))))))

(pe (* ((D F-sphere-evolution) 
	(up 'theta 'phi 'theta-dot 'phi-dot))
       (up 'dtheta 'dphi 'dtheta-dot 'dphi-dot)))
(up
 dtheta-dot
 dphi-dot
 (+ (* 2 dtheta (expt phi-dot 2) (expt (cos theta) 2))
    (* 2 dphi-dot phi-dot (cos theta) (sin theta))
    (* -1 dtheta (expt phi-dot 2)))
 (/
  (+ (* -2 dphi-dot theta-dot (cos theta) (sin theta))
     (* -2 dtheta-dot phi-dot (cos theta) (sin theta))
     (* 2 dtheta phi-dot theta-dot))
  (expt (sin theta) 2)))

(define (variational-sphere-evolution vstate)
  (let ((theta (ref state 0))
	(phi (ref state 1))
	(theta-dot (ref state 2))
	(phi-dot (ref state 3))
	(dtheta (ref state 4))
	(dphi (ref state 5))
	(dtheta-dot (ref state 6))
	(dphi-dot (ref state 7)))
    (up theta-dot
	phi-dot
	(* (cos theta) (sin theta) (expt phi-dot 2))
	(* -1 (/ (* 2 (cos theta) theta-dot phi-dot) (sin theta)))
	dtheta-dot
	dphi-dot
	(+ (* 2 dtheta (expt phi-dot 2) (expt (cos theta) 2))
	   (* 2 dphi-dot phi-dot (cos theta) (sin theta))
	   (* -1 dtheta (expt phi-dot 2)))
	(/
	 (+ (* -2 dphi-dot theta-dot (cos theta) (sin theta))
	    (* -2 dtheta-dot phi-dot (cos theta) (sin theta))
	    (* 2 dtheta phi-dot theta-dot))
	 (expt (sin theta) 2)))))

(pe (((((Riemann-curvature
	 (covariant-derivative-over-map sphere-Cartan gamma))
	d/dt d/dt)
       ((differential gamma) d/dt))
      (literal-manifold-function 'f M))
     ((the-real-line '->point) 't)))
0

(define ((geodesic-deviation nabla) gamma n w)
  (let ((u ((differential gamma) d/dt)))
    (- (w ((nabla d/dt) ((nabla d/dt) n)))
       ((Riemann nabla) w u n u))))

(define (n basis gamma)
  (procedure->vector-field
   (lambda (f)
     (lambda (mt)
       (* (((basis->vector-basis basis) f) (gamma mt))
	  (up ((literal-function 'n0) ((the-real-line '->coords) mt))
	      ((literal-function 'n1) ((the-real-line '->coords) mt))))))))

 
(pe (((n M-basis gamma) (literal-manifold-function 'f M))
     ((the-real-line '->point) 't)))
(+ (* (n0 t) (((partial 0) f) (up (theta t) (phi t))))
   (* (n1 t) (((partial 1) f) (up (theta t) (phi t)))))
  
(pe (((geodesic-deviation
       (covariant-derivative-over-map sphere-Cartan gamma))
      gamma
      (n M-basis gamma) 
      ((form-field->form-field-over-map gamma)
       (literal-1form-field 'w M)))
     ((the-real-line '->point) 't)))

|#
