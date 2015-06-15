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

;;; find eigenvalues and eigenvectors for 2d map fixed points

;;; lets assume the form of the map is the same as in explore map
;;; (map x y cont fail)  (cont xp yp)

#|
(define (fixed-point-eigen T xe ye eps cont)
  (let ((M00 ((richardson-derivative 
	       (lambda (dx)
		 (T (+ xe dx) ye 
		    (lambda (x y) ((principal-value pi) (- x xe)))
		    'failure))
	       eps) 0.0))
	(M01 ((richardson-derivative 
	       (lambda (dx)
		 (T xe (+ ye dx) 
		    (lambda (x y) ((principal-value pi) (- x xe)))
		    'failure))
	       eps) 0.0))
	(M10 ((richardson-derivative 
	       (lambda (dx)
		 (T (+ xe dx) ye (lambda (x y) y) 'failure))
	       eps) 0.0))
	(M11 ((richardson-derivative 
	       (lambda (dx)
		 (T xe (+ ye dx) (lambda (x y) y) 'failure))
	       eps) 0.0)))
    (let ((trace (+ M00 M11))
	  (determinant (- (* M00 M11) (* M01 M10))))
      (quadratic 
       1. (- trace) determinant 
       (lambda (root1 root2)
	 (cont root1 M01 (- root1 M00)
	       root2 M01 (- root2 M00)))))))

|#
#|

(let* ((g 9.8)
       (l 1.0)
       (m 1.0)
       (A 0.01)
       (ml^2 (* m l l))
       (mlg (* m l g))
       (omega (* 4.2 (sqrt (/ g l))))
       (Aw^2/g (* A omega omega (/ 1 g))))
  (fixed-point-eigen
   (dp-map ml^2 mlg Aw^2/g omega)
   0. 0. 1.e-12 list))
;Value 24: 
(.07387928185655412+.9972671917356725i 
 .31028395146610266 6.245004513516506e-16+.9972671917356725i
 .07387928185655411-.9972671917356724i
 .31028395146610266 6.106226635438361e-16-.9972671917356724i)
;Value 16: (.07387928185655412+.9972671917356725i .07387928185655411-.9972671917356724i)			   
(magnitude .07387928185655412+.9972671917356725i)
;Value: .9999999999999974
(let ((omega (* 4.2 (sqrt 9.8))))
  (* (/ (angle .07387928185655412+.9972671917356725i)
	2pi)
     omega))
;Value: 3.132280497149872

(let* ((g 9.8)
       (l 1.0)
       (m 1.0)
       (A 0.00)
       (ml^2 (* m l l))
       (mlg (* m l g))
       (omega (* 4.2 (sqrt (/ g l))))
       (Aw^2/g (* A omega omega (/ 1 g))))
  (fixed-point-eigen
   (dp-map ml^2 mlg Aw^2/g omega)
   0. 0. 1.e-12 list))
;Value 25: (.07473009358642436+.997203797181181i .3185450682740048 4.440892098500626e-16+.997203797181181i .07473009358642436-.997203797181181i .3185450682740048 4.440892098500626e-16-.997203797181181i)
(let ((omega (* 4.2 (sqrt 9.8))))
  (* (/ (angle .07473009358642436+.997203797181181i)
	2pi)
     omega))
;Value: 3.1304951684997055
(sqrt 9.8)
;Value: 3.1304951684997055

(let* ((g 9.8)
       (l 1.0)
       (m 1.0)
       (A 0.01)
       (ml^2 (* m l l))
       (mlg (* m l g))
       (omega (* 4.2 (sqrt (/ g l))))
       (Aw^2/g (* A omega omega (/ 1 g))))
  (fixed-point-eigen
   (dp-map ml^2 mlg Aw^2/g omega)
   pi 0. 1.e-12 list))
;Value 26: 
(4.461381889207846 .6881802386332591 2.1186180460230517
 .22414579716196714 .6881802386332591 -2.1186180460228274)

(let* ((g 9.8)
       (l 1.0)
       (m 1.0)
       (A 0.00)
       (ml^2 (* m l l))
       (mlg (* m l g))
       (omega (* 4.2 (sqrt (/ g l))))
       (Aw^2/g (* A omega omega (/ 1 g))))
  (fixed-point-eigen
   (dp-map ml^2 mlg Aw^2/g omega)
   pi 0. 1.e-12 list))

;Value 27: 
(4.463782504176754 .6771703874392959 2.1198786261298044 
.2240252519167095 .6771703874392959 -2.11987862613024)
(* 4.463782504176754 .2240252519167095)
;Value: .9999999999995977
		     
|#

(define ((unstable-manifold T xe ye dx dy lam eps) param)
  (let ((n (floor->exact (/ (log (/ param eps)) (log lam)))))
    ((iterated-map T n) (+ xe (* dx (/ param (expt lam n))))
			(+ ye (* dy (/ param (expt lam n))))
			cons
			list)))

#| 

;;; faster way to compute homoclinic5

(set! *ode-integration-method* 'bulirsch-stoer)

(define win (frame -pi pi -10 10 1000 1000))
(graphics-clear win)
(graphics-close win)

(let* ((g 9.8)
       (l 1.0)
       (m 1.0)
       (A 0.05)
       (ml^2 (* m l l))
       (mlg (* m l g))
       (omega (* 4.2 (sqrt (/ g l))))
       (Aw^2/g (* A omega omega (/ 1 g)))
       (T (dp-map ml^2 mlg Aw^2/g omega))
       (T^-1 (dp-map ml^2 mlg Aw^2/g (- omega))))
  (fixed-point-eigen T pi 0. 1.e-12 
   (lambda (root1 dx1 dy1 root2 dx2 dy2)
     (plot-parametric-fill
      win
      (unstable-manifold T pi 0. dx1 dy1 root1 1.e-10)
      0.01 500. (cylinder-near? 0.01)
      )
     (plot-parametric-fill
      win
      (unstable-manifold T^-1 pi 0. dx2 dy2 root1 1.e-10)
      0.01 500. (cylinder-near? 0.01)
      )
     (plot-parametric-fill
      win
      (unstable-manifold T pi 0. (- dx1) (- dy1) root1 1.e-10)
      0.01 500. (cylinder-near? 0.01)
      )
     (plot-parametric-fill
      win
      (unstable-manifold T^-1 pi 0. (- dx2) (- dy2) root1 1.e-12)
      .01 500. (cylinder-near? 0.01)
      )
     )))



(define win (frame -pi pi -10 10 1000 1000))
(graphics-clear win)
(graphics-close win)

;;; alternate way
;;; homoclinic5

(let* ((g 9.8)
       (l 1.0)
       (m 1.0)
       (A 0.05)
       (ml^2 (* m l l))
       (mlg (* m l g))
       (omega (* 4.2 (sqrt (/ g l))))
       (Aw^2/g (* A omega omega (/ 1 g)))
       (T (dp-map ml^2 mlg Aw^2/g omega))
       (T^-1 (dp-map ml^2 mlg Aw^2/g (- omega))))
  (fixed-point-eigen T pi 0. 1.e-12 
   (lambda (root1 dx1 dy1 root2 dx2 dy2)
     (plot-parametric
      win
      (unstable-manifold T pi 0. dx1 dy1 root1 1.e-10)
      0.01 500. 0.01
      )
     (plot-parametric
      win
      (unstable-manifold T^-1 pi 0. dx2 dy2 root1 1.e-10)
      0.01 500. 0.01
      )
     (plot-parametric
      win
      (unstable-manifold T pi 0. (- dx1) (- dy1) root1 1.e-10)
      0.01 500. .01
      )
     (plot-parametric
      win
      (unstable-manifold T^-1 pi 0. (- dx2) (- dy2) root1 1.e-12)
      .01 500. .01
      )
     )))

;;; homoclinic5-section

(let* ((g 9.8)
       (l 1.0)
       (m 1.0)
       (A 0.05)
       (ml^2 (* m l l))
       (mlg (* m l g))
       (omega (* 4.2 (sqrt (/ g l))))
       (Aw^2/g (* A omega omega (/ 1 g)))
       (T (dp-map ml^2 mlg Aw^2/g omega)))
  (explore-map win T 1000))

(graphics-clear win)


;;; homoclinic02

(let* ((g 9.8)
       (l 1.0)
       (m 1.0)
       (A 0.002)
       (ml^2 (* m l l))
       (mlg (* m l g))
       (omega (* 4.2 (sqrt (/ g l))))
       (Aw^2/g (* A omega omega (/ 1 g)))
       (T (dp-map ml^2 mlg Aw^2/g omega))
       (T^-1 (dp-map ml^2 mlg Aw^2/g (- omega))))
  (fixed-point-eigen T pi 0. 1.e-12 
   (lambda (root1 dx1 dy1 root2 dx2 dy2)
     (plot-parametric
      win
      (unstable-manifold T pi 0. dx1 dy1 root1 1.e-10)
      0.01 500. 0.01
      )
     (plot-parametric
      win
      (unstable-manifold T^-1 pi 0. dx1 dy1 root1 1.e-10)
      0.01 500. 0.01
      )
     (plot-parametric
      win
      (unstable-manifold T pi 0. (- dx1) (- dy1) root1 1.e-10)
      0.01 500. .01
      )
     (plot-parametric
      win
      (unstable-manifold T^-1 pi 0. (- dx1) (- dy1) root1 1.e-12)
      .01 500. .01
      )
     )))

;;; homoclinic02-section

(let* ((g 9.8)
       (l 1.0)
       (m 1.0)
       (A 0.002)
       (ml^2 (* m l l))
       (mlg (* m l g))
       (omega (* 4.2 (sqrt (/ g l))))
       (Aw^2/g (* A omega omega (/ 1 g)))
       (T (dp-map ml^2 mlg Aw^2/g omega)))
  (explore-map win T 1000))


;;; how about the standard map

(define (fixed-point-eigen T xe ye eps cont)
  (let ((M00 ((richardson-derivative 
	       (lambda (dx)
		 (T (+ xe dx) ye 
		    (lambda (x y) ((principal-value pi) (- x xe)))
		    'failure))
	       eps) 0.0))
	(M01 ((richardson-derivative 
	       (lambda (dx)
		 (T xe (+ ye dx) 
		    (lambda (x y) ((principal-value pi) (- x xe)))
		    'failure))
	       eps) 0.0))
	(M10 ((richardson-derivative 
	       (lambda (dx)
		 (T (+ xe dx) ye 
		    (lambda (x y) ((principal-value pi) (- y ye))) 
		    'failure))
	       eps) 0.0))
	(M11 ((richardson-derivative 
	       (lambda (dx)
		 (T xe (+ ye dx) 
		    (lambda (x y) ((principal-value pi) (- y ye)))
		    'failure))
	       eps) 0.0)))
    (let ((trace (+ M00 M11))
	  (determinant (- (* M00 M11) (* M01 M10))))
      (quadratic 
       1. (- trace) determinant 
       (lambda (root1 root2)
	 (cont root1 M01 (- root1 M00)
	       root2 M01 (- root2 M00)))))))

I' = I + K sin th
th' th + I'

th = th' - I'
I = I' - K sin th

(named-lambda (standard-map K)
  (lambda (x y continue fail)
    (let ((yp (flo:pv (flo:+ y (flo:* K (flo:sin x))))))
      (continue (flo:pv (flo:+ x yp)) yp))))

(define (standard-map-inverse K) 
  (lambda (x y continue fail)
    (let ((thp ((principal-value 2pi) (- x y))))
      (continue thp ((principal-value 2pi) (- y (* K (sin thp))))))))
	 

(define win (frame 0. 2pi 0. 2pi))
(graphics-clear win)

|#

(define (do-it win K len eps1 eps2)
  (graphics-clear win)
  (let ((T (standard-map K))
	(T^-1 (standard-map-inverse K)))
    (fixed-point-eigen T 0. 0. eps1 
		       (lambda (root1 dx1 dy1 root2 dx2 dy2)
			 (plot-parametric-fill
			  win
			  (unstable-manifold T 0. 0. dx1 dy1 root1 eps2)
			  0.01 len (torus-near? 0.01)
			  )
			 (plot-parametric-fill
			  win
			  (unstable-manifold T^-1 0. 0. dx2 dy2 root1 eps2)
			  0.01 len (torus-near? 0.01)
			  )
			 (plot-parametric-fill
			  win
			  (unstable-manifold T 0. 0. (- dx1) (- dy1) root1 eps2)
			  0.01 len (torus-near? 0.01)
			  )
			 (plot-parametric-fill
			  win
			  (unstable-manifold T^-1 0. 0. (- dx2) (- dy2) root1 eps2)
			  .01 len (torus-near? 0.01)
			  )
			 ))))

#|

(do-it win 0.8 100. 1.e-10 1.e-10)
(do-it win 1. 100. 1.e-10 1.e-10)
(graphics-clear win)
(do-it win 2. 100. 1.e-10 1.e-10)
 (do-it win 2. 1000. 1.e-10 1.e-10) (do-it win 2. 1000. 1.e-10 .01) 
(do-it win 2. 100. 1.e-10 1.e-4)
(do-it win 2. 100. 1.e-10 .01)

  (do-it win 2. 1000. 1.e-10 1.e-10)
  (do-it win 2. 10000. 1.e-10 1.e-10)

(do-it win 2.1 100. 1.e-10 1.e-10)
(do-it win 3. 100. 1.e-10 1.e-10)
(do-it win 5. 100. 1.e-10 1.e-10)

|#