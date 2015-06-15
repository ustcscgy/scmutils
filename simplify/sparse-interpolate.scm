#| -*-Scheme-*-

$Id$

Copyright (c) 2002 Massachusetts Institute of Technology

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
|#

;;;;         Sparse Multivariate Polynomial Interpolation 
;;;         a probabilistic method based on Richard Zippel's
;;;          "Interpolating Polynomials From Their Values"
;;;   TR89-963, Department of CS, Cornell University, January 1989
;;;      coded and debugged by Gerald Jay Sussman and Dan Zuras  
;;;                        June 1998                  

;;; This code differs from Zippel's in that it does not use modular
;;; arithmetic or do anything special with the Vandermonde matrices
;;; that arise in the problem.  This makes the idea stand out in stark
;;; contrast without the confusing complications introduced by those
;;; optimizations.

(declare (usual-integrations))

;;; Given a polynomial function f of arity n and maximum degree d, to
;;; find a representation for the terms of the polynomial.

(define (sparse-interpolate f n d)
  (let* ((rargs0 (generate-list (- n 1) interpolate-random))
	 (f1 (lambda (x) (apply f x rargs0)))
	 (p1 (univariate-interpolate f1 d)))
    (let stagelp			;p has k vars interpolated
	((k 1) (p p1) (rargs rargs0))
      (if (= k n)
	  p
	  (let* ((fk
		  (lambda (xk+1)
		    (lambda x1-xk
		      (apply f (append x1-xk (list xk+1) (cdr rargs))))))
		 (xk+1s
		  (generate-list (+ d 1) interpolate-random))
		 (ps
		  (map (lambda (xk+1)
			 (interpolate-skeleton (fk xk+1) p))
		       xk+1s))
		 (css
		  (list-transpose
		   (map (lambda (p)
			  (map sparse-coefficient p))
			ps))))
	    (let ((cps
		   (let clp ((css css))
		     (if (null? css)
			 '()
			 (univariate-interpolate-values xk+1s (car css)
			  (lambda (cp) (cons cp (clp (cdr css))))
			  (lambda () (stagelp k p rargs)))))))
	      (stagelp (+ k 1) (expand-poly p cps) (cdr rargs))))))))

#|
(sparse-interpolate
 (lambda (x y z) (+ (* 3 (square x) (cube y)) (* x y z) (* 4 z) 1))
 3
 4)
;Value: (((2 3 0) . 3) ((1 1 1) . 1) ((0 0 1) . 4) ((0 0 0) . 1))
|#

(define (interpolate-skeleton f skeleton-polynomial)
  (let ((skeleton (map sparse-exponents skeleton-polynomial))
	(nterms (length skeleton-polynomial))
	(arity (length (sparse-exponents (car skeleton-polynomial)))))
    (define (new-args)
      (generate-list nterms
		     (lambda (i)
		       (generate-list arity interpolate-random))))
    (let try-again ((trial-arglists (new-args)))
      (let ((matrix
	     (matrix-by-row-list
	      (map (lambda (argument-list)
		     (map (lambda (exponent-list)
			    (apply * (map expt argument-list exponent-list)))
			  skeleton))
		   trial-arglists)))
	    (values
	     (map (lambda (argl) (apply f argl))
		  trial-arglists)))
	(lu-solve matrix
		  (list->vector values)
		  (lambda (coefficients)
		    (filter (lambda (term)
			      (not (zero? (sparse-coefficient term))))
			    (map (lambda (exponent-list coefficient)
				   (sparse-term exponent-list coefficient))
				 skeleton
				 (vector->list coefficients))))
		  (lambda (ignore) (try-again (new-args))))))))

#|
(interpolate-skeleton
 (lambda (x) (+ (* 3 (expt x 5)) (expt x 2) x 4))
 '(((5) . 1) ((2) . 1) ((1) . 1) ((0) . 1)))
;Value: (((5) . 3) ((2) . 1) ((1) . 1) ((0) . 4))
|#

(define (expand-poly p cps)
  (sort
   (apply append
	  (map (lambda (skel-term cp)
		 (let ((old-exponents (sparse-exponents skel-term)))
		   (map (lambda (coeff-term)
			  (sparse-term
			   (append old-exponents (sparse-exponents coeff-term))
			   (sparse-coefficient coeff-term)))
			cp)))
	       p cps))
   sparse-term->))

#|
(pp (expand-poly '(((5) . 3) ((2) . 1) ((1) . 1) ((0) . 4))
		 '( (((1) . 1) ((0) . 3))
		    (((1) . 1))
		    (((3) . 2) ((0) . 4))
		    (((1) . 2) ((0) . 5)) )))
(((5 1) . 1) ((5 0) . 3) ((1 3) . 2) ((2 1) . 1) ((1 0) . 4) ((0 1) . 2) ((0 0) . 5))
|#

;;; f is a univariate polynomial function.  
;;; d+1 is the number of unknown coefficients.
;;; (usually d is the degree of the polynomial)

(define (univariate-interpolate f d)
  (let* ((xs (generate-list (+ d 1) interpolate-random))
	 (fs (map f xs)))
    (univariate-interpolate-values
     xs
     fs
     (lambda (poly) poly)
     (lambda () (univariate-interpolate f d)))))

(define (univariate-interpolate-values xs fs succeed fail)
  (let ((n (length xs)))
    (assert (= n (length fs)))
    (let* ((exponents (iota n))
	   (matrix
	    (matrix-by-row-list
	     (map (lambda (x)
		    (map (lambda (e) (expt x e))
			 exponents))
		  xs))))
      (lu-solve matrix
		(list->vector fs)
		(lambda (coefficients)
		  (succeed (reverse
			    (filter (lambda (term)
				      (not (zero? (sparse-coefficient term))))
				    (map (lambda (exponent coefficient)
					   (sparse-term (list exponent)
							coefficient))
					 exponents
					 (vector->list coefficients))))))
		(lambda (ignore) (fail))))))


(define *interpolate-size* 1000000)

(define (interpolate-random i)
  (+ (random *interpolate-size*) 1))

#|
(univariate-interpolate
 (lambda (x) (+ (* 3 (expt x 5)) (expt x 2) x 4))
 6)
;Value: (((5) . 3) ((2) . 1) ((1) . 1) ((0) . 4))
|#
