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

;;;; Top-level optimization defaults

(declare (usual-integrations))

(define (minimize f lowx highx)
  (brent-min f lowx highx brent-error))
  
(define brent-error 1.0e-5)


;;; f is a function of the parameters.

(define (multidimensional-minimize f parameters)
  (let ((f (compose f vector->list)))
    (let ((result
	   (nelder-mead f
			(list->vector parameters)
			nelder-start-step
			nelder-epsilon
			nelder-maxiter)))
      (if (eq? 'ok (car result))
	  (vector->list (caadr result))
	  (error "Minimizer did not converge")))))

(define nelder-start-step .01)
(define nelder-epsilon 1.0e-10)
(define nelder-maxiter 1000)

#| ;;; Historical nonsense
(define (multidimensional-minimize f x0 cont)
  ;; cont=(lambda (status minimum-point minimum-value) ...)
  (let* ((bundle?
	  (cond ((vector? x0) #f)
		((list-of-vectors? x0) #t)
		(else
		 (error "Bad initial point -- MINIMIZE"
			x0))))
	 (result
	  (nelder-mead (if bundle?
			   (compose f (bundle-vectors (length x0)))
			   f)
		       (if bundle?
			   (flatten-list-of-vectors x0)
			   x0)
		       nelder-start-step
		       nelder-epsilon
		       nelder-maxiter)))
    (cont (eq? 'OK (car result))
	  (if bundle?
	      ((bundle-vectors (length x0)) (caadr result))
	      (caadr result))
	  (cdadr result))))


(define ((bundle-vectors n) qs)
  (let ((dimension (quotient (vector-length qs) n)))
    (let lp ((i 0) (ans '()))
      (if (fix:= i n)
	  (reverse ans)
	  (lp (fix:+ i 1)
	      (cons (subvector qs
			       (fix:* i dimension)
			       (fix:* (fix:+ i 1) dimension))
		    ans))))))

(define (flatten-list-of-vectors l)
  (list->vector (apply append (map vector->list l))))

  (define (list-of-vectors? l)
    (or (null? l)
	(and (pair? l)
	     (vector? (car l))
	     (list-of-vectors? (cdr l)))))
|#

