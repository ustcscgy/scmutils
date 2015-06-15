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

(declare (usual-integrations))

(define (pcf->sparse p)
  (let lp ((p p) (arity (poly/arity p)))
    ;(pp `((p ,p) (arity ,arity)))
    (if (base? p)
	(if (zero? p)
	    '()
	    (list (sparse-term (make-list arity 0) p)))
	(let ((degree (poly/degree p))
	      (c (poly/leading-coefficient p))
	      (r (poly/except-leading-term arity p)))
	  ;(pp (list degree c r))
	  (combine-like-terms
	   (append
	    (map (lambda (s-term)
		   (sparse-term (cons degree (sparse-exponents s-term))
				(sparse-coefficient s-term)))
		 (lp c (- arity 1)))
	    (lp r arity)))))))


(define (sparse->pcf s)
  (if (null? s)
      poly/zero
      (let ((v (poly:new-variables (length (sparse-exponents (car s))))))
	(a-reduce poly/add
		  (map (lambda (sterm)
			 (poly/mul (sparse-coefficient sterm)
			  (a-reduce poly/mul
				    (map poly/expt
					 v
					 (sparse-exponents sterm)))))
		       s)))))
				 
(define (pcf:gcd-sparse u v)
  (cond ((poly/zero? u) v)
	((poly/zero? v) u)
	((poly/one? u) u)
	((poly/one? v) v)
	((base? u)
	 (if (base? v)
	     (base/gcd u v)
	     (base/gcd u (univariate-content (pcf->sparse v)))))
	((base? v)
	 (base/gcd (univariate-content (pcf->sparse u)) v))
	(else
	 (let ((arity (poly/check-same-arity u v))
	       (su (pcf->sparse u))
	       (sv (pcf->sparse v)))
	   (if (or (exists
		    (lambda (term)
		      (let ((c (sparse-coefficient term)))
			(or (not (real? c)) (inexact? c))))
		    su)
		   (exists
		    (lambda (term)
		      (let ((c (sparse-coefficient term)))
			(or (not (real? c)) (inexact? c))))
		    sv))
	       poly/one
	       (sparse->pcf (sparse-gcd su sv)))))))

(set! pcf:gcd pcf:gcd-sparse)
(set! poly:gcd pcf:gcd-sparse)
