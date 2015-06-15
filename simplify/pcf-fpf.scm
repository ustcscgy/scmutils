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

;;;; This is the top level of sparse polynomial gcd stuff

(declare (usual-integrations))
				 
(define (poly:gcd-sparse u v)
  (cond ((poly:zero? u) v)
	((poly:zero? v) u)
	((poly:one? u) u)
	((poly:one? v) v)
	((base? u)
	 (if (base? v)
	     (base/gcd u v)
	     (base/gcd u (sparse-base-content (poly->sparse v)))))
	((base? v)
	 (base/gcd (sparse-base-content (poly->sparse u)) v))
	(else
	 (let ((arity (gcd-check-same-arity u v))
	       (tt (gcd-target-type u)))
	   (if (fix:< arity *euclid-breakpoint-arity*)
	       (cond ((explicit-pcf? u)
		      (cond ((explicit-pcf? v)
			     (poly/gcd-memoized u v))
			    ((explicit-fpf? v)
			     (poly/gcd-memoized u (fpf->pcf v)))
			    (else (error "What do I do here?"))))
		     ((explicit-fpf? u)
		      (cond ((explicit-pcf? v)
			     (pcf->fpf (poly/gcd-memoized (fpf->pcf u) v)))
			    ((explicit-fpf? v)
			     (pcf->fpf
			      (poly/gcd-memoized (fpf->pcf u)
						 (fpf->pcf v))))
			    (else (error "What do I do here?")))))
	       (let ((su (poly->sparse u))
		     (sv (poly->sparse v)))
		 (if (or (there-exists? su
			   (lambda (term)
			     (let ((c (sparse-coefficient term)))
			       (or (not (real? c)) (inexact? c)))))
			 (there-exists? sv
			   (lambda (term)
			     (let ((c (sparse-coefficient term)))
			       (or (not (real? c)) (inexact? c))))))
		     poly:one
		     (sparse->poly (sparse-gcd su sv) tt))))))))

(define *euclid-breakpoint-arity* 3)

(define poly:gcd poly:gcd-sparse)

(define (gcd-check-same-arity u v)
  (let ((au (poly:arity u)))
    (if (not (fix:= au (poly:arity v)))
	(error "Unequal arities -- poly:gcd-sparse" u v))
    au))
    
(define (gcd-target-type u)
  (cond ((explicit-pcf? u) '*pcf*)
	((explicit-fpf? u) '*fpf*)
	(else
	 (error "Unknown type: gcd-target-type" u))))

(define (poly->sparse p)
  (cond ((explicit-pcf? p) (pcf->sparse p))
	((explicit-fpf? p) (fpf:terms p))
	(else
	 (error "Unknown type: poly->sparse" p))))

(define (sparse->poly s type)
  (cond ((eq? type '*pcf*)
	 (sparse->pcf s))
	((eq? type '*fpf*)
	 (fpf:make s))
	(else
	 (error "Unknown type: sparse->poly" s type))))  

(define (fpf->pcf p)
  (sparse->pcf (fpf:terms v)))

(define (pcf->fpf p)
  (fpf:make (pcf->sparse p)))

(define (pcf->sparse p)
  (let lp ((p p) (arity (poly:arity p)))
    ;;(pp `((p ,p) (arity ,arity)))
    (if (base? p)
	(if (zero? p)
	    '()
	    (list (sparse-term (make-list arity 0) p)))
	(let ((degree (poly:degree p))
	      (c (poly:leading-coefficient p))
	      (r (poly:except-leading-term arity p)))
	  ;;(pp (list degree c r))
	  (sparse-combine-like-terms
	   (append
	    (map (lambda (s-term)
		   (sparse-term (cons degree (sparse-exponents s-term))
				(sparse-coefficient s-term)))
		 (lp c (fix:- arity 1)))
	    (lp r arity)))))))

(define (sparse->pcf s)
  (if (null? s)
      poly:zero
      (let ((v (poly:new-variables (length (sparse-exponents (car s))))))
	(a-reduce poly:+
		  (map (lambda (sterm)
			 (poly:* (sparse-coefficient sterm)
				 (a-reduce poly:*
					   (map poly:expt
						v
						(sparse-exponents sterm)))))
		       s)))))


