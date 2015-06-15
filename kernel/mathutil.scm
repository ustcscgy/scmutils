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

;;;;    Derived Generic Operators

(declare (usual-integrations))

(define (g:cube x)
  (g:* x x x))


(define (g:log10 x)
  (g:/ (g:log x) (g:log 10)))
(define (g:log2 x)
  (g:/ (g:log x) (g:log 2)))

(define (g:exp10 x)
  (g:expt 10 x))
(define (g:exp2 x)
  (g:expt 2 x))


(define (g:tan x)
  (g:/ (g:sin x) (g:cos x)))
(define (g:sec x)
  (g:/ :one (g:cos x)))
(define (g:csc x)
  (g:/ :one (g:sin x)))


(define (g:tanh x)
  (g:/ (g:sinh x) (g:cosh x)))
(define (g:sech x)
  (g:/ :one (g:cosh x)))
(define (g:csch x)
  (g:/ :one (g:sinh x)))


(define (g:asinh z)
  (g:log (g:+ z (g:sqrt (g:+ :one (g:square z))))))

(define (g:acosh z)
  (g:* :two
       (g:log (g:+ (g:sqrt (g:/ (g:+ z :one) :two))
		   (g:sqrt (g:/ (g:- z :one) :two))))))

(define (g:atanh z)
  (g:/ (g:- (g:log (g:+ :one z)) (g:log (g:- :one z))) :two))


(define (g:arg-shift f . shifts)
  (define (g . xs)
    (g:apply f (map g:+ xs shifts)))
  g)

(define (g:arg-scale f . scales)
  (define (g . xs)
    (g:apply f (map g:* xs scales)))
  g)


(define (g:sigma f low high)
  (if (fix:> low high)
      0
      (let lp ((i (fix:+ low 1)) (sum (f low)))
	(if (fix:> i high)
	    sum
	    (lp (fix:+ i 1) (g:+ sum (f i)))))))

;;; The generalized selector:


(define (g:ref x . selectors)
  (ref-internal x selectors))

(define ((component . selectors) x)
  (ref-internal x selectors))

(define (ref-internal x selectors)
  (if (null? selectors)
      x
      (let ((i (car selectors)) (js (cdr selectors)))
	(cond ((vector? x)
	       (if (or (fix:< i 0)
		       (not (fix:< i (vector-length x))))
		   (error "Bad index -- REF" x i js))		 
	       (ref-internal (vector-ref x i) js))
	      ((matrix? x)
	       (if (null? js)
		   (cond ((column-matrix? x)
			  (matrix-ref x i 0))
			 ((row-matrix? x)
			  (matrix-ref x 0 i))
			 (else
			  (error "Not enuf indices -- REF" x i js)))
		   (ref-internal (matrix-ref x i (car js)) (cdr js))))
	      ((structure? x)
	       (if (or (fix:< i 0)
		       (not (fix:< i (s:length x))))
		   (error "Bad index -- REF" x i js))
	       (ref-internal (s:ref x i) js))
	      ((series? x)
	       (ref-internal (stream-ref (series->stream x) i) js))
	      ((stream-pair? x)
	       (ref-internal (stream-ref x i) js))
	      ((list? x)
	       (ref-internal (list-ref x i) js))
	      ((string? x)
	       (if (or (fix:< i 0)
		       (not (fix:< i (string-length x))))
		   (error "Bad index -- REF" x i js))
	       (if (not (null? js))
		   (error "String has no substructure -- REF" x i js))
	       (string-ref x i))
	      ((procedure? x)
	       (if (operator? x)
		   (make-operator (compose (lambda (y)
					     (ref-internal y selectors))
					   x)
				  `(compose (component ,@selectors)
					    ,(operator-name x))
				  (operator-subtype x))
		   (compose (lambda (y)
			      (ref-internal y selectors))
			    x)))
	      (else
	       (error "Unknown compound -- G:REF" x i))))))


(define (g:size x)
  (cond ((vector? x)      (vector-length x))
	((matrix? x)      (matrix-size x))
	((structure? x)   (s:length x))
	((series? x)      #f)
	((stream-pair? x) #f)
	((list? x)        (length x))
	((string? x)      (string-length x))
	(else
	 (error "Unknown compound -- G:size" x))))

;;; Generic composition duplicates composition in utils

(define (g:compose . fs)
  (define (lp fs)
    (cond ((null? (cdr fs)) (car fs))
	  (else (g:compose-2 (car fs) (lp (cdr fs))))))
  (cond ((null? fs) g:identity)
	((null? (cdr fs)) (car fs))
	(else
	 (g:compose-bin (lp (butlast fs))
			(car (last-pair fs))))))

(define (g:identity x) x)

(define (g:compose-2 f g)
  (cond ((pair? g)
	 (lambda x
	   (g:apply f
		    (map (lambda (gi)
			   (g:apply gi x))
			 g))))
	(else
	 (lambda x
	   (f (g:apply g x))))))

(define (g:compose-bin f g)
  (cond ((pair? g)
	 (let ((a
		(a-reduce joint-arity
			  (map g:arity g))))
	   (cond ((equal? a *at-least-zero*)
		  (lambda x
		    (g:apply f
			   (map
			    (lambda (gi)
			      (g:apply gi x))
			    g))))
		 ((equal? a *exactly-zero*)
		  (lambda ()
		    (g:apply f
			   (map (lambda (gi)
				  (gi))
				g))))
		 ((equal? a *at-least-one*)
		  (lambda (x . y)
		    (g:apply f
			   (map (lambda (gi)
				  (g:apply gi x y))
				g))))
		 ((equal? a *exactly-one*)
		  (lambda (x)
		    (g:apply f
			   (map (lambda (gi)
				  (gi x))
				g))))

		 ((equal? a *at-least-two*)
		  (lambda (x y . z)
		    (g:apply f
			   (map (lambda (gi)
				  (g:apply gi x y z))
				g))))
		 ((equal? a *exactly-two*)
		  (lambda (x y)
		    (g:apply f
			   (map (lambda (gi)
				  (gi x y))
				g))))
		 ((equal? a *at-least-three*)
		  (lambda (u x y . z)
		    (g:apply f
			   (map (lambda (gi)
				  (g:apply gi u x y z))
				g))))
		 ((equal? a *exactly-three*)
		  (lambda (x y z)
		    (g:apply f
			   (map (lambda (gi)
				  (gi x y z))
				g))))
		 ((equal? a *one-or-two*)
		  (lambda (x #!optional y)
		    (if (default-object? y)
			(g:apply f
			       (map (lambda (gi)
				      (gi x))
				    g))
			(g:apply f
			       (map (lambda (gi)
				      (gi x y))
				    g)))))
		 (else
		  (lambda x
		    (g:apply f
			   (map
			    (lambda (gi)
			      (g:apply gi x))
			    g)))))))
	(else
	 (let ((a (g:arity g)))
	   (cond ((equal? a *at-least-zero*)
		  (lambda x
		    (g:apply f
			     (list (g:apply g x)))))
		 ((equal? a *exactly-zero*)
		  (lambda ()
		    (g:apply f
			     (list (g:apply g '())))))
		 ((equal? a *at-least-one*)
		  (lambda (x . y)
		    (g:apply f
			     (list (g:apply g x y)))))
		 ((equal? a *exactly-one*)
		  (lambda (x)
		    (g:apply f
			     (list (g:apply g (list x))))))
		 ((equal? a *at-least-two*)
		  (lambda (x y . z)
		    (g:apply f
			     (list (g:apply g x y z)))))
		 ((equal? a *exactly-two*)
		  (lambda (x y)
		    (g:apply f
			     (list (g:apply g (list x y))))))
		 ((equal? a *at-least-three*)
		  (lambda (u x y . z)
		    (g:apply f
			     (list (g:apply g u x y z)))))
		 ((equal? a *exactly-three*)
		  (lambda (x y z)
		    (g:apply f
			     (list (g:apply g (list x y z))))))
		 ((equal? a *one-or-two*)
		  (lambda (x #!optional y)
		    (if (default-object? y)
			(g:apply f
				 (list (g:apply g (list x))))
			(g:apply f
				 (list (g:apply g (list x y)))))))
		 (else
		  (lambda x
		    (g:apply f
			     (list (g:apply g x))))))))))


