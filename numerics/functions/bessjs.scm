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

;;;; Bessel Functions:

;;; Jn+1(x) = 2n/x Jn(x) - Jn-1(x)   ------   unstable!

;;;   But Jn-1(x) = 2n/x Jn(x) - Jn+1(x) is stable, so use Miller trick.

;;;  (bessjs n x) returns a list
;;;    ( J0(x) ... Jn(x) ) 
;;;   that is good to machine precision for x < 2 and large n.

(define (bessjs nmax x)
  (let lp ((n (round-to-even (+ nmax (sqrt (* bessjs-accur (+ nmax 3)))))) 
	   (ans '())
	   (miller-sum 0.0)
	   (Jn+1 1.0)
	   (Jn 1.0))
    (define (next)
      (- (/ (* 2 n Jn) x) Jn+1))
    (define (ms)
      (if (even? n) (+ Jn miller-sum) miller-sum))
    (if (> Jn bessjs-bigno)
	(begin (set! ans (map (lambda (x) (* x bessjs-bigni)) ans))
	       (set! miller-sum (* miller-sum bessjs-bigni))
	       (set! Jn (* Jn bessjs-bigni))
	       (set! Jn+1 (* Jn+1 bessjs-bigni))))
    (cond ((= n 0)
	   (let ((miller-sum (+ Jn (* 2 miller-sum))))
	     (map (lambda (x) (/ x miller-sum))
		  (cons Jn ans))))
	  ((<= n nmax)
	   (let ((Jn-1 (next)))
	     (lp (- n 1) (cons Jn ans) (ms) Jn Jn-1)))
	  (else
	   (lp (- n 1) ans (ms) Jn (next))))))

(define bessjs-bigno (exact->inexact (expt 2 36)))
(define bessjs-bigni (exact->inexact (expt 2 -36)))
(define bessjs-accur 50)

(define (round-to-even x)
  (let ((xn (inexact->exact (round x))))
    (if (odd? xn)
	(+ xn 1)
	xn)))