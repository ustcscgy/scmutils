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

;;; LEGENDRE.SCM -- the Legendre Polynomials returned as a stream or singly

;;; The following defines a stream whose nth term (0-based) is
;;; P[n]. We use the recurrence relation:
;;;      P[0](x) = 1, P[1](x) = x
;;;   and for n > 1
;;;      P[n](x) = ((2n-1)/n)*x*P[n-1](x) - ((n-1)/n)*P[n-2](x)

(define legendre-polynomials
  (cons-stream poly:one
    (cons-stream poly:identity
      (map-streams (lambda (p1 p2)
		     (let* ((n (+ (poly:degree p1) 1))
			    (a (/ (- (* 2 n) 1) n))
			    (b (/ (- n 1) n)))
		       (poly:- (poly:* (poly:scale a poly:identity) p1)
			       (poly:scale b p2))))
		   (tail legendre-polynomials)
		   legendre-polynomials))))

(define (legendre-polynomial n)
  (stream-ref legendre-polynomials n))
