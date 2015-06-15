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

(declare (usual-integrations))

(define (first-elliptic-integral k)
  (if (flo:< k 1.0)
      (let loop ((a 1.0) (b (flo:sqrt (flo:- 1.0 (flo:* k k)))) (c k))
	(if (flo:< (flo:abs c) (* 2.0 *machine-epsilon*))
	    (flo:/ pi/2 a)
	    (loop (flo:/ (flo:+ a b) 2.0) (flo:sqrt (flo:* a b)) (flo:/ (flo:- a b) 2.0))))
      (error "first-elliptic-integral k>=1" k)))

(define (elliptic-integrals k)
  (let loop ((a 1.0) (b (flo:sqrt (flo:- 1.0 (flo:* k k)))) (c k) (d 0.0) (powers-2 1.0))
    (if (flo:< (flo:abs c) (* 2.0 *machine-epsilon*))
	(let ((first-elliptic-integral (flo:/ pi/2 a)))
	  (cons first-elliptic-integral
		(flo:* first-elliptic-integral (flo:- 1.0 (flo:/ d 2.0)))))
	(loop (flo:/ (flo:+ a b) 2.0) (flo:sqrt (flo:* a b)) (flo:/ (flo:- a b) 2.0) 
	      (flo:+ d (flo:* (flo:* c c) powers-2))
	      (flo:* powers-2 2.0)))))

(define (second-elliptic-integral k)
  (cdr (elliptic-integrals k)))