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

;;; (load "/usr/local/scmutils/src/simplify/test" rule-environment)

(define test
  (rule-system
    ( (atan (? y) (? x))

      (let ((s (simplify `(gcd ,y ,x))))
	(if (equal? s 1)
	    #f
	    (begin (match-assign! 'temp *dictionary* s)
		   #t)))

      (atan (: (simplify `(/ ,y ,temp)))
	    (: (simplify `(/ ,x ,temp)))) )
   ))




#|

(pe ((access test rule-environment)
     '(atan (* (sin (theta_0 t)) (sin (phi_0 t)))
	    (* (cos (phi_0 t)) (sin (theta_0 t))))))
(phi_0 t)
|#