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

(define fconstant
  (sigfun:make (constant 1)
	       (sigfun:make-span -25.6 25.6)))

(plot-trace 1 (magnitude fconstant))

(define tdelta
  (signal->time-function (frequency-function->signal fconstant)))

(plot-trace 2 (magnitude tdelta))

(define (train f n)
  (if (= n 0)
      f
      (let* ((span (sigfun:span f))
	     (period (- (sigfun:max span) (sigfun:min span)))
	     (shift (/ period (expt 2 (+ n 1)))))
	(sigfun:make (+ (sigfun:shift (train f (- n 1)) (- shift))
			(sigfun:shift (train f (- n 1)) (+ shift)))
		     span))))

(plot-trace 3 (train tdelta 3))