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


(define (inverting-amplifier rs rf op)
  (cascade (series->2-port (resistor rs))
	   (parallel-parallel op
			      (series->2-port (resistor rf)))))

(define LF357
  (opamp 1e5				;Avol
	 10.0				;Rout
	 1e12				;Rin
	 (* :2pi 60)			;principal-pole
	 ))

#|
(cpp
 (rcf:simplify-and-flatten
  ((voltage-transfer-ratio
    (inverting-amplifier '1000 '10000 LF357))
   's)))

(/ (+ 90826.86366957377 (* 9.082652126165545e-4 s))
   (+ -8705.661007734769 (* 1. s)))
;Unspecified return value
|#