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

(define *suppressed-argument-list* '())

(define (suppress-arguments arguments)
  (set! *suppressed-argument-list*
	(cons arguments *suppressed-argument-list*))
  (length *suppressed-argument-list*))

(define (clear-arguments)
  (set! *suppressed-argument-list* '())
  0)


(define (arg-suppressor expression)
  (if (pair? expression)
      (if (member (cdr expression) *suppressed-argument-list*)
	  (arg-suppressor (car expression))
	  (cons (arg-suppressor (car expression))
		(arg-suppressor (cdr expression))))
      expression))