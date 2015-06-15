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

;;; We define the standard singular-matrix failure continuations as follows.

(define (barf-on-zero-pivot dismiss)
  (singular-matrix-error))

(define (allow-zero-pivot dismiss)
  (dismiss))

;;; Rebind this to catch errors
(define singular-matrix-error)

;;; default value
(define (default-singular-matrix-error)
  (error "Singular matrix - zero pivot"))

(set! singular-matrix-error
      default-singular-matrix-error)

(define (with-singular-matrix-handler handler thunk)
  (fluid-let ((singular-matrix-error handler))
    (thunk)))

(define (handle-singularity-errors-with error-handler)
  (set! singular-matrix-error error-handler))
