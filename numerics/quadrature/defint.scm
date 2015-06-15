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

;;;; Numerical definite integration system interface.
;;;   The integrand is a procedure that computes a function.

;;; Use: (definite-integral <integrand> <from> <to> [compile? #t/#f])
#|
(definite-integral (lambda (x) (* (exp (- x)) (log x))) 0.0 :+infinity)
;Value: -.5772156647120303

(define (foo n)
  (definite-integral (lambda (x) (expt (log (/ 1 x)) n)) 0.0 1.0))

(foo 1)
;Value: .9999999998053075

(foo 2)
;Value: 1.9999999997719113

(foo 3)
;Value: 5.999999999805274

(foo 4)
;Value: 23.999999999815316

(foo 5)
;Value: 119.99999999980271

(foo 6)
;Value: 719.9999999997759
|#

(declare (usual-integrations))

;;; Sets the default for whether the integrand should be compiled.
(define *compile-integrand? #t)

;;; Caches the last compiled function, to avoid recompiling.
(define memoized-compiler
  (linear-memoize-1arg generic-procedure->numerical-procedure 10 weak-find-eq?))

;;; Default error specification
(define *definite-integral-allowable-error* 1.0e-11)

(define (definite-integral f t1 t2 #!optional epsilon compile?)
  (if (default-object? epsilon)
      (set! epsilon *definite-integral-allowable-error*))
  (if (default-object? compile?)
      (set! compile? *compile-integrand?))
  (cond ((and (number? t1) (number? t2) (= t1 t2)) 0)
	((not compile?)
	 ((make-definite-integrator f t1 t2 epsilon)
	  'integral))
	(else
	 ((make-definite-integrator (memoized-compiler f) t1 t2 epsilon)
	  'integral))))

(define (definite-integral-with-tol f x1 x2 tol)
  ((make-definite-integrator (memoized-compiler f) x1 x2 tol) 'integral))

(set! *quadrature-neighborhood-width* #f)
