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

(define (total-time-derivative f)
  (define (Df-on-path q)
    (D (compose f (path->state q))))
  (abstract-to-state-function Df-on-path))

#|

(print-expression
 ((total-time-derivative 
   (lambda (state)
     (let ((t (state->t state))
	   (q (state->q state)))
       (square q))))
  (->state
   't
   (vector 'x 'y) 
   (vector 'vx 'vy))))
(+ (* 2 vx x) (* 2 vy y))
;No value

(print-expression
 ((total-time-derivative 
   (literal-function 'f (-> (^ Real 5) Real) ))
  (->state
   't
   (vector 'x 'y) 
   (vector 'vx 'vy))))

(+ (* vx (((partial 1) f) (vector t x y vx vy)))
   (* vy (((partial 2) f) (vector t x y vx vy)))
   (((partial 0) f) (vector t x y vx vy)))
;No value


|#


#|

;Value: abstract-to-state-function

(pp abstract-to-state-function)
(named-lambda (abstract-to-state-function f)
  (lambda (state)
    (let ((t (state->t state)) (q (state->q state)) (qdot (state->qdot state)))
      (let ((osc-q (osculating-path t q qdot)))
        ((f osc-q) t)))))
;Unspecified return value

(define ((Lpendulum alpha beta) state)
  (+ (* (square (state->qdot state)) alpha)
     (* beta (cos (state->q state)))))
;Value: Lpendulum

(print-expression
 ((compose (Lpendulum 'alpha 'beta)
	   (path->state (literal-function 'q))) 't))
(+ (* alpha (expt ((D q) t) 2)) (* beta (cos (q t))))
;Unspecified return value


(define (Lq q)
  (compose (Lpendulum 'alpha 'beta)
	   (path->state q)))
;Value: Lq

(print-expression
 ((Lq (literal-function 'q)) 't))
(+ (* alpha (expt ((D q) t) 2)) (* beta (cos (q t))))
;Unspecified return value

(print-expression 
 ((abstract-to-state-function Lq)
  (->state 't 'q 'qdot)))
(+ (* alpha (expt qdot 2)) (* beta (cos q)))
;Unspecified return value

(define (*D* state-function)
  (define (f q) (D (compose state-function (path->state q))))
  (abstract-to-state-function f))
;Value: *D*

(define (F state)
  (let ((x (state->q state)))
    (* x x x)))
;Value: F

(print-expression
 ((*D* F) (->state 't 'x 'xdot)))
(* 3 (expt x 2) xdot)
;Unspecified return value

(print-expression
 ((*D* F) (->state 'a 'b 'c)))
(* 3 (expt b 2) c)
;Unspecified return value

(print-expression
 ((D F) (->state 't 'x 'xdot)))
(vector 0 (* 3 (expt x 2)) 0)
;Unspecified return value

(print-expression
 ((partial_t F) (->state 't 'x 'xdot)))
0
;Unspecified return value


|#