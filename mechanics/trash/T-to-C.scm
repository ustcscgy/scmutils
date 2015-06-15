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





(define ((T->C T) state)
  (let ((t (state->t state))
	(q (state->q state))
	(qdot (state->qdot state)))
    (->state t
	     (T t q)
	     (* ((derivative T) t q)
		(vector 1 qdot)))))


#|  tricky way to do it if T takes states for input

(define ((T->C T) state)
  (->state (state->t state)
	   (T state)
	   ((D* T) state)))

|#

#|
(define (p->r polar-vector)
  (let ((r (vector-ref polar-vector 0))
        (phi (vector-ref polar-vector 1)))
    (let ((x (* r (cos phi)))
          (y (* r (sin phi))))
      (vector x y))))

(print-expression
 ((T->C (lambda (t q) (p->r q)))
  (->state 't (vector 'r 'th) (vector 'rdot 'thdot))) )

(vector t
	(* r (cos th))
	(* r (sin th))
	(+ (* -1 r thdot (sin th)) (* rdot (cos th)))
	(+ (* r thdot (cos th)) (* rdot (sin th))))
|#

(define ((T->C T) state)
  (let ((t (state->t state))
	(q (state->q state))
	(qdot (state->qdot state)))
    (let ((qosc (osculating-path t q qdot)))
      (let ((qp (lambda (t) (T t (qosc t)))))
	((path->state qp) t)))))

#| with same result

(print-expression
 ((T->C (lambda (t q) (p->r q)))
  (->state 't (vector 'r 'th) (vector 'rdot 'thdot))) )
(vector t 
	(* r (cos th)) 
	(* r (sin th))
	(+ (* -1 r thdot (sin th)) (* rdot (cos th)))
	(+ (* r thdot (cos th)) (* rdot (sin th))))

|#

#|
another idea

(print-expression
 ((path->state 
   (compose p->r 
	    (vector (literal-function 'r)
		    (literal-function 'th))))
  't))
(vector t
        (* (r t) (cos (th t)))
        (* (sin (th t)) (r t))
        (+ (* ((D r) t) (cos (th t))) (* -1 ((D th) t) (sin (th t)) (r t)))
        (+ (* ((D r) t) (sin (th t))) (* ((D th) t) (r t) (cos (th t)))))

(define ((abstract-to-state-function f) state)
  (let ((t (state->t state))
        (q (state->q state))
        (qdot (state->qdot state)))
    (let ((osc-q (osculating-path t q qdot)))
      ((f osc-q) t))))

;;; or should it be  (to keep higher derivatives if needed)
#|
;;; Not working
(define ((abstract-to-state-function f) state)
  (let ((osc-q (apply osculating-path (vector->list state))))
    ((f osc-q) (state->t state))))
|#

(print-expression
 ((abstract-to-state-function
   (lambda (polar-path) 
     (path->state (compose p->r polar-path))))
  (->state 't (vector 'r 'th) (vector 'rdot 'thdot))))
(vector t
        (* r (cos th))
        (* r (sin th))
        (+ (* -1 r thdot (sin th)) (* rdot (cos th)))
        (+ (* r thdot (cos th)) (* rdot (sin th))))

(define (T->C T)
  (abstract-to-state-function
   (lambda (q)
     (path->state (lambda (t) (T t (q t)))))))

(print-expression
 ((T->C (lambda (t q) (p->r q)))
  (->state 't (vector 'r 'th) (vector 'rdot 'thdot))) )

(vector t
        (* r (cos th))
        (* r (sin th))
        (+ (* -1 r thdot (sin th)) (* rdot (cos th)))
        (+ (* r thdot (cos th)) (* rdot (sin th))))

|#


#|
(define ((central-Lagrangian-rectangular m V) state)
  (let ((q (state->q state))
	(qdot (state->qdot state)))
    (- (* 1/2 m (square qdot))
       (V (sqrt (square q))))))

(define ((central-Lagrangian-polar m V) state)
  (let ((q (state->q state))
        (qdot (state->qdot state)))
    (let ((r (vector-ref q 0))
          (phi (vector-ref q 1))
          (rdot (vector-ref qdot 0))
          (phidot (vector-ref qdot 1)))
      (- (* 1/2 m
           (+ (square rdot)
              (square (* r phidot))) )
         (V r)))))
 
;;; or

(print-expression
 ((compose (central-Lagrangian-rectangular 'm (literal-function 'V))
	   (T->C (lambda (t q) (p->r q))))
  (->state 't (vector 'r 'th) (vector 'rdot 'thdot))))
(+ (* 1/2 m (expt r 2) (expt thdot 2))
   (* 1/2 m (expt rdot 2))
   (* -1 (V r)))

;;; same as ....
(print-expression
 ((central-Lagrangian-polar 'm (literal-function 'V))
  (->state 't (vector 'r 'th) (vector 'rdot 'thdot))))
(+ (* 1/2 m (expt r 2) (expt thdot 2))
   (* 1/2 m (expt rdot 2))
   (* -1 (V r)))

|#
#|

(print-expression
 (((lagrange-equations
    (compose (central-lagrangian-rectangular 'm (literal-function 'v))
	   (t->c (lambda (t q) (p->r q)))))
   (vector (literal-function 'r) (literal-function 'th)))
  't))
(vector
 (+ (* m (((expt D 2) r) t))
    (* -1 m (r t) (expt ((D th) t) 2))
    ((D V) (r t)))
 (+ (* m (((expt D 2) th) t) (expt (r t) 2))
    (* 2 m ((D r) t) (r t) ((D th) t))))

(print-expression
 (((Lagrange-equations
    (central-Lagrangian-polar 'm (literal-function 'V)))
   (vector (literal-function 'r) (literal-function 'th)))
  't))
(vector
 (+ (* m (((expt D 2) r) t))
    (* -1 m (r t) (expt ((D th) t) 2))
    ((D V) (r t)))
 (+ (* 2 m ((D r) t) (r t) ((D th) t))
    (* m (((expt D 2) th) t) (expt (r t) 2))))

|#

#| driven pendulum as coordinate transformation

(define ((dpend-coordinates l y_s) t th)
  (let ((x (* l (sin th)))
	(y (- (y_s t) (* l (cos th)))))
    (vector x y)))

(define ((L-uniform-gravity m g) state)
  (let ((q (state->q state))
	(v (state->qdot state)))
    (let ((height (vector-ref q 1)))
      (- (* 1/2 m (square v))
	 (* m g height)))))

(define (L-dpend m g l y_s)
  (compose (L-uniform-gravity m g)
	   (T->C (dpend-coordinates l y_s))))

(print-expression
 ((L-dpend 'm 'g 'l (literal-function 'y_s))
  (->state 't 'th 'thdot)))

(+ (* g l m (cos th))
   (* -1 g m (y_s t))
   (* 1/2 (expt l 2) m (expt thdot 2))
   (* l m thdot ((D y_s) t) (sin th))
   (* 1/2 m (expt ((D y_s) t) 2)))

;;; Lagrange equations

(print-expression
 (((Lagrange-equations
    (L-dpend 'm 'g 'l (literal-function 'y_s)))
   (literal-function 'th)) 
  't))

(+ (* g l m (sin (th t)))
   (* (expt l 2) m (((expt D 2) th) t))
   (* l m (((expt D 2) y_s) t) (sin (th t))))

;;; without the (t)

(print-expression
 ((abstract-to-state-function
   (Lagrange-equations
    (L-dpend 'm 'g 'l (literal-function 'y_s))))
  (vector 't 'theta 'thetadot 'thetadotdot)))


|#

#| double pendulum as coordinate transformation

(define ((double-pend-coordinates l1 l2) t q)
  (let ((th1 (vector-ref q 0))
	(th2 (vector-ref q 1)))
    (let ((x1 (* l1 (sin th1)))
	  (y1 (* -1 l1 (cos th1))))
      (let ((x2 (+ x1 (* l2 (sin th2))))
	    (y2 (+ y1 (* -1 l2 (cos th2)))))
	(vector x1 y1 x2 y2)))))

(define ((L-uniform-gravity-2 m1 m2 g) state)
  (let ((q (state->q state))
	(v (state->qdot state)))
    (let ((h1 (vector-ref q 1))
	  (h2 (vector-ref q 3))
	  (v1 (vector (vector-ref v 0) (vector-ref v 1)))
	  (v2 (vector (vector-ref v 2) (vector-ref v 3))))
      (let ((T (+ (* 1/2 m1 (square v1)) 
		  (* 1/2 m2 (square v2))))
	    (V (+ (* m1 g h1)
		  (* m2 g h2))))
	(- T V)))))

(define (L-double m1 m2 l1 l2 g)
  (compose (L-uniform-gravity-2 m1 m2 g)
	   (T->C (double-pend-coordinates l1 l2))))

(print-expression
 ((L-double 'm1 'm2 'l1 'l2 'g)
  (->state 't (vector 'th1 'th2) (vector 'th1dot 'th2dot))))

(+ (* g l1 m1 (cos th1))
   (* g l1 m2 (cos th1))
   (* g l2 m2 (cos th2))
   (* 1/2 (expt l1 2) m1 (expt th1dot 2))
   (* 1/2 (expt l1 2) m2 (expt th1dot 2))
   (* l1 l2 m2 th1dot th2dot (sin th2) (sin th1))
   (* l1 l2 m2 th1dot th2dot (cos th2) (cos th1))
   (* 1/2 (expt l2 2) m2 (expt th2dot 2)))

(print-expression
 (((Lagrange-equations
    (L-double 'm1 'm2 'l1 'l2 'g))
   (vector (literal-function 'th1) (literal-function 'th2))) 
  't))
(vector
 (+ (* g l1 m1 (sin (th1 t)))
    (* g l1 m2 (sin (th1 t)))
    (* (expt l1 2) m1 (((expt D 2) th1) t))
    (* (expt l1 2) m2 (((expt D 2) th1) t))
    (* l1 l2 m2 (cos (th2 t)) (((expt D 2) th2) t) (cos (th1 t)))
    (* l1 l2 m2 (cos (th2 t)) (expt ((D th2) t) 2) (sin (th1 t)))
    (* l1 l2 m2 (sin (th2 t)) (((expt D 2) th2) t) (sin (th1 t)))
    (* -1 l1 l2 m2 (sin (th2 t)) (expt ((D th2) t) 2) (cos (th1 t))))
 (+ (* g l2 m2 (sin (th2 t)))
    (* -1 l1 l2 m2 (expt ((D th1) t) 2) (cos (th2 t)) (sin (th1 t)))
    (* l1 l2 m2 (expt ((D th1) t) 2) (sin (th2 t)) (cos (th1 t)))
    (* l1 l2 m2 (cos (th2 t)) (cos (th1 t)) (((expt D 2) th1) t))
    (* l1 l2 m2 (sin (th2 t)) (((expt D 2) th1) t) (sin (th1 t)))
    (* (expt l2 2) m2 (((expt D 2) th2) t))))

(print-expression
 ((Lagrange-explicit
   (L-double 'm1 'm2 'l1 'l2 'g))
  (->state 't (vector 'th1 'th2) (vector 'th1dot 'th2dot))))

... mess, but succeeds.

|#
