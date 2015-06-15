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

;;;;           Variational Mechanics
  
;;; Caution... This file is case sensitive!

(if (and (not
	  (lexical-unreferenceable? user-initial-environment
				    '*parser-canonicalize-symbols?*))
	 (not *parser-canonicalize-symbols?*))
    (local-assignment (the-environment)
		      'Lagrange-interpolation-function
		      lagrange-interpolation-function))


;;; However, there are alternative names for the actual data types.

(define coordinate-tuple column)
(define velocity-tuple column)
(define acceleration-tuple column)
(define momentum-tuple row)

;;; Lagrangian mechanics requires a configuration
;;;   space Q, and a function L:RxQxQ' --> R

;;; Mechanical systems have state at each instant.  The state is the
;;;  information required, along with the equations of motion, to
;;;  determine the future of the system.

;;; At every instant a system has a kinematic state, which has the
;;;  time, the configuration, and the rate of change of the
;;;  configuration.  Lagrangian mechanics is formulated in terms of
;;;  the kinematic state.

;;; Kinematic states and their derivatives are represented as Scheme
;;; vectors, with components time, configuration, and derivatives.

(define (->local t q qdot . derivs)
  (apply vector t q qdot derivs))

(define ->state ->local)

(define ->L-state ->local)

(define (state->n-dof state)
  (let ((q (vector-ref state 1)))
    (if (column? q)
	(s:length q)
	1)))


;;; Selectors are provided for the components of a state.

(define (state->t state)
  (if (not (and (vector? state) (fix:> (vector-length state) 0)))
      (error "Cannot extract time from" state))
 (ref state 0))

(define (state->q state)
  (if (not (and (vector? state) (fix:> (vector-length state) 1)))
      (error "Cannot extract coordinate from" state))
  (ref state 1))

(define (state->qdot state)
  (if (not (and (vector? state) (fix:> (vector-length state) 2)))
      (error "Cannot extract velocity from" state))
  (ref state 2))

(define (state->qddot state)
  (if (not (and (vector? state) (fix:> (vector-length state) 3)))
      (error "Cannot extract acceleration from" state))
  (ref state 3))
    
(define time state->t)
(define coordinate state->q)
(define velocity state->qdot)
(define acceleration state->qddot)

;;; Hamiltonian mechanics requires a phase
;;;   space QxP, and a function H:RxQxP --> R

;;; A system has a dynamic state, which has the time, the
;;; configuration, and the momenta.  Hamiltonian mechanics is
;;; formulated in terms of the dynamic state.

(define (->H-state t q p)
  (vector t q p))

(define (state->p state)
  (if (not (and (vector? state) (fix:> (vector-length state) 2)))
      (error "Cannot extract momentum from" state))
  (ref state 2))

(define momentum state->p)


;;; It is sometimes convenient to be able to split a state.

(define (->parts state cont)
  ;; cont = (lambda (t q p-or-qdot) ...)
  (apply cont (vector->list state)))


(define (with-state cont)
  (lambda (state)
    (->parts state cont)))

(define (with-dynamic-state cont)
  (lambda (state)
    (->parts state
	     (lambda (t q pv)
	       (cont t
		     q
		     (if (vector? pv)
			 (vector->row pv)
			 pv))))))

(define (flatten-state state)
  (list->vector
   (apply append
	  (map (lambda (x)
		 (cond ((vector? x)
			(vector->list x))
		       ((column? x)
			(vector->list (column->vector x)))
		       ((row? x)
			(vector->list (row->vector x)))
		       (else (list x))))
	       (vector->list state)))))

(define local->istate flatten-state)
(define H-state->istate flatten-state)

(define (->istate . args) (flatten-state (apply ->state args)))

(define (flat-state->t fstate) (vector-ref fstate 0))

;;; Assumes k state vector components: [t q q' ... q^(k)]
;;;  Needs to know n (degrees of freedom).

(define (unflatten-L-state flat-state #!optional n)
  (if (default-object? n)
      (set! n (quotient (fix:- (vector-length flat-state) 1) 2)))
  (let* ((kn+1 (vector-length flat-state))
	 (kn (fix:- kn+1 1))
	 (k (quotient kn n)))
    (assert (fix:= 0 (remainder kn n)))
    (if (fix:= n 1)
	flat-state
	(v:generate (fix:+ k 1)		; ->L-state
		    (lambda (i)
		      (if (fix:= i 0)
			  (vector-ref flat-state 0)
			  (vector->column
			   (subvector flat-state
				      (fix:+ (fix:* (fix:- i 1) n) 1)
				      (fix:+ (fix:* i n) 1)))))))))

(define istate->local unflatten-L-state)
(define istate->t time)


;;; Assumes that only t,q,p are present.

(define (unflatten-H-state flat-state)
  (let* ((2n+1 (vector-length flat-state))
	 (2n (fix:- 2n+1 1))
	 (n (quotient 2n 2)))
    (assert (odd? 2n+1))
    (if (fix:= n 1)
	flat-state
	(->H-state
	 (vector-ref flat-state 0)
	 (vector->column (subvector flat-state 1 (fix:+ n 1)))
	 (vector->row (subvector flat-state (fix:+ n 1) 2n+1))))))

(define istate->H-state unflatten-H-state)

;;; Partial derivatives with respect to parts
;;;  of the state are defined for L functions
;;;  (on kinematic states), for H functions (on 
;;;  dynamic states), and for LH functions (the 
;;;  union of these).


;;; Aliases

(define partial_q (partial 1))

(define partial_qdot (partial 2))

(define partial_t (partial 0))

(define partial_p (partial 2))

(define D     derivative)
(define pDq    partial_q)
(define pDqdot partial_qdot)
(define pDp    partial_p)
(define pDt    partial_t)

(define Q     state->q)
(define Qdot  state->qdot)
(define P     state->p)

;;;; Chapter 1

;;; Paths in the configuration manifold are functions that give a
;;; configuration for each time.  From such a path we can construct a
;;; path in the kinematic state space.  If such a path is described 
;;; in terms of generalized coordinates, we have

#|
(define (path->state-path q)
  (lambda (t)
    (->local t
	     (q t)
	     ((D q) t))))
|#

(define (path->state-path q #!optional n)
  (if (default-object? n)
      (set! n 3)
      (assert (fix:> n 1)))
  (lambda (t)
    (list->vector
     (cons t
	   (cons (q t)
		 (let lp ((i (fix:- n 2)) (fi (D q)))
		   (if (fix:= i 0)
		       '()
		       (cons (fi t)
			     (lp (- i 1)
				 (D fi))))))))))

(define Gamma path->state-path)


#|
;;; Another way to make Gamma

(define ((path->state-path q #!optional n) t)
  (if (default-object? n) (set! n 3))
  (list->vector
   (stream-head
    (cons-stream t
		 (cons-stream (q t)
			      (map-stream (lambda (e) (((expt D e) q) t)) 
					  natural-number-stream)))
    n)))
|#

#|
;;; Can we do it this way?  No...  
;;;  We don't know number of degrees of freedom when we build state vector.

(define (path->state q)
  (->local identity
	   q
	   (D q)))
|#

;;; A Lagrangian is an example of an L-function.
;;; An L-function takes  a scalar argument and 2 vector arguments
;;; (t, q, q-dot).  An L-function produces a scalar result.

(define (make-Lagrangian kinetic-energy potential-energy)
  (- kinetic-energy potential-energy))

#|
(define ((L-free-particle mass) local)
  (let ((v (velocity local)))
    (* 1/2 mass (square v))))

(show-expression
 ((L-free-particle 'm)
  (->local 't
	   (coordinate-tuple 'x 'y 'z)
	   (velocity-tuple 'xdot 'ydot 'zdot))))
(+ (* 1/2 m (expt xdot 2))
   (* 1/2 m (expt ydot 2))
   (* 1/2 m (expt zdot 2)))

(show-expression
 ((compose
   (L-free-particle 'm)
   (Gamma (coordinate-tuple (literal-function 'x)
			    (literal-function 'y)
			    (literal-function 'z))))
  't))
(+ (* 1/2 (expt ((D x) t) 2) m)
   (* 1/2 (expt ((D y) t) 2) m)
   (* 1/2 (expt ((D z) t) 2) m))
|#

(define (Lagrangian-action L q t1 t2)
  (definite-integral (compose L (Gamma q))
                     t1 t2))

#|
(define (test-path t)
  (coordinate-tuple (+ (* 4 t) 7)
		    (+ (* 3 t) 5)
		    (+ (* 2 t) 1)))

(Lagrangian-action (L-free-particle 3) test-path 0 10)
;Value: 435.

(define ((variation nu t1 t2 h) t)
  (* h (- t t1) (- t t2) (nu t)))

(define ((varied-free-particle-action mass path nu t1 t2) h)
  (let ((dpath (variation nu t1 t2 h)))
    (Lagrangian-action (L-free-particle mass)
                       (+ path dpath)
                       t1
                       t2)))

((varied-free-particle-action 3.0 test-path 
                              (coordinate-tuple sin cos square)
                              0.0 10.0)
 0.001)
;Value: 436.29121428571153

(minimize
 (varied-free-particle-action 3.0 test-path 
                              (coordinate-tuple sin cos square) 
                              0.0 10.0)
 -2.0 1.0)
;Value: (-1.5987211554602254e-14 435.0000000000237 5)
|#

;;; Equal-spaced times.

(define (make-path t0 q0 t1 q1 qs)
  (let ((n (length qs)))
    (let ((ts (linear-interpolants t0 t1 n)))
      (Lagrange-interpolation-function
       (append (list q0) qs (list q1))
       (append (list t0) ts (list t1))))))

(define ((parametric-path-action Lagrangian t0 q0 t1 q1) qs)
  (let ((path (make-path t0 q0 t1 q1 qs)))
    (Lagrangian-action Lagrangian path t0 t1)))

(define (find-path Lagrangian t0 q0 t1 q1 n)
  (let ((initial-qs (linear-interpolants q0 q1 n)))
    (let ((minimizing-qs
	   (multidimensional-minimize
	    (parametric-path-action Lagrangian t0 q0 t1 q1)
	    initial-qs)))
      (make-path t0 q0 t1 q1 minimizing-qs))))

(define (linear-interpolants x0 x1 n)
  (let ((dx (- x1 x0)) (n+1 (fix:+ n 1)))
    (let lp ((i 1) (xs '()))
      (if (fix:> i n)
	  (reverse xs)
	  (lp (fix:+ i 1)
	      (cons (+ x0 (/ (* i dx) n+1)) xs))))))

#|
;;; For example, consider the Harmonic oscillator with
;;;  spring constant, k, and mass, m.


(define ((L-harmonic m k) local)
  (let ((q (coordinate local)) 
        (v (velocity local)))
    (- (* 1/2 m (square v))
       (* 1/2 k (square q)))))

(define q
  (find-path (L-harmonic 1.0 1.0) 0. 1. pi/2 0. 3))


(define p (frame 0.0 pi/2 -3e-4 3e-4))

(plot-function p 
 (lambda (t)
   (- (q t) (cos t)))
 0.0
 pi/2
 .01)

(graphics-clear p)
(graphics-close p)

;;; Error curve is in actionerror3.xwd

;;; ********************** too slow...

(define q1
  (show-time
   (lambda ()
     (find-path (L-harmonic 1.0 1.0) 0. 1. pi/2 0. 4))))
;;;maharal process time: 22310 (20740 RUN + 1570 GC); real time: 22316
;;;PPA process time: 160700 (153800 RUN + 6900 GC); real time: 161913


(define p (frame 0.0 pi/2 -1e-4 1e-4))

(plot-function p 
 (lambda (t)
   (- (q1 t) (cos t)))
 0.0
 pi/2
 .01)

;;; Error curve is in actionerror4.xwd
|#

#|
(graphics-close win2)
(define win2 (frame 0. pi/2 0. 1.2))

(define ((parametric-path-action Lagrangian t0 q0 t1 q1) 
         intermediate-qs)
    (let ((path (make-path t0 q0 t1 q1 intermediate-qs)))
      ;; display path
      (graphics-clear win2)
      (plot-function win2 path t0 t1 (/ (- t1 t0) 100))
      ;; compute action
      (Lagrangian-action Lagrangian path t0 t1)))

(define ((parametric-path-action Lagrangian t0 q0 t1 q1) 
         intermediate-qs)
    (let ((path (make-path t0 q0 t1 q1 intermediate-qs)))
      ;; display path
      ;(graphics-clear win2)
      (plot-function win2 path t0 t1 (/ (- t1 t0) 100))
      ;; compute action
      (Lagrangian-action Lagrangian path t0 t1)))

(find-path (L-harmonic 1. 1.) 0. 1. pi/2 0. 3)

(graphics-close win2)
|#

;;; Given a Lagrangian, we can obtain Lagrange's equations of motion.

(define ((Lagrange-equations Lagrangian #!optional dissipation-function) q)
  (let ((state-path (Gamma q)))
    (if (default-object? dissipation-function)
	(- (D (compose ((partial 2) Lagrangian) state-path))
	   (compose ((partial 1) Lagrangian) state-path))
	(- (D (compose ((partial 2) Lagrangian) state-path))
	   (compose ((partial 1) Lagrangian) state-path)
	   (- (compose ((partial 2) dissipation-function) state-path))))))


#|
(define ((Lagrange-equations Lagrangian) q)
  (let ((local-path (Gamma q)))
    (- (D (compose ((partial 2) Lagrangian) local-path))
       (compose ((partial 1) Lagrangian) local-path))))
|#


#|
(define (test-path t)
  (coordinate-tuple (+ (* 'a t) 'a0)
		    (+ (* 'b t) 'b0)
		    (+ (* 'c t) 'c0)))

(print-expression
 (((Lagrange-equations (L-free-particle 'm))
   test-path)
  't))
(down 0 0 0)

(show-expression
 (((Lagrange-equations (L-free-particle 'm))
   (literal-function 'x))
  't))
(* m (((expt D 2) x) t))

(show-expression
 (((Lagrange-equations (L-harmonic 'm 'k))
   (literal-function 'x))
  't))
(+ (* k (x t)) (* m (((expt D 2) x) t)))

(show-expression
 (((Lagrange-equations (L-harmonic 'm 'k))
   (lambda (t) (* 'a (cos (+ (* 'omega t) 'phi)))))
  't))
(+ (* a k (cos (+ (* omega t) phi)))
   (* -1 a m (expt omega 2) (cos (+ (* omega t) phi))))
|#

#|
(define ((L-uniform-acceleration m g) local)
  (let ((q (coordinate local))
        (v (velocity local)))
    (let ((y (ref q 1)))
      (- (* 1/2 m (square v)) (* m g y)))))

(show-expression
 (((Lagrange-equations
    (L-uniform-acceleration 'm 'g))
   (coordinate-tuple (literal-function 'x)
		     (literal-function 'y)))
  't))
(down (* m (((expt D 2) x) t))
      (+ (* g m) (* m (((expt D 2) y) t))))


(define ((L-central-rectangular m V) local)
  (let ((q (coordinate local))
        (v (velocity local)))
    (- (* 1/2 m (square v))
       (V (sqrt (square q))))))

(show-expression
 (((Lagrange-equations 
    (L-central-rectangular 'm (literal-function 'V)))
   (coordinate-tuple (literal-function 'x) (literal-function 'y)))
  't))
(down
 (+ (* m (((expt D 2) x) t))
    (/ (* ((D V) (sqrt (+ (expt (x t) 2) (expt (y t) 2)))) (x t))
       (sqrt (+ (expt (x t) 2) (expt (y t) 2)))))
 (+ (* m (((expt D 2) y) t))
    (/ (* ((D V) (sqrt (+ (expt (x t) 2) (expt (y t) 2)))) (y t))
       (sqrt (+ (expt (x t) 2) (expt (y t) 2))))))
|#

#|
;;; Consider planar motion in a central force field, with an arbitrary
;;; potential, U, depending only on the radius.  The generalized
;;; coordinates are polar. 


(define ((L-central-polar m V) local)
  (let ((q (coordinate local))
        (qdot (velocity local)))
    (let ((r (ref q 0))
          (phi (ref q 1))
          (rdot (ref qdot 0))
          (phidot (ref qdot 1)))
      (- (* 1/2 m
           (+ (square rdot)
              (square (* r phidot))) )
         (V r)))))

(show-expression
 (((Lagrange-equations
    (L-central-polar 'm (literal-function 'V)))
   (coordinate-tuple (literal-function 'r)
		     (literal-function 'phi)))
  't))
(down
 (+ (* -1 m (r t) (expt ((D phi) t) 2))
    (* m (((expt D 2) r) t))
    ((D V) (r t)))
 (+ (* 2 m ((D r) t) (r t) ((D phi) t))
    (* m (((expt D 2) phi) t) (expt (r t) 2))))
|#

;;; Coordinate Transformation to State Transformation


#|

;;; if defined F as state function F(t, q, v); (partial 2) F = 0

(define (F->C F)
  (define (C state)
    (->local (time state)
	     (F state)
	     (+ (((partial 0) F) state)
		(* (((partial 1) F) state) 
		   (velocity state)))))
  C)

;;; with q = F(t, q')

(define ((F->C F) local)
  (let ((t (time local))
        (q (coordinate local))
        (v (velocity local)))
    (->local t
             (F t q)
             (+ (((partial 0) F) t q)
                (* (((partial 1) F) t q) v)))))


;;; alternate in text

(define (F->C F)
  (Gamma-bar
   (lambda (q)
     (Gamma
      (lambda (t) (F t (q t)))))))

(define ((F->C F) local)
  (let ((n (vector-length local)))
    ((Gamma-bar
      (lambda (q)
	(Gamma
	 (lambda (t) (F t (q t)))
	 n)))
     local)))

|#

;;; current definition of F->C

(define ((F->C F) local)
  (let ((n (vector-length local)))
    ((Gamma-bar
      (lambda (qp)
	(Gamma
          (compose F (Gamma qp))
	 n)))
     local)))

#|

;;; version for display in text

(define (F->C F)
  (define (f-bar q-prime)
    (define q
      (compose F (Gamma q-prime)))
    (Gamma q))
  (Gamma-bar f-bar))

|#


;;; The following transformations are applicable to 
;;;  configuration coordinates. 

#|
(define (r->p t rectangular-tuple)
  (let ((x (ref rectangular-tuple 0))
        (y (ref rectangular-tuple 1)))
    (let ((r (sqrt (+ (square x) (square y))))
          (phi (atan y x)))
      (up r phi))))
|#

(define (r->p tqv)
  (let ((rectangular-tuple (coordinate tqv)))
    (let ((x (ref rectangular-tuple 0))
	  (y (ref rectangular-tuple 1)))
      (let ((r (sqrt (+ (square x) (square y))))
	    (phi (atan y x)))
	(up r phi)))))

#|
(define (p->r t polar-tuple)
  (let ((r (ref polar-tuple 0)) 
        (phi (ref polar-tuple 1)))
    (let ((x (* r (cos phi))) 
          (y (* r (sin phi))))
      (up x y))))
|#

(define (p->r tqv)
  (let ((polar-tuple (coordinate tqv)))
    (let ((r (ref polar-tuple 0)) 
          (phi (ref polar-tuple 1)))
      (let ((x (* r (cos phi))) 
            (y (* r (sin phi))))
        (up x y)))))

#|

(show-expression 
 (velocity
  ((F->C p->r)
   (->local 't 
	    (coordinate-tuple 'r 'phi) 
	    (velocity-tuple 'rdot 'phidot)))))
(up (+ (* -1 r phidot (sin phi)) (* rdot (cos phi)))
    (+ (* r phidot (cos phi)) (* rdot (sin phi))))


(define (L-central-polar m V)
  (compose (L-central-rectangular m V)
	   (F->C p->r)))

(show-expression
  ((L-central-polar 'm (literal-function 'V))
   (->local 't (coordinate-tuple 'r 'phi) 
               (velocity-tuple 'rdot 'phidot))))
(+ (* 1/2 m (expt phidot 2) (expt r 2))
   (* 1/2 m (expt rdot 2))
   (* -1 (V r)))
|#

(define ((Gamma-bar f) local)
  ((f (osculating-path local)) (time local)))

;;; An alternative method allows taking derivatives in the
;;; construction of the Lagrangian.

(define ((osculating-path state0) t)
  (let ((t0 (time state0))
	(q0 (coordinate state0))
	(k (vector-length state0)))
    (let ((dt (- t t0)))
      (let loop ((n 2) (sum q0) (dt^n/n! dt))
	(if (fix:= n k)
	    sum
	    (loop (+ n 1)
		  (+ sum (* (vector-ref state0 n) dt^n/n!))
		  (/ (* dt^n/n! dt) n)))))))

#|
(define ((L-pend-1 m l g ys) state)
  (let ((t0 (time state)))
    (let ((osc-theta (osculating-path state)))
      (let ((x (* l (sin osc-theta)))
	    (y (- ys (* l (cos osc-theta)))))
	(let ((xdot0 ((D x) t0))
	      (ydot0 ((D y) t0)))
	  (- (* 1/2 m
		(+ (square xdot0)
		   (square ydot0)))
	     (* m g (y t0))))))))

(show-expression
 (((Lagrange-equations (L-pend-1 'm 'l 'g (literal-function 'y_s)))
   (literal-function 'theta))
  't))
(+ (* g l m (sin (theta t)))
   (* (expt l 2) m (((expt D 2) theta) t))
   (* l m (sin (theta t)) (((expt D 2) y_s) t)))
|#

#|
;;; Driven pendulum example

(define ((T-pend m l g ys) local)
  (let ((t (time local))
        (theta (coordinate local))
        (thetadot (velocity local)))
    (let ((ysdot (D ys)))
      (* 1/2 m
         (+ (square (* l thetadot))
            (square (ysdot t))
            (* 2 (ysdot t) l (sin theta) thetadot))))))

(define ((V-pend m l g ys) local)
  (let ((t (time local))
        (theta (coordinate local)))
    (* m g (- (ys t) (* l (cos theta))))))

(define L-pend (- T-pend V-pend))

(show-expression
 ((L-pend 'm 'l 'g (literal-function 'y_s))
  (->local 't 'theta 'thetadot)))
(+ (* 1/2 (expt l 2) m (expt thetadot 2))
   (* l m thetadot ((D y_s) t) (sin theta))
   (* g l m (cos theta))
   (* -1 g m (y_s t))
   (* 1/2 m (expt ((D y_s) t) 2)))

(show-expression
 (((Lagrange-equations
    (L-pend 'm 'l 'g (literal-function 'y_s)))
   (literal-function 'theta))
  't))
(+ (* g l m (sin (theta t)))
   (* (expt l 2) m (((expt D 2) theta) t))
   (* l m (((expt D 2) y_s) t) (sin (theta t))))
|#

#|
;;; Same driven pendulum by coordinate transformation

(define ((Lf m g) local)
  (let ((q (coordinate local))
        (v (velocity local)))
    (let ((h (ref q 1)))
      (- (* 1/2 m (square v)) (* m g h)))))

(define ((dp-coordinates l y_s) local)
  (let ((t (time local))
	(theta (coordinate local)))
    (let ((x (* l (sin theta)))
	  (y (- (y_s t) (* l (cos theta)))))
      (coordinate-tuple x y))))

(define (L-pend m l g y_s)
  (compose (Lf m g) 
           (F->C (dp-coordinates l y_s))))

(show-expression
 ((L-pend 'm 'l 'g (literal-function 'y_s))
  (->local 't 'theta 'thetadot)))
(+ (* 1/2 (expt l 2) m (expt thetadot 2))
   (* l m thetadot (sin theta) ((D y_s) t))
   (* g l m (cos theta))
   (* -1 g m (y_s t))
   (* 1/2 m (expt ((D y_s) t) 2)))

(show-expression
 (((Lagrange-equations
    (L-pend 'm 'l 'g (literal-function 'y_s)))
   (literal-function 'theta))
  't))
(+ (* g l m (sin (theta t)))
   (* (expt l 2) m (((expt D 2) theta) t))
   (* l m (((expt D 2) y_s) t) (sin (theta t))))
|#

;;; "total time derivative"

#|
(define (Dt-procedure F)
  (define (DF-on-path q)
    (D (compose F (Gamma q))))
  (Gamma-bar DF-on-path))
|#

(define ((Dt-procedure F) state)
  (let ((n (vector-length state)))
    (define (DF-on-path q)
      (D (compose F (Gamma q n))))
    ((Gamma-bar DF-on-path) state)))

(define Dt
  (make-operator Dt-procedure 'Dt))

#|
(print-expression
 ((Dt
   (lambda (state)
     (let ((t (time state))
	   (q (coordinate state)))
       (square q))))
  (->local
   't
   (coordinate-tuple 'x 'y) 
   (velocity-tuple 'vx 'vy))))
(+ (* 2 vx x) (* 2 vy y))


(print-expression
 ((Dt (Dt (lambda (state) (coordinate state))))
  (->local 't 'x 'v 'a 'j)))
a

(print-expression
 ((Dt (Dt (lambda (state)
	    (square (coordinate state)))))
  (->local 't 'x 'v 'a 'j)))
(+ (* 2 a x) (* 2 (expt v 2)))
|#

(define (Euler-Lagrange-operator Lagrangian)
  (- (Dt ((partial 2) Lagrangian))
     (compose ((partial 1) Lagrangian) trim-last-argument)))

(define (trim-last-argument local)
  (apply up (except-last-pair (vector->list (up->vector local)))))

(define LE Euler-Lagrange-operator)
(define Lagrange-equations-operator LE)

#|
(print-expression
 ((LE (L-harmonic 'm 'k))
  (->local 't 'x 'v 'a)))
(+ (* a m) (* k x))

(print-expression
 ((LE (L-harmonic 'm 'k))
  (->local 't
	   #(x y)
	   #(vx vy)
	   #(ax ay))))
(down (+ (* ax m) (* k x))
      (+ (* ay m) (* k y)))


;;; Adding extra state components is harmless.

(print-expression
 ((LE (L-harmonic 'm 'k))
  (->local 't 'x 'v 'a 'j)))
(+ (* a m) (* k x))

;;; But watch out.  If not enuf local componenents
;;;  are specified we lose.

(print-expression
 ((LE (L-harmonic 'm 'k))
  (->local 't 'x 'v)))
;;; error

(print-expression
 ((compose (LE (L-central-polar 'm (literal-function 'V)))
	   (Gamma
	    (coordinate-tuple (literal-function 'r)
			      (literal-function 'phi))
	    4))
  't))
(down
 (+ (* -1 m (r t) (expt ((D phi) t) 2)) (* m (((expt D 2) r) t)) ((D V) (r t)))
 (+ (* m (((expt D 2) phi) t) (expt (r t) 2))
    (* 2 m ((D r) t) (r t) ((D phi) t))))

#|
(define ((generalized-LE Lagrangian) state)
  (let lp ((i 1))
    (if (= i (s:length state))
	Lagrangian
	(- (Dt (lp (+ i 1)))
	   (compose ((partial i) Lagrangian)
		    trim-last-argument)))))
|#

(define ((generalized-LE Lagrangian) state)
  (let ((m (s:length state)))
    (assert (and (fix:> m 3) (even? m))
	    "Incorrect state size for Lagrange Equations")
    (let lp ((i (quotient m 2)))
      (if (fix:= i 0)
	  0
	  (- ((compose ((expt Dt (fix:- i 1))
			((partial i) Lagrangian))
		       (iterated trim-last-argument
				 (fix:- (quotient m 2) i)))
	      state)
	     (lp (fix:- i 1)))))))

  
(define ((L2harmonic m k) state)
  (let ((x (coordinate state))
	(a (acceleration state)))
    (+ (* 1/2 m x a) (* 1/2 k (square x)))))

(print-expression
 ((generalized-LE (L2harmonic 'm 'k))
  (->local 't 'x 'v 'a 'j 'p)))
(+ (* -1 a m) (* -1 k x))

(define L (literal-function 'L (-> (UP* Real) Real)))

(pe ((generalized-LE L) (->local 't 'x 'v 'a)))
(+ (* a (((partial 2) ((partial 2) L)) (up t x v)))
   (* v (((partial 1) ((partial 2) L)) (up t x v)))
   (((partial 0) ((partial 2) L)) (up t x v))
   (* -1 (((partial 1) L) (up t x v))))

(pe ((generalized-LE L) (->local 't 'x 'v 'a 'j 'p)))
(+ (* (expt a 2) (((partial 2) ((partial 2) ((partial 3) L))) (up t x v a)))
   (* 2 a j (((partial 2) ((partial 3) ((partial 3) L))) (up t x v a)))
   (* 2 a v (((partial 1) ((partial 2) ((partial 3) L))) (up t x v a)))
   (* (expt j 2) (((partial 3) ((partial 3) ((partial 3) L))) (up t x v a)))
   (* 2 j v (((partial 1) ((partial 3) ((partial 3) L))) (up t x v a)))
   (* (expt v 2) (((partial 1) ((partial 1) ((partial 3) L))) (up t x v a)))
   (* 2 a (((partial 0) ((partial 2) ((partial 3) L))) (up t x v a)))
   (* a (((partial 1) ((partial 3) L)) (up t x v a)))
   (* -1 a (((partial 2) ((partial 2) L)) (up t x v a)))
   (* 2 j (((partial 0) ((partial 3) ((partial 3) L))) (up t x v a)))
   (* p (((partial 3) ((partial 3) L)) (up t x v a)))
   (* 2 v (((partial 0) ((partial 1) ((partial 3) L))) (up t x v a)))
   (* -1 v (((partial 1) ((partial 2) L)) (up t x v a)))
   (((partial 0) ((partial 0) ((partial 3) L))) (up t x v a))
   (* -1 (((partial 0) ((partial 2) L)) (up t x v a)))
   (((partial 1) L) (up t x v a)))
|#

#|
;;; Coupled harmonic oscillators.

(define ((L-coupled-harmonic m k) state)
  (let ((q (coordinate state))
	(qdot (velocity state)))
    (- (* 1/2 qdot m qdot)
       (* 1/2 q k q))))

(show-expression
 (((Lagrange-equations
    (L-coupled-harmonic (down (down 'm_1 0) (down 0 'm_2))
			(down (down 'k_1 'c) (down 'c 'k_2))))
   (coordinate-tuple (literal-function 'x)
		     (literal-function 'y)))
  't))
(down (+ (* c (y t)) (* k_1 (x t)) (* m_1 (((expt D 2) x) t)))
      (+ (* c (x t)) (* k_2 (y t)) (* m_2 (((expt D 2) y) t))))
|#

#|
;;; Pendulum of mass m2 and length b, hanging from a support of mass
;;; m1 that is free to move horizontally (from Groesberg, Advanced
;;; Mechanics, p. 72) 

(define ((L-sliding-pend m1 m2 b g) state)
  (let ((q (coordinate state))
	(qdot (velocity state)))
    (let* ((x (ref q 0))
	   (xdot (ref qdot 0))
	   (theta (ref q 1))
	   (thetadot (ref qdot 1))
	   (rel-pend-vel
	    (* b thetadot (velocity-tuple (cos theta) (sin theta))))
	   (pend-vel (+ rel-pend-vel (velocity-tuple xdot 0)))
	   (Tpend (* 1/2 m2 (square pend-vel)))
	   (Tsupport (* 1/2 m1 (square xdot)))
	   (V (- (* m2 g b (cos theta)))))
      (+ Tpend Tsupport (- V)))))

(show-expression
 (((Lagrange-equations (L-sliding-pend 'm_1 'm_2 'b 'g))
   (coordinate-tuple (literal-function 'x)
		     (literal-function 'theta)))
  't))
(down
 (+ (* -1 b m_2 (sin (theta t)) (expt ((D theta) t) 2))
    (* b m_2 (((expt D 2) theta) t) (cos (theta t)))
    (* m_1 (((expt D 2) x) t))
    (* m_2 (((expt D 2) x) t)))
 (+ (* (expt b 2) m_2 (((expt D 2) theta) t))
    (* b g m_2 (sin (theta t)))
    (* b m_2 (((expt D 2) x) t) (cos (theta t)))))
|# 

#|
;;; Consider a simple pendulum with Rayleigh dissipation:

(define ((L-pendulum g m l) state)
  (let ((theta (coordinate state))
	(thetadot (velocity state)))
    (+ (* 1/2 m (square (* l thetadot)))
       (* g m l (cos theta)))))

(define ((Rayleigh-dissipation k) state)
  (let ((qdot (velocity state)))
    (* qdot k qdot)))

(show-expression
 (((Lagrange-equations (L-pendulum 'g 'm 'l)
		       (Rayleigh-dissipation 'k))
   (literal-function 'theta))
  't))
(+ (* 2 k ((D theta) t))
   (* g l m (sin (theta t)))
   (* (expt l 2) m (((expt D 2) theta) t)))
|#

#|
;;; Can group coordinates.  Procedures don't care.

(define ((L-two-particle m1 m2) local)
  (let ((x (coordinate local))
	(v (velocity local))
	(V (literal-function 'V (-> (X (^ Real 2) (^ Real 2)) Real))))
    (let ((x1 (ref x 0)) (x2 (ref x 1))
          (v1 (ref v 0)) (v2 (ref v 1)))
      (- (+ (* 1/2 m1 (square v1))
	    (* 1/2 m2 (square v2)))
	 (V x1 x2)))))

(se
 (((Lagrange-equations (L-two-particle 'm_1 'm_2))
   (coordinate-tuple
    (coordinate-tuple (literal-function 'x_1) (literal-function 'y_1))
    (coordinate-tuple (literal-function 'x_2) (literal-function 'y_2))))
  't))
(down
 (down
  (+ (* m_1 (((expt D 2) x_1) t))
     (((partial 0 0) V) (up (x_1 t) (y_1 t)) (up (x_2 t) (y_2 t))))
  (+ (* m_1 (((expt D 2) y_1) t))
     (((partial 0 1) V) (up (x_1 t) (y_1 t)) (up (x_2 t) (y_2 t)))))
 (down
  (+ (* m_2 (((expt D 2) x_2) t))
     (((partial 1 0) V) (up (x_1 t) (y_1 t)) (up (x_2 t) (y_2 t))))
  (+ (* m_2 (((expt D 2) y_2) t))
     (((partial 1 1) V) (up (x_1 t) (y_1 t)) (up (x_2 t) (y_2 t))))))
|#

;;; An alternative way to obtain Lagrange's equations arises from 
;;;  expanding the derivative of the momentum by the chain rule to
;;;  get the Lagrange operator.  Lagrange's equations are then
;;;  obtained by calling the Lagrange operator on the objects
;;;  q, qdot, qddot, all functions of time.

;;; ******* This stuff only works for L(t,q,qd)

(define (on-jet q qdot)
  (lambda (lfun)
    (compose lfun
	     (lambda (t)
	       (->local t
			(q t)
			(qdot t))))))

(define (Lagrange-operator Lagrangian)
  (let ((P ((partial 2) Lagrangian))
	(F ((partial 1) Lagrangian)))
    (let ((Pq ((partial 1) P))
	  (Pqdot ((partial 2) P))
	  (Pt ((partial 0) P)))
      (lambda (q qdot qddot)
	(let ((lift (on-jet q qdot)))
	  (+ (* (lift Pqdot) qddot)
	     (* (lift Pq) qdot)
	     (lift Pt)
	     (- (lift F))))))))

(define (Lagrange-equations-from-operator Lagrangian)
  (let ((lop (Lagrange-operator Lagrangian)))
    (lambda (q)
      (let ((qdot (D q)))
	(let ((qddot (D qdot)))
	  (lop q qdot qddot))))))

#|
(show-expression
 (((Lagrange-equations-from-operator
    (L-sliding-pend 'm_1 'm_2 'b 'g))
   (coordinate-tuple (literal-function 'x)
		     (literal-function 'theta)))
  't))
(down
 (+ (* -1 b m_2 (sin (theta t)) (expt ((D theta) t) 2))
    (* b m_2 (((expt D 2) theta) t) (cos (theta t)))
    (* m_1 (((expt D 2) x) t))
    (* m_2 (((expt D 2) x) t)))
 (+ (* (expt b 2) m_2 (((expt D 2) theta) t))
    (* b g m_2 (sin (theta t)))
    (* b m_2 (((expt D 2) x) t) (cos (theta t)))))
|#

;;; For integrating Lagrange's equations we need them in a form which 
;;; has the highest derivative isolated.
;;; The following is an explicit solution for the second-derivative 
;;; from Lagrange's equations, based on the operator form above:

#|
(define (Lagrange-explicit Lagrangian)
  (let ((P ((partial 2) Lagrangian))
        (F ((partial 1) Lagrangian)))
    (let ((dP/dq ((partial 1) P))
          (dP/dqdot ((partial 2) P))
          (dP/dt ((partial 0) P))
	  (qdot state->qdot))
      (/ (- F (+ (* dP/dq qdot) dP/dt))
	 dP/dqdot))))

(define (Lagrangian->acceleration L)
  (let ((P ((partial 2) L))
        (F ((partial 1) L)))
    (/ (- F
          (+ ((partial 0) P) 
             (* ((partial 1) P) velocity)))
       ((partial 2) P))))
|#

(define ((Lagrangian->acceleration L) state)
  (let ((P ((partial 2) L))
        (F ((partial 1) L)))
    (* (s:inverse (velocity state)
		  (((partial 2) P) state)
		  (velocity state))
       ((- F
	   (+ ((partial 0) P) 
	      (* ((partial 1) P) velocity)))
	state))))

(define Lagrange-explicit Lagrangian->acceleration)

#|
;;; Thus, for example, we can obtain the general form of the vector
;;; of accelerations as a function of the positions, and velocities:

(show-expression
 ((Lagrange-explicit (L-sliding-pend 'm_1 'm_2 'b 'g))
  (->local 't
	   (coordinate-tuple 'x 'theta)
	   (velocity-tuple 'xdot 'thetadot))))
(up
 (+
  (/ (* b m_2 (expt thetadot 2) (sin theta))
     (+ (* m_2 (expt (sin theta) 2)) m_1))
  (/ (* g m_2 (sin theta) (cos theta))
     (+ (* m_2 (expt (sin theta) 2)) m_1)))
 (+
  (/ (* -1 m_2 (expt thetadot 2) (sin theta) (cos theta))
     (+ (* m_2 (expt (sin theta) 2)) m_1))
  (/ (* -1 g m_1 (sin theta))
     (+ (* b m_2 (expt (sin theta) 2)) (* b m_1)))
  (/ (* -1 g m_2 (sin theta))
     (+ (* b m_2 (expt (sin theta) 2)) (* b m_1)))))
|#

(define ((Lagrange-equations-1 L) q v)
  (let ((local-path (qv->local-path q v)))
    (- (D local-path)
       (compose (local-state-derivative L)
		local-path))))

(define ((local-state-derivative L) local)
  (->local 1
	   (velocity local)
	   ((Lagrangian->acceleration L) local)))

(define Lagrange-equations-first-order
        Lagrange-equations-1)
#|
(define ((Lagrange-equations-1 L) q v)
  (let ((local-path (qv->local-path q v)))
    (- (D local-path)
       (->local 1
                (compose velocity local-path)
                (compose (Lagrangian->acceleration L)
                         local-path)))))
|#

(define ((qv->local-path q v) t)
  (->local t (q t) (v t)))

#|
(show-expression
 (((Lagrange-equations-1 (L-harmonic 'm 'k))
   (coordinate-tuple (literal-function 'x)
                     (literal-function 'y))
   (velocity-tuple (literal-function 'v_x)
                   (literal-function 'v_y)))
  't))
(up 0
    (up (+ ((D x) t) (* -1 (v_x t))) (+ ((D y) t) (* -1 (v_y t))))
    (up (+ (/ (* k (x t)) m) ((D v_x) t)) (+ (/ (* k (y t)) m) ((D v_y) t))))
|#

;;; Structured states

(define (Lagrangian->state-derivative L)
  (let ((acceleration (Lagrange-explicit L)))
    (lambda (state)
      (up
       1
       (velocity state)
       (acceleration state)))))



#|
;;; Going back to the driven pendulum...

(define ((periodic-drive amplitude frequency phase) t)
  (* amplitude (cos (+ (* frequency t) phase))))

(define (L-periodically-driven-pendulum m l g a omega)
  (let ((ys (periodic-drive a omega 0)))
    (L-pend m l g ys)))

(show-expression
 (((Lagrange-equations
    (L-periodically-driven-pendulum 'm 'l 'g 'a 'omega))       
   (literal-function 'theta))
  't))
(+ (* -1 a l m (expt omega 2) (sin (theta t)) (cos (* omega t)))
   (* g l m (sin (theta t)))
   (* (expt l 2) m (((expt D 2) theta) t)))

(show-expression
 ((Lagrange-explicit (L-periodically-driven-pendulum 'm 'l 'g 'a 'omega))
  (->local 't 'theta 'thetadot)))
(+ (/ (* a (expt omega 2) (sin theta) (cos (* omega t))) l)
   (/ (* -1 g (sin theta)) l))

(show-expression
 ((Lagrangian->state-derivative
   (L-periodically-driven-pendulum 'm 'l 'g 'a 'omega))
  (up 't 'theta 'thetadot)))
(up
 1
 thetadot
 (+ (/ (* a (expt omega 2) (cos (* omega t)) (sin theta)) l)
    (/ (* -1 g (sin theta)) l)))

(define (pend-sysder m l g a omega)
  (Lagrangian->state-derivative
    (L-periodically-driven-pendulum m l g a omega)))
|#

#|
;;; Using these we can do graphing of trajectories:

(define plot-win (frame 0. 100. -pi pi))

(define ((monitor-theta win) state)
  (let ((theta ((principal-value pi) (coordinate state))))
    (plot-point win (time state) theta)))
  
((evolve pend-sysder 
         1.0                    ;m=1kg
         1.0                    ;l=1m
         9.8                    ;g=9.8m/s$^2$
         0.1                    ;a=1/10 m
         (* 2.0 (sqrt 9.8)) )
 (up 0.0             ;t$_0$=0
     1.22            ;theta$_0$=1 radian
     1e-10)          ;thetadot$_0$=0 radians/s
 (monitor-theta plot-win)
 0.01                     ;step between plotted points
 100.0                    ;final time
 1.0e-13)

;;;(set! *ode-integration-method* 'bulirsch-stoer)
;;;(set! *ode-integration-method* 'qcrk4)

((evolve pend-sysder 
         1.0                    ;m=1kg
         1.0                    ;l=1m
         9.8                    ;g=9.8m/s$^2$
         0.1                    ;a=1/10 m
         (* 2.0 (sqrt 9.8)) )
 (up 0.0             ;t$_0$=0
     1.22            ;theta$_0$=1 radian
     0.0)          ;thetadot$_0$=0 radians/s
 (monitor-theta plot-win)
 0.01                     ;step between plotted points
 100.0                    ;final time
 1.0e-13)

(show-time 
 (lambda ()
   ((evolve pend-sysder 
         1.0                    ;m=1kg
         1.0                    ;l=1m
         9.8                    ;g=9.8m/s$^2$
         0.1                    ;a=1/10 m
         (* 2.0 (sqrt 9.8)) )
 (up 0.0             ;t$_0$=0
     1.22            ;theta$_0$=1 radian
     0.0)          ;thetadot$_0$=0 radians/s
 (monitor-theta plot-win)
 0.01                     ;step between plotted points
 100.0                    ;final time
 1.0e-13)))

;;; qcrk4 -- process time: 47220 (41550 RUN + 5670 GC); real time: 47208
;;; bs    -- process time: 6440 (5680 RUN + 760 GC); real time: 6435

(graphics-clear plot-win)
(graphics-close plot-win)

;;; Picture stored in mechanics/driven-pend-theta.xwd
|#

;;; Given a Lagrangian, we can make an energy function on (t, Q, Qdot). 

(define (Lagrangian->energy L)
  (let ((P ((partial 2) L)))
    (- (* P velocity) L)))


;;; On a trajectory there may be power lost (if dissipation)
;;;  The following produces the power lost.  

(define ((Lagrangian->power-loss L) q)
  (D (compose (Lagrangian->energy L)
	      (Gamma q))))

#|
;;; Alternatively

(define ((Lagrangian->power-loss L) q)
  (- (* ((Lagrange-equations L) q) (D q))
     (compose ((partial 0) L)
	      (Gamma q))))
|#

#|
;;; For example, on a specified trajectory, we can compute the energy,
;;; which turns out to be T+V.

(show-expression
 ((compose
   (Lagrangian->energy (L-central-polar 'm (literal-function 'U)))
   (Gamma 
    (coordinate-tuple (literal-function 'r) (literal-function 'phi))))
  't))
(+ (* 1/2 m (expt (r t) 2) (expt ((D phi) t) 2))
   (* 1/2 m (expt ((D r) t) 2))
   (U (r t)))


;;; In fact, we can see how the energy is conserved:

(show-expression
 (((Lagrangian->power-loss (L-central-polar 'm (literal-function 'U)))
   (coordinate-tuple (literal-function 'r) (literal-function 'phi)))
  't))
(+ (* m (((expt D 2) phi) t) ((D phi) t) (expt (r t) 2))
   (* m (expt ((D phi) t) 2) (r t) ((D r) t))
   (* m (((expt D 2) r) t) ((D r) t))
   (* ((D U) (r t)) ((D r) t)))

;;; This last expression is (nontrivially!) zero on any trajectory
;;; which satisfies Lagrange's equations.
|#

#|
;;; Note, this can be implemented in terms of T-CURVILINEAR.

(define ((T3-spherical m) local)
  (let ((t (time local))
        (q (coordinate local))
        (qdot (velocity local)))
    (let ((r (ref q 0))
          (theta (ref q 1))
          (phi (ref q 2))
          (rdot (ref qdot 0))
          (thetadot (ref qdot 1))
          (phidot (ref qdot 2)))
      (* 1/2 m
        (+ (square rdot)
           (square (* r thetadot))
           (square (* r (sin theta) phidot)))))))

(define (L3-central m Vr)
  (define (Vs local)
    (let ((r (ref (coordinate local) 0)))
      (Vr r)))
  (- (T3-spherical m) Vs))


(show-expression
 (((partial 1) (L3-central 'm (literal-function 'V)))
  (->local 't
           (coordinate-tuple 'r 'theta 'phi)
           (velocity-tuple 'rdot 'thetadot 'phidot))))
(down
 (+ (* m r (expt phidot 2) (expt (sin theta) 2))
    (* m r (expt thetadot 2))
    (* -1 ((D V) r)))
 (* m (expt r 2) (expt phidot 2) (cos theta) (sin theta))
 0)

(show-expression
 (((partial 2) (L3-central 'm (literal-function 'V)))
  (->local 't 
           (coordinate-tuple 'r 'theta 'phi)
           (velocity-tuple 'rdot 'thetadot 'phidot))))
(down (* m rdot)
      (* m (expt r 2) thetadot)
      (* m (expt r 2) phidot (expt (sin theta) 2)))
|#

;;; Another coordinate transformation 

(define (s->r local)
  (let ((q (coordinate local)))
    (let ((r (ref q 0))
	  (theta (ref q 1))
	  (phi (ref q 2)))
      (let ((x (* r (sin theta) (cos phi)))
	    (y (* r (sin theta) (sin phi)))
	    (z (* r (cos theta))))
	(coordinate-tuple x y z)))))

#|
(define ((ang-mom-z m) local)
  (let ((q (coordinate local))
        (v (velocity local)))
     (ref (cross-product q (* m v)) 2)))

(show-expression
  ((compose (ang-mom-z 'm) (F->C s->r))
   (->local 't 
            (coordinate-tuple 'r 'theta 'phi)
            (velocity-tuple 'rdot 'thetadot 'phidot))))
(* m (expt r 2) phidot (expt (sin theta) 2))

(show-expression
 ((Lagrangian->energy
   (L3-central 'm (literal-function 'V)))
  (->local 't
           (coordinate-tuple 'r 'theta 'phi)
           (velocity-tuple 'rdot 'thetadot 'phidot))))
(+ (* 1/2 m (expt r 2) (expt phidot 2) (expt (sin theta) 2))
   (* 1/2 m (expt r 2) (expt thetadot 2))
   (* 1/2 m (expt rdot 2))
   (V r))
|#

;;; Noether Theorem Support

(define ((Rx angle) q)
  (let ((ca (cos angle))
	(sa (sin angle)))
    (let ((x (ref q 0))
	  (y (ref q 1))
	  (z (ref q 2)))
      (up
       x
       (- (* ca y) (* sa z))
       (+ (* ca z) (* sa y))))))

(define ((Ry angle) q)
  (let ((ca (cos angle))
	(sa (sin angle)))
    (let ((x (ref q 0))
	  (y (ref q 1))
	  (z (ref q 2)))
      (up
       (+ (* ca x) (* sa z))
       y
       (- (* ca z) (* sa x))))))

(define ((Rz angle) q)
  (let ((ca (cos angle))
	(sa (sin angle)))
    (let ((x (ref q 0))
	  (y (ref q 1))
	  (z (ref q 2)))
      (up
       (- (* ca x) (* sa y))
       (+ (* ca y) (* sa x))
       z))))

;;; F-tilde is a parametric coordinate transformation that given
;;; parameters takes a state and returns transformed coordinates.   
;;; F-tilde may take an arbitrary number of real-valued parameters.
;;; F-tilde applied to zeros is the coordinate selector:  It takes a
;;; state and returns the coordinates.  The hypothesis of Noether's
;;; theorem is that the Lagrangian is invariant under the
;;; transformation for all values of the parameters.

;;; (D (lambda parms (compose L (F->C (apply F-tilde parms))))) = 0

(define (Noether-integral L F-tilde)
  (let ((zero-parameters (make-list (car (arity F-tilde)) 0)))
    (* ((partial 2) L) (apply (D F-tilde) zero-parameters))))

#|

(define (F-tilde theta phi psi)
  (compose (Rx theta)
	   (Ry phi)
	   (Rz psi)
	   coordinate))

(pe ((Noether-integral
      (L-central-rectangular 'm (literal-function 'Vr))
      F-tilde)
     (up 't
	 (up 'x 'y 'z)
	 (up 'vx 'vy 'vz))))
(down (+ (* -1 m vy z) (* m vz y))
      (+ (* m vx z) (* -1 m vz x))
      (+ (* -1 m vx y) (* m vy x)))

|#


;;;; Chapter 2

(define (antisymmetric->column-matrix A)
  ;; Should check for antisymmetry of A.
  (column-matrix (matrix-ref A 2 1)
		 (matrix-ref A 0 2)
		 (matrix-ref A 1 0)))

(define (3vector-components->antisymmetric v)
  (matrix-by-rows
   (list 0 (- (ref v 2)) (ref v 1))
   (list (ref v 2) 0 (- (ref v 0)))
   (list (- (ref v 1)) (ref v 0) 0)))

(define (((M-of-q->omega-of-t M-of-q) q) t)
  (define M-on-path (compose M-of-q q))
  (define (omega-cross t)
    (* ((D M-on-path) t)
       (m:transpose (M-on-path t))))
  (antisymmetric->column-matrix (omega-cross t)))

(define (((M-of-q->omega-body-of-t M-of-q) q) t)
  (* (m:transpose (M-of-q (q t)))
     (((M-of-q->omega-of-t M-of-q) q) t)))

(define (M->omega M-of-q)
  (Gamma-bar (M-of-q->omega-of-t M-of-q)))

(define (M->omega-body M-of-q)
  (Gamma-bar (M-of-q->omega-body-of-t M-of-q)))

(define (rotate-z-matrix angle)
  (matrix-by-rows
    (list (cos angle) (- (sin angle))               0)
    (list (sin angle)     (cos angle)               0)
    (list           0               0               1)))

(define (rotate-x-matrix angle)
  (matrix-by-rows 
    (list           1               0               0)
    (list           0     (cos angle) (- (sin angle)))
    (list           0     (sin angle)     (cos angle))))

(define (Euler->M angles)
  (let ((theta (ref angles 0))
	(phi (ref angles 1))
	(psi (ref angles 2)))
    (* (rotate-z-matrix phi)
       (* (rotate-x-matrix theta)
	  (rotate-z-matrix psi)))))

(define ((Euler->omega angles-path) t)
  (define (M-on-path t)
    (Euler->M (angles-path t)))
  (define (w-cross t)
    (* ((D M-on-path) t)
       (m:transpose (M-on-path t))))
  (antisymmetric->column-matrix (w-cross t)))

(define ((Euler->omega-body angles-path) t)
  (* (m:transpose (Euler->M (angles-path t)))
     ((Euler->omega angles-path) t)))

#|
(show-expression
  ((Euler->omega-body
    (coordinate-tuple (literal-function 'theta)
		      (literal-function 'phi)
		      (literal-function 'psi)))
   't))
(matrix-by-rows
 (list
  (+ (* ((D phi) t) (sin (theta t)) (sin (psi t)))
     (* (cos (psi t)) ((D theta) t))))
 (list
  (+ (* ((D phi) t) (sin (theta t)) (cos (psi t)))
     (* -1 (sin (psi t)) ((D theta) t))))
 (list (+ (* (cos (theta t)) ((D phi) t)) ((D psi) t))))
|#

#|
(show-expression
 (((M-of-q->omega-body-of-t Euler->M)
   (coordinate-tuple (literal-function 'theta)
		     (literal-function 'phi)
		     (literal-function 'psi)))
  't))
(matrix-by-rows
 (list
  (+ (* ((D phi) t) (sin (theta t)) (sin (psi t)))
     (* (cos (psi t)) ((D theta) t))))
 (list
  (+ (* ((D phi) t) (sin (theta t)) (cos (psi t)))
     (* -1 (sin (psi t)) ((D theta) t))))
 (list (+ (* (cos (theta t)) ((D phi) t)) ((D psi) t))))

(show-expression
 ((M->omega-body Euler->M)
  (->local 't 
           (coordinate-tuple 'theta 'phi 'psi)
           (velocity-tuple 'thetadot 'phidot 'psidot))))
(matrix-by-rows
 (list (+ (* phidot (sin psi) (sin theta)) (* thetadot (cos psi))))
 (list (+ (* phidot (sin theta) (cos psi)) (* -1 thetadot (sin psi))))
 (list (+ (* phidot (cos theta)) psidot)))
|#

(define (Euler-state->omega-body local)
  (let ((t (time local))
        (q (coordinate local))
        (qdot (velocity local)))
    (let ((theta (ref q 0))
          (psi (ref q 2))
          (thetadot (ref qdot 0))
          (phidot (ref qdot 1))
          (psidot (ref qdot 2)))
      (let ((omega-a (+ (* thetadot (cos psi))
                        (* phidot (sin theta) (sin psi))))
            (omega-b (+ (* -1 thetadot (sin psi))
                        (* phidot (sin theta) (cos psi))))
            (omega-c (+ (* phidot (cos theta)) psidot)))
        (column-matrix omega-a omega-b omega-c)))))

(define ((T-rigid-body A B C) local)
  (let ((omega-body (Euler-state->omega-body local)))
    (* 1/2
       (+ (* A (square (ref omega-body 0)))
	  (* B (square (ref omega-body 1)))
	  (* C (square (ref omega-body 2)))))))

(define ((Euler-state->L-body A B C) local)
  (let ((omega-body (Euler-state->omega-body local)))
    (column-matrix (* A (ref omega-body 0))
		   (* B (ref omega-body 1))
		   (* C (ref omega-body 2)))))

(define ((Euler-state->L-space A B C) local)
  (let ((angles (coordinate local)))
    (* (Euler->M angles)
       ((Euler-state->L-body A B C) local))))

#|
(define an-Euler-state
  (->local 't
           (coordinate-tuple 'theta 'phi 'psi)
           (velocity-tuple 'thetadot 'phidot 'psidot)))

(show-expression
 (ref
   (((partial 2) (T-rigid-body 'A 'B 'C))
    an-Euler-state)
   1))
(+ (* A phidot (expt (sin theta) 2) (expt (sin psi) 2))
   (* B phidot (expt (cos psi) 2) (expt (sin theta) 2))
   (* A thetadot (cos psi) (sin theta) (sin psi))
   (* -1 B thetadot (cos psi) (sin theta) (sin psi))
   (* C phidot (expt (cos theta) 2))
   (* C psidot (cos theta)))

(print-expression
 (- (ref ((Euler-state->L-space 'A 'B 'C) an-Euler-state) 2)        ;$L_z$
    (ref (((partial 2) (T-rigid-body 'A 'B 'C)) an-Euler-state) 1)  ;$p_\phi$
    ))
0

(show-expression
 (determinant
  (((compose (partial 2) (partial 2)) 
    (T-rigid-body 'A 'B 'C))
   an-Euler-state)))
(* A B C (expt (sin theta) 2))
|#

(define (relative-error value reference-value)
  (if (zero? reference-value)
      (error "Zero reference value -- RELATIVE-ERROR")
      (/ (- value reference-value) reference-value)))

#|
(define (rigid-sysder A B C)
  (Lagrangian->state-derivative (T-rigid-body A B C)))

(define ((monitor-errors win A B C L0 E0) state)
  (let ((t (time state))
	(L ((Euler-state->L-space A B C) state))
	(E ((T-rigid-body A B C) state)))
    (plot-point win t (relative-error (ref L 0) (ref L0 0)))
    (plot-point win t (relative-error (ref L 1) (ref L0 1)))
    (plot-point win t (relative-error (ref L 2) (ref L0 2)))
    (plot-point win t (relative-error E E0))))


;;; rkqc4
;;;(set! *ode-integration-method* 'qcrk4)
(define win (frame 0. 100. -1.e-12 1.e-12))

;;; bulirsch-stoer
;;;(set! *ode-integration-method* 'bulirsch-stoer)
(define win (frame 0. 100. -2.e-13 2.e-13))

(let ((A 1.) (B (sqrt 2.)) (C 2.)
      (state0 (up 0.0
		  (up 1. 0. 0.)
		  (up 0.1 0.1 0.1))))
  (let ((L0 ((Euler-state->L-space A B C) state0))
	(E0 ((T-rigid-body A B C) state0)))
    ((evolve rigid-sysder A B C)
     state0
     (monitor-errors win A B C L0 E0)
     0.1
     100.0
     1.0e-12)))

(graphics-close win)
(graphics-clear win)
|#

#|
(show-expression
 ((T-rigid-body 'A 'A 'C) 
   (->local 't 
            (coordinate-tuple 'theta 'phi 'psi)
            (velocity-tuple 'thetadot 'phidot 'psidot))))
(+ (* 1/2 A (expt phidot 2) (expt (sin theta) 2))
   (* 1/2 C (expt phidot 2) (expt (cos theta) 2))
   (* C phidot psidot (cos theta))
   (* 1/2 A (expt thetadot 2))
   (* 1/2 C (expt psidot 2)))

;;; Transformation of A(v):
;;;  M^T A(Mv) M = A(v) for arbitrary v orthogonal M

(print-expression
  (let ((Euler (coordinate-tuple 'theta 'phi 'psi))
	(v (coordinate-tuple 'x 'y 'z)))
    (let ((M (Euler->M Euler)))
      (- (* (3vector-components->antisymmetric (* M v))
	    M)
	 (* M
	    (3vector-components->antisymmetric v))))))
(matrix-by-rows (list 0 0 0) (list 0 0 0) (list 0 0 0))
|#

#| 
;;; Configuration equations for Euler's equations with Euler angles
   
(print-expression
  (let ((Euler (coordinate-tuple (literal-function 'theta)
				 (literal-function 'phi)
				 (literal-function 'psi))))
    (antisymmetric->column-matrix 
     (* (m:transpose ((Euler->M Euler) 't))
	((D (Euler->M Euler)) 't)))))
(matrix-by-rows
 (list
  (+ (* ((D phi) t) (sin (psi t)) (sin (theta t)))
     (* ((D theta) t) (cos (psi t)))))
 (list
  (+ (* ((D phi) t) (sin (theta t)) (cos (psi t)))
     (* -1 (sin (psi t)) ((D theta) t))))
 (list (+ (* (cos (theta t)) ((D phi) t)) ((D psi) t))))

(define foo
  (matrix-by-rows (list (cos 'psi)  (* (sin 'psi) (sin 'theta)) 0)
		  (list (- (sin 'psi)) (* (cos 'psi) (sin 'theta)) 0)
		  (list 0 (cos 'theta) 1)))


;;; funny -1 here in tex displayer
(show-expression (/ foo))
(matrix-by-rows
 (list (cos psi) (* -1 (sin psi)) 0)
 (list (/ (sin psi) (sin theta)) (/ (cos psi) (sin theta)) 0)
 (list (/ (* -1 (cos theta) (sin psi)) (sin theta))
       (/ (* -1 (cos theta) (cos psi)) (sin theta))
       1))

;;; do it directly

(define (foo Eulerdot)
  (let ((thetadot (ref Eulerdot 0))
	(phidot (ref Eulerdot 1))
	(psidot (ref Eulerdot 2)))
    (let ((Euler
	   (coordinate-tuple
	    (osculating-path (->L-state 't 'theta thetadot))
	    (osculating-path (->L-state 't 'phi phidot))
	    (osculating-path (->L-state 't 'psi psidot)))))
      (column-matrix->up
       (antisymmetric->column-matrix 
	(* (m:transpose ((Euler->M Euler) 't))
	   ((D (Euler->M Euler)) 't)))))))

(print-expression
 (/ 1 ((D foo) (coordinate-tuple 'x 'y 'z))))
(down
 (up (cos psi)
     (/ (sin psi) (sin theta))
     (/ (* -1 (cos theta) (sin psi)) (sin theta)))
 (up (* -1 (sin psi))
     (/ (cos psi) (sin theta))
     (/ (* -1 (cos theta) (cos psi)) (sin theta)))
 (up 0 0 1))

(print-expression
 ((D foo) (coordinate-tuple 'x 'y 'z)))
(down
 (up (cos psi)
     (* -1 (sin psi))
     0)
 (up (* (sin theta) (sin psi))
     (* (sin theta) (cos psi))
     (cos theta))
 (up 0 0 1))

(print-expression
  (/ (matrix-by-rows
      (list (cos 'psi) (* (sin 'theta) (sin 'psi)) 0)
      (list (* -1 (sin 'psi)) (* (sin 'theta) (cos 'psi)) 0)
      (list 0 (cos 'theta) 1))))
(matrix-by-rows
 (list (cos psi) (* -1 (sin psi)) 0)
 (list (/ (sin psi) (sin theta)) (/ (cos psi) (sin theta)) 0)
 (list (* -1 (cot theta) (sin psi)) (* -1 (cot theta) (cos psi)) 1))

|#

;;;; Chapter 3

(define ((Hamilton-equations Hamiltonian) q p)
  (let ((H-state-path (qp->H-state-path q p)))
    (- (D H-state-path)
       (compose (phase-space-derivative Hamiltonian)
                H-state-path))))

(define ((qp->H-state-path q p) t)
  (->H-state t (q t) (p t)))

(define ((phase-space-derivative Hamiltonian) H-state)
  (->H-state 1
             (((partial 2) Hamiltonian) H-state)
             (- (((partial 1) Hamiltonian) H-state))))


;;; If we express the energy in terms of t,Q,P we have the Hamiltonian.
;;; A Hamiltonian is an example of an H-function: an H-function takes
;;; 2 vector arguments and a scalar argument (t, Q, P).  It produces a
;;; scalar result.   

#|
(define ((H-rectangular m V) H-state)
  (let ((q (coordinate H-state))
        (p (momentum H-state)))
    (+ (/ (square p) (* 2 m))
       (V (ref q 0) (ref q 1)))))

(show-expression 
 (((Hamilton-equations
      (H-rectangular 
          'm
          (literal-function 'V (-> (X Real Real) Real))))
     (coordinate-tuple (literal-function 'x)
                       (literal-function 'y))
     (momentum-tuple (literal-function 'p_x)
                     (literal-function 'p_y)))
    't))
(up
 0
 (up (+ ((D x) t) (/ (* -1 (p_x t)) m))
     (+ ((D y) t) (/ (* -1 (p_y t)) m)))
 (down (+ ((D p_x) t) (((partial 0) V) (x t) (y t)))
       (+ ((D p_y) t) (((partial 1) V) (x t) (y t)))))
|#


(define (make-Hamiltonian kinetic-energy potential-energy)
  (+ kinetic-energy potential-energy))

;;; If we express the energy in terms of t,Q,P we have the Hamiltonian
;;;        H(t,Q,P) = P*Qdot - L(t, Q, Qdot(t, Q, P))
;;; To do this we need to invert P(t, Q, Qdot) to get Qdot(t, Q, P).
;;; This is easy when L is a quadratic form in Qdot:
;;;        L(t, Q, Qdot) = 1/2*Qdot*M*Qdot + B*Qdot - V
;;; Fortunately this is the case in almost all of Newtonian mechanics,
;;; otherwise the P(t,Q,Qdot) function would be much more difficult to
;;; invert to obtain Qdot(t,Q,P).

;;; Assume that F is quadratic in its arguments
;;;  F(u) = 1/2 A u u + b u + c
;;;  then v = A u + b, so u = A^(-1) (v - b)

#|
(define (Legendre-transform-procedure F)
  (let ((w-of-v (D F)))
    (define (G w)
      (let ((z (dual-zero w)))
        (let ((M ((D w-of-v) z))
              (b (w-of-v z)))
          (let ((v (/ (- w b) M)))
            (- (* w v) (F v))))))
    G))

(define (dual-zero v)
  (if (or (column? v) (row? v))
      (s:generate (s:length v) (s:opposite v)
		  (lambda (i) :zero))
      :zero))
|#

;;; A better definition of Legendre transform that works for
;;;   structured coordinates that have substructure

(define (Legendre-transform-procedure F)
  (let ((w-of-v (D F)))
    (define (G w)
      (let ((z (compatible-zero w)))
        (let ((M ((D w-of-v) z))
              (b (w-of-v z)))
          (let ((v (* (s:inverse z M z) (- w b))))
            (- (* w v) (F v))))))
    G))

(define Legendre-transform
  (make-operator Legendre-transform-procedure
                 'Legendre-transform))

;;; Notice that Lagrangians and Hamiltonians are symmetrical with
;;; respect to the Legendre transform.

(define ((Lagrangian->Hamiltonian-procedure the-Lagrangian) H-state)
  (let ((t (time H-state))
	(q (coordinate H-state))
	(p (momentum H-state)))
    (define (L qdot)
      (the-Lagrangian (->L-state t q qdot)))
    ((Legendre-transform-procedure L) p)))

(define Lagrangian->Hamiltonian
  (make-operator Lagrangian->Hamiltonian-procedure
                 'Lagrangian->Hamiltonian))


(define ((Hamiltonian->Lagrangian-procedure the-Hamiltonian) L-state)
  (let ((t (time L-state))
	(q (coordinate L-state))
	(qdot (velocity L-state)))
    (define (H p)
      (the-Hamiltonian (->H-state t q p)))
    ((Legendre-transform-procedure H) qdot)))

(define Hamiltonian->Lagrangian
  (make-operator Hamiltonian->Lagrangian-procedure
                 'Hamiltonian->Lagrangian))

(define ((Lstate->Hstate L) Ls)
  (up (time Ls)
      (coordinate Ls)
      (((partial 2) L) Ls)))

(define ((Hstate->Lstate H) Hs)
  (up (time Hs)
      (coordinate Hs)
      ((Poisson-bracket coordinate H) Hs)))

#|
(define ((L-rectangular m V) local)
  (let ((q (coordinate local))
        (qdot (velocity local)))
    (- (* 1/2 m (square qdot))
       (V (ref q 0) (ref q 1)))))

(show-expression 
 ((Lagrangian->Hamiltonian
   (L-rectangular 
      'm
      (literal-function 'V (-> (X Real Real) Real))))
  (->H-state 't
             (coordinate-tuple 'x 'y)
             (momentum-tuple 'p_x 'p_y))))
(+ (V x y)
   (/ (* 1/2 (expt p_x 2)) m)
   (/ (* 1/2 (expt p_y 2)) m))
|#

#|
(show-expression
 ((Lagrangian->Hamiltonian 
    (L-central-polar 'm (literal-function 'V)))
  (->H-state 't
             (coordinate-tuple 'r 'phi)
             (momentum-tuple 'p_r 'p_phi))))
(+ (V r)
   (/ (* 1/2 (expt p_r 2)) m)
   (/ (* 1/2 (expt p_phi 2)) (* m (expt r 2))))


(show-expression
 (((Hamilton-equations
     (Lagrangian->Hamiltonian 
       (L-central-polar 'm (literal-function 'V))))
   (coordinate-tuple (literal-function 'r)
                     (literal-function 'phi))
   (momentum-tuple (literal-function 'p_r)
                   (literal-function 'p_phi)))
  't))
(up
 0
 (up (+ ((D r) t) (/ (* -1 (p_r t)) m))
     (+ ((D phi) t) (/ (* -1 (p_phi t)) (* m (expt (r t) 2)))))
 (down
  (+ ((D p_r) t)
     ((D V) (r t))
     (/ (* -1 (expt (p_phi t) 2)) (* m (expt (r t) 3))))
  ((D p_phi) t)))

;;; If we substitute a Coulomb potential in for V we get the equations 
;;;  for satellite motion around a spherical primary.

(show-expression
 (((Hamilton-equations
    (Lagrangian->Hamiltonian
     (L-central-polar 'm
		      (lambda (r)
			(- (/ (* 'GM 'm) r))))))
   (coordinate-tuple (literal-function 'r)
		     (literal-function 'phi))
   (momentum-tuple
    (literal-function 'p_r)
    (literal-function 'p_phi)))
  't))
(up
 0
 (up (+ ((D r) t) (/ (* -1 (p_r t)) m))
     (+ ((D phi) t) (/ (* -1 (p_phi t)) (* m (expt (r t) 2)))))
 (down
  (+ ((D p_r) t)
     (/ (* GM m) (expt (r t) 2))
     (/ (* -1 (expt (p_phi t) 2)) (* m (expt (r t) 3))))
  ((D p_phi) t)))
|#

#|
;;; Continuing with our coupled harmonic oscillators
;;;  we obtain the Hamiltonian:

(show-expression
 ((Lagrangian->Hamiltonian
   (L-coupled-harmonic (down (down 'm_1 0)
			     (down 0 'm_2))
		       (down (down 'k_1 'c)
			     (down 'c 'k_2))))
  (->H-state 't
	      (coordinate-tuple 'x_1 'x_2)
	      (momentum-tuple 'p_1 'p_2))))
(+ (* c x_1 x_2)
   (* 1/2 k_1 (expt x_1 2))
   (* 1/2 k_2 (expt x_2 2))
   (/ (* 1/2 (expt p_2 2)) m_2)
   (/ (* 1/2 (expt p_1 2)) m_1))
|#

;;; phase-space-derivative returns #( 1 qdot pdot )

(define ((Hamilton-equation-1 H) q p)
  (let ((state-path (qp->H-state-path q p)))
    (- (D (compose state->q state-path))
       (compose ((partial 2) H) state-path))))

(define ((Hamilton-equation-2 H) q p)
  (let ((state-path (qp->H-state-path q p)))
    (+ (D (compose momentum state-path))
       (compose ((partial 1) H) state-path))))

#|
;;; This formulation, although clearer and prettier than the previous
;;; one forces a distinction between row and column vectors.

(define (((Hamilton-equation-1 H) q p) t)
  (- ((D q) t) 
     (((partial 2) H) (->H-state t (q t) (p t)))))

(define (((Hamilton-equation-2 H) q p) t)
  (+ ((D p) t) 
     (((partial 1) H) (->H-state t (q t) (p t)))))
|#

#|
(show-expression
 (((Hamilton-equation-1
    (Lagrangian->Hamiltonian
     (L-periodically-driven-pendulum 'm 'l 'g 'a 'omega)))
   (literal-function 'q)
   (literal-function 'p))
  't))
(+ ((D q) t)
   (/ (* -1 a omega (sin (* omega t)) (sin (q t))) l) 
   (/ (* -1 (p t)) (* (expt l 2) m)))

(show-expression
 (((Hamilton-equation-2
    (Lagrangian->Hamiltonian
     (L-periodically-driven-pendulum 'm 'l 'g 'a 'omega)))
   (literal-function 'q)
   (literal-function 'p))
  't))
(+
 (* (expt a 2)
    m
    (expt omega 2)
    (sin (q t))
    (cos (q t))
    (expt (sin (* omega t)) 2))
 (/ (* a omega (p t) (sin (* omega t)) (cos (q t))) l)
 (* g l m (sin (q t)))
 ((D p) t))
|#

#|
(define ((H-central-polar m V) state)
  (let ((q (coordinate state))
        (p (momentum state)))
    (let ((r ((component 0) q))
          (phi ((component 1) q))
          (pr ((component 0) p))
          (pphi ((component 1) p)))
      (+ (/ (+ (square pr)
	       (square (/ pphi r)))
	    (* 2 m))
         (V r)))))

(show-expression
 (((Hamilton-equation-1
    (H-central-polar 'm (literal-function 'V)))
   (coordinate-tuple (literal-function 'r) (literal-function 'phi))
   (momentum-tuple (literal-function 'p_r) (literal-function 'p_phi)))
  't))
(up
 (+ ((D r) t) (/ (* -1 (p_r t)) m))
 (+ ((D phi) t) (/ (* -1 (p_phi t)) (* m (expt (r t) 2)))))

(show-expression
 (((Hamilton-equation-2
    (H-central-polar 'm (literal-function 'V)))
   (coordinate-tuple (literal-function 'r) (literal-function 'phi))
   (momentum-tuple (literal-function 'p_r) (literal-function 'p_phi)))
  't))
(down
 (+ ((D p_r) t)
    ((D V) (r t))
    (/ (* -1 (expt (p_phi t) 2)) (* m (expt (r t) 3))))
 ((D p_phi) t))
|#

#|
(define (state->qp dynamic-state)
  (vector-tail dynamic-state 1))

(define ((Hamilton-equations Hamiltonian) q p)
  (let ((H-state-path (qp->H-state-path q p))
	(dH (phase-space-derivative Hamiltonian)))
    (compose state->qp
	     (- (D H-state-path)
		(compose dH H-state-path)))))
|#

#|
;;; Going on to obtain Hamilton's equations here:

(show-expression
 (((Hamilton-equations
    (Lagrangian->Hamiltonian
     (L-harmonic 'm 'k)))
   (coordinate-tuple
    (literal-function 'x_1)
    (literal-function 'x_2))
   (momentum-tuple
    (literal-function 'p_1)
    (literal-function 'p_2)))
  't))
(up
 0
 (up (+ ((D x_1) t) (/ (* -1 (p_1 t)) m_1))
     (+ ((D x_2) t) (/ (* -1 (p_2 t)) m_2)))
 (down (+ (* c (x_2 t)) (* k_1 (x_1 t)) ((D p_1) t))
       (+ (* c (x_1 t)) (* k_2 (x_2 t)) ((D p_2) t))))
|#

#|
;;; Continuing our coupled oscillators example:

(show-expression
 (((Hamilton-equations
    (Lagrangian->Hamiltonian
     (L-coupled-harmonic (down (down 'm_1 0)
			       (down 0 'm_2))
			 (down (down 'k_1 'c)
			       (down 'c 'k_2)))))
   (coordinate-tuple (literal-function 'x_1)
		     (literal-function 'x_2))
   (momentum-tuple (literal-function 'p_1)
		   (literal-function 'p_2)))
  't))
(up
 0
 (up (+ ((D x_1) t) (/ (* -1 (p_1 t)) m_1))
     (+ ((D x_2) t) (/ (* -1 (p_2 t)) m_2)))
 (down (+ (* c (x_2 t)) (* k_1 (x_1 t)) ((D p_1) t))
       (+ (* c (x_1 t)) (* k_2 (x_2 t)) ((D p_2) t))))
|#

#|
;;; Continuation of demonstration of bundled coordinates.

(se
 (((Hamilton-equations (Lagrangian->Hamiltonian (L-two-particle 'm_1 'm_2)))
   (coordinate-tuple 
    (coordinate-tuple (literal-function 'x_1) (literal-function 'y_1))
    (coordinate-tuple (literal-function 'x_2) (literal-function 'y_2)))
   (momentum-tuple 
    (momentum-tuple (literal-function 'p_x_1) (literal-function 'p_y_1))
    (momentum-tuple (literal-function 'p_x_2) (literal-function 'p_y_2))))
  't))

(up
 0
 (up
  (up (+ ((D x_1) t) (/ (* -1 (px_1 t)) m_1))
      (+ ((D y_1) t) (/ (* -1 (py_1 t)) m_1)))
  (up (+ ((D x_2) t) (/ (* -1 (px_2 t)) m_2))
      (+ ((D y_2) t) (/ (* -1 (py_2 t)) m_2))))
 (down
  (down
   (+ ((D px_1) t)
      (((partial 0 0) V) (up (x_1 t) (y_1 t)) (up (x_2 t) (y_2 t))))
   (+ ((D py_1) t)
      (((partial 0 1) V) (up (x_1 t) (y_1 t)) (up (x_2 t) (y_2 t)))))
  (down
   (+ ((D px_2) t)
      (((partial 1 0) V) (up (x_1 t) (y_1 t)) (up (x_2 t) (y_2 t))))
   (+ ((D py_2) t)
      (((partial 1 1) V) (up (x_1 t) (y_1 t)) (up (x_2 t) (y_2 t)))))))
|#

;;; The Poisson Bracket is a differential operator on H-functions:

(define (Poisson-bracket f g)
  (if (structure? f)
      (if (structure? g)
	  (error "Poisson bracket of two structures" f g)
	  (s:generate (s:length f) (s:same f)
		      (lambda (i)
			(Poisson-bracket (ref f i) g))))
      (if (structure? g)
	  (s:generate (s:length g) (s:same g)
		      (lambda (i)
			(Poisson-bracket f (ref g i))))
	  (- (* ((partial 1) f) ((partial 2) g))
	     (* ((partial 2) f) ((partial 1) g))))))

#| 
;;; This can give WRONG ANSWERS if handed structure valued functions...

(define (Poisson-bracket f g)
  (- (* ((partial 1) f) ((partial 2) g))
     (* ((partial 2) f) ((partial 1) g))))

;;; for example L = [Lx, Ly, Lz] 
;;; where Li are components of angular momentum


|#

#|

(define F (literal-function 'F (Hamiltonian 2)))
(define G (literal-function 'G (Hamiltonian 2)))
(define H (literal-function 'H (Hamiltonian 2)))


;;; Jacobi identity

(pe ((+ (Poisson-bracket F (Poisson-bracket G H))
	(Poisson-bracket G (Poisson-bracket H F))
	(Poisson-bracket H (Poisson-bracket F G)))
     (up 't (up 'x 'y) (down 'px 'py))))
0

;;; for Lie derivatives

(define L_F (Lie-derivative F))
(define L_G (Lie-derivative G))

(pe (((+ (commutator L_F L_G)
	 (Lie-derivative (Poisson-bracket F G)))
      H)
     (up 't (up 'x 'y) (down 'px 'py))))
0

|#


(define Sx (compose (component 0) coordinate))
(define Sy (compose (component 1) coordinate))
(define Sz (compose (component 2) coordinate))

(define Spx (compose (component 0) momentum))
(define Spy (compose (component 1) momentum))
(define Spz (compose (component 2) momentum))

(define Lx (- (* Sy Spz) (* Spy Sz)))
(define Ly (- (* Sz Spx) (* Spz Sx)))
(define Lz (- (* Sx Spy) (* Spx Sy)))

#|

(define L (down Lx Ly Lz))

(define 3-state 
  (->H-state 't 
	     (coordinate-tuple 'x 'y 'z)
	     (momentum-tuple 'p_x 'p_y 'p_z)))

(pe ((Poisson-bracket Lx L) 3-state))
(down 0 (+ (* -1 p_x y) (* p_y x)) (+ (* -1 p_x z) (* p_z x)))
|#

;;; Actually, this is not the vector field, because a vector field
;;; takes a function as an argument, and here the function is given
;;; as the identity.

;;;(define (Hamiltonian-vector-field H)
;;;  (Poisson-bracket identity H))

;;; name conflicts with cartesian product
;;; (define X Hamiltonian-vector-field) 


#|
;;; From 3.3 -- Phase-space reduction


(define ((L-axisymmetric-top A C gMR) local)
  (let ((q (coordinate local))
        (qdot (velocity local)))
    (let ((theta (ref q 0))
          (thetadot (ref qdot 0))
          (phidot (ref qdot 1))
          (psidot (ref qdot 2)))
      (+ (* 1/2 A
            (+ (square thetadot)
               (square (* phidot (sin theta)))))
         (* 1/2 C
            (square (+ psidot (* phidot (cos theta)))))
         (* -1 gMR (cos theta))))))


(show-expression
 ((Lagrangian->Hamiltonian (L-axisymmetric-top 'A 'C 'gMR)) 
  (->H-state 't
             (vector 'theta 'phi 'psi)
             (vector 'p_theta 'p_phi 'p_psi))))
(+ (* gMR (cos theta))
   (/ (* 1/2 (expt p_psi 2)) C)
   (/ (* 1/2 (expt p_psi 2) (expt (cos theta) 2)) (* A (expt (sin theta) 2)))
   (/ (* 1/2 (expt p_theta 2)) A)
   (/ (* -1 p_phi p_psi (cos theta)) (* A (expt (sin theta) 2)))
   (/ (* 1/2 (expt p_phi 2)) (* A (expt (sin theta) 2))))
|#

#|
;;; ********** This got ugly *************

(show-expression
 (((Hamilton-equations
    (Lagrangian->Hamiltonian
     (L-axisymmetric-top 'A 'C 'gMR)))
   (coordinate-tuple (literal-function 'theta)
		     (literal-function 'phi)
		     (literal-function 'psi))
   (momentum-tuple (literal-function 'p_theta)
		   (literal-function 'p_phi)
		   (literal-function 'p_psi)))
  't))
(up
 0
 (up
  (+ ((D theta) t) (/ (* -1 (p_theta t)) A))
  (+ ((D phi) t)
     (/ (* (cos (theta t)) (p_psi t)) (* A (expt (sin (theta t)) 2)))
     (/ (* -1 (p_phi t)) (* A (expt (sin (theta t)) 2))))
  (+
   ((D psi) t)
   (/ (* -1 (p_psi t)) C)
   (/ (* -1 (expt (cos (theta t)) 2) (p_psi t)) (* A (expt (sin (theta t)) 2)))
   (/ (* (p_phi t) (cos (theta t))) (* A (expt (sin (theta t)) 2)))))
 (down
  (+
   (/ (* -1 gMR (expt (cos (theta t)) 4)) (expt (sin (theta t)) 3))
   ((D p_theta) t)
   (/ (* 2 gMR (expt (cos (theta t)) 2)) (expt (sin (theta t)) 3))
   (/ (* (p_phi t) (expt (cos (theta t)) 2) (p_psi t))
      (* A (expt (sin (theta t)) 3)))
   (/ (* -1 (cos (theta t)) (expt (p_psi t) 2)) (* A (expt (sin (theta t)) 3)))
   (/ (* -1 (expt (p_phi t) 2) (cos (theta t))) (* A (expt (sin (theta t)) 3)))
   (/ (* -1 gMR) (expt (sin (theta t)) 3))
   (/ (* (p_phi t) (p_psi t)) (* A (expt (sin (theta t)) 3))))
  ((D p_phi) t)
  ((D p_psi) t)))
|#

#|
;;; Substituting p for p_psi and p_phi I do not get 3.87!
;;;  But, using (->poisson-form '(square (tan (/ theta 2))))
;;;  solves my problem.
(pp
 ((compose simplify expression)
  (substitute
   'p 'p_psi
   (substitute
    'p 'p_phi
    ((Lagrangian->Hamiltonian
      (L-axisymmetric-top 'A 'C 'gMR)) 
     (->H-state 't
		 (coordinate-tuple 'theta 'phi 'psi)
		 (momentum-tuple 'p_theta 'p_phi 'p_psi)))))))
(+ (/ (* gMR (expt (cos theta) 2)) (+ 1 (cos theta)))
   (/ (* gMR (cos theta)) (+ 1 (cos theta)))
   (/ (* 1/2 (expt p 2) (cos theta)) (+ (* C (cos theta)) C))
   (/ (* 1/2 (expt p 2)) (+ (* C (cos theta)) C))
   (/ (* -1/2 (expt p 2) (cos theta)) (+ (* A (cos theta)) A))
   (/ (* 1/2 (expt p 2)) (+ (* A (cos theta)) A))
   (/ (* 1/2 (expt p_theta 2) (cos theta)) (+ (* A (cos theta)) A))
   (/ (* 1/2 (expt p_theta 2)) (+ (* A (cos theta)) A)))
|#

#|
(define ((ueff p A C gMR) theta)
  (+ (/ (square p) (* 2 C))
     (* (/ (square p) (* 2 A))
	(square (tan (/ theta 2))))
     (* gMR (cos theta))))

|#
#|
;;; Critical value of bifurcation when D^2 Ueff (0) = 0

(print-expression
 (((square derivative) (ueff 'p_c 'A 'C 'gMR)) 0))
(+ (* -1 gMR) (/ (* 1/4 (expt p_c 2)) A))

;;; critical angular speed in RPM is:
(* (/ 60 2pi) (/ 7.734804457773965e-3 6.6e-5))
;Value: 1119.1203302763215
|#

(define Hamiltonian->state-derivative phase-space-derivative)

#|
(show-expression
 ((Lagrangian->Hamiltonian
    (L-periodically-driven-pendulum 'm 'l 'g 'a 'omega))
  (->H-state 't 'theta 'p_theta)))
(+
 (* -1/2
    (expt a 2)
    m
    (expt omega 2)
    (expt (cos theta) 2)
    (expt (sin (* omega t)) 2))
 (* a g m (cos (* omega t)))
 (/ (* a omega p_theta (sin theta) (sin (* omega t))) l)
 (* -1 g l m (cos theta))
 (/ (* 1/2 (expt p_theta 2)) (* (expt l 2) m)))

(define (H-pend-sysder m l g a omega)
  (Hamiltonian->state-derivative
    (Lagrangian->Hamiltonian
      (L-periodically-driven-pendulum m l g a omega))))

;;; for driven-pendulum-phase-space

(define window (frame -pi pi -10.0 10.0))

(start-gnuplot "pendulum-2x")

(define ((monitor-p-theta win) state)
  (let ((q ((principal-value pi) (coordinate state)))
	(p (momentum state)))
    (plot-point win q p)))
  
(let ((m 1.)		                ;m=1kg
      (l 1.)				;l=1m
      (g 9.8)				;g=9.8m/s$^2$
      (A 0.1)				;A=1/10 m
      (omega (* 2 (sqrt 9.8))))
  ((evolve H-pend-sysder m l g A omega)
   (up 0.0				;t$_0$=0
       1.0				;theta$_0$=1 radian
       0.0)				;thetadot$_0$=0 radians/s
   (monitor-p-theta window)
   0.01					;step between plotted points
   100.0				;final time
   1.0e-12))

(stop-gnuplot)

(graphics-clear window)

;;; for driven-pend-nuniq1 and 2
(let ((m 1.)		                ;m=1kg
      (l 1.)				;l=1m
      (g 9.8)				;g=9.8m/s$^2$
      (A 0.1)				;A=1/10 m
      (omega (* 2 (sqrt 9.8))))
  ((evolve H-pend-sysder m l g A omega)
   (up 0.0				;t$_0$=0
       1.0				;theta$_0$=1 radian
       0.0)				;thetadot$_0$=0 radians/s
   (monitor-p-theta window)
   0.01					;step between plotted points
   10.0					;final time
   1.0e-12))

(define ((monitor-pprime-theta mlA omega win) state)
  (let ((t (time state))
	(q ((principal-value pi) (coordinate state)))
	(p (momentum state)))
    (plot-point 
     win
     q 
     (+ p (* mlA omega (sin q) (sin (* omega t)))))))

(let ((m 1.)		                ;m=1kg
      (l 1.)				;l=1m
      (g 9.8)				;g=9.8m/s$^2$
      (A 0.1)				;A=1/10 m
      (omega (* 2 (sqrt 9.8))))
  ((evolve H-pend-sysder m l g A omega)
   (up 0.0				;t$_0$=0
       1.0				;theta$_0$=1 radian
       0.0)				;thetadot$_0$=0 radians/s
   (monitor-pprime-theta (* m l A) omega window)
   0.01					;step between plotted points
   10.0					;final time
   1.0e-12))

(graphics-close window)
|#

;;; Now to do the Poincare section.  
;;;  Map explorer:
;;;   Left button starts a trajectory.
;;;   Middle button continues a trajectory.
;;;   Right button interrogates coordinates.

(define (explore-map window poincare-map #!optional n)
  (define (iterate-map i x y)
    (if (fix:> i 0) 
	(poincare-map 
	 x y
	 (lambda (nx ny)
	   (plot-point window x y)
	   (iterate-map (fix:- i 1) nx ny))
	 (lambda ()  
	   (newline)
	   (display "Illegal point: ")
	   (write (list  x y))
	   (button-loop x y)))
	(button-loop x y)))
  (define (button-loop ox oy)
    (pointer-coordinates 
     window 
     (lambda (x y button)
       (case button
	 ((0)
	  (write-line (list x y))
	  (display " started.")
	  (iterate-map n x y))
	 ((1) 
	  (write-line (list ox oy))
	  (display " continued.")
	  (iterate-map n ox oy))
	 ((2)
	  (write-line (list x y))
	  (display " hit.")
	  (button-loop ox oy))))))
  (if (default-object? n) (set! n 1000))
  (newline)
  (display "Left button starts a trajectory.")
  (newline)
  (display "Middle button continues a trajectory.")
  (newline)
  (display "Right button interrogates coordinates.")
  (button-loop 9. 9.))

(define (pointer-coordinates window continue)
  (beep)
  (get-pointer-coordinates window continue))
  
;;; For iterating maps of the shape used by explore-map.

(define ((iterated-map map n) x y continue fail)
  (if (fix:< n 0) (error "iterated-map: cannot invert map"))
  (let loop ((x x) (y y) (i n))
    (if (fix:= i 0) 
	(continue x y)
	(map x y
	     (lambda (nx ny)
	       (loop nx ny (fix:- i 1)))
	     fail))))


(define (display-map window poincare-map x y n)
  (plot-point window x y)
  (if (fix:> n 0)
      (poincare-map 
       x y
       (lambda (nx ny)
	 (display-map window poincare-map nx ny (- n 1)))
       (lambda ()  
	   (newline)
	   (display "Illegal point: ")
	   (write (list  x y))))))


;;; For example, first a simple problem.

(define ((standard-map K) x y continue fail)
  (let ((yp (flo:pv (flo:+ y (flo:* K (flo:sin x))))))
    (continue (flo:pv (flo:+ x yp)) yp)))

;;; This is the 0-2pi principal value:

(define (flo:pv x)
  (flo:- x (flo:* 2pi (flo:floor (flo:/ x 2pi)))))

#|
(define win (frame 0.0 2pi 0.0 2pi))
(explore-map win (standard-map 1.0) 1000)
(graphics-close win)
|#

#| 
;;; The driven pendulum map is very interesting.
(set! *ode-integration-method* 'bulirsch-stoer)

(define (driven-pendulum-map mass l g a omega)
  (let ((map-period (/ 2pi omega)))
    (lambda (theta p-theta continue fail)
      (let ((ns
	     ((state-advancer H-pend-sysder
			      mass l g a omega)
	      (up			;state to start from
	       0.0			;t0 = phase/drive-freq.
	       theta
	       p-theta)
	      map-period			
	      1.0e-10)))
	(continue 
	 ((principal-value pi) (vector-ref ns 1))
	 (vector-ref ns 2))))))

;;; driven at twice the small-oscillation resonance.
;;; driven-pend-fig6.ps
(define win (frame -pi pi -10.0 10.0))

(let* ((m 1.0)				;m=1kg
       (l 1.0)				;l=1m
       (g 9.8)				;g=9.8m/s^2
       (small-amplitude-frequency (sqrt (/ g l)))
       (drive-amplitude 0.1)		;a=1/10 m
       (drive-frequency
	(* 2.0 small-amplitude-frequency)))
  (explore-map
   win
   (driven-pendulum-map m l g drive-amplitude drive-frequency)
   1000))

(graphics-close win)

;;; driven off resonance, at 4.2 times the natural frequency.
;;; driven-pend-fig7.ps

(define win (frame -pi pi -20 20))

(let* ((m 1.0)				;m=1kg
       (l 1.0)				;l=1m
       (g 9.8)				;g=9.8m/s^2
       (small-amplitude-frequency (sqrt (/ g l)))
       (drive-amplitude 0.05)		;a=1/20 m
       (drive-frequency
	(* 4.2 small-amplitude-frequency)))
  (explore-map
   win
   (driven-pendulum-map m l g drive-amplitude drive-frequency)
   1000))

(graphics-close win)

;;; Driven at high frequency, with larger amplitude.
;;; This shows the stable inverted pendulum island.
;;; driven-pend-fig8.ps

(define win (frame -pi pi -20 20))

(let* ((m 1.0)				;m=1kg
       (l 1.0)				;l=1m
       (g 9.8)				;g=9.8m/s^2
       (small-amplitude-frequency (sqrt (/ g l)))
       (drive-amplitude 0.2)		;a=1/20 m
       (drive-frequency
	(* 10.1 small-amplitude-frequency)))
  (explore-map
   win
   (driven-pendulum-map m l g drive-amplitude drive-frequency)
   1000))

(graphics-close win)
|#

#|
;;; An alternative, using construction and selection rather than
;;; continuations Here a map returns a new vector of x,y or an object
;;; describing the reason why it cannot.

(define (explore-map window poincare-map #!optional n)
  (define (iterate-map i x y)
    (if (fix:> i 0) 
	(let ((nxy (poincare-map x y)))
	  (if (vector? nxy)
	      (begin (plot-point window x y)
		     (iterate-map (fix:- i 1)
				  (vector-ref nxy 0)
				  (vector-ref nxy 1)))
	      (begin (newline)
		     (display "Illegal point: ")
		     (write (list  x y))
		     (button-loop x y))))
	(button-loop x y)))
  (define (button-loop ox oy)
    (pointer-coordinates 
     window 
     (lambda (x y button)
       (case button
	 ((0)
	  (write-line (list x y))
	  (display " started.")
	  (iterate-map n x y))
	 ((1) 
	  (write-line (list ox oy))
	  (display " continued.")
	  (iterate-map n ox oy))
	 ((2)
	  (write-line (list x y))
	  (display " hit.")
	  (button-loop ox oy))))))
  (if (default-object? n) (set! n 1000))
  (newline)
  (display "Left button starts a trajectory.")
  (newline)
  (display "Middle button continues a trajectory.")
  (newline)
  (display "Right button interrogates coordinates.")
  (button-loop 9. 9.))

;;; In this form the map is slightly different, but everything 
;;; else is the same.

(define ((standard-map K) x y)
  (let ((yp (flo:pv (flo:+ y (flo:* K (flo:sin x))))))
    (vector (flo:pv (flo:+ x yp)) yp)))

(define (driven-pendulum-map mass l g a omega)
  (let ((map-period (/ 2pi omega)))
    (lambda (theta p-theta)
      (let ((ns
	     ((state-advancer H-pend-sysder
			      mass l g a omega)
	      (up			;state to start from
	       0.0			;t0 = phase/drive-freq.
	       theta
	       p-theta)
	      map-period			
	      1.0e-10)))
	(vector
	 ((principal-value pi) (vector-ref ns 1))
	 (vector-ref ns 2))))))

|#


;;;; Chapter 4

(define ((unstable-manifold T xe ye dx dy A eps) param)
  (let ((n (floor->exact (/ (log (/ param eps)) (log A)))))
    ((iterated-map T n) (+ xe (* dx (/ param (expt A n))))
                        (+ ye (* dy (/ param (expt A n))))
                        cons
			(lambda ()
			  (error "Failed")))))

#| in open.scm

(define (plot-parametric-fill win f a b near?)
  (let loop ((a a) (xa (f a)) (b b) (xb (f b)))
    (let ((m (/ (+ a b) 2)))
      (let ((xm (f m)))
        (plot-point win (car xm) (cdr xm))
        (if (not (and (near? xa xm) (near? xb xm)))
            (begin (loop a xa m xm)
                   (loop m xm b xb)))))))

(define (cylinder-near? eps)
  (let ((eps2 (square eps)))
    (lambda (x y)
      (< (+ (square ((principal-value pi)
                     (- (car x) (car y))))
            (square (- (cdr x) (cdr y))))
         eps2))))

|#

(define (radially-mapping-points map Jmin Jmax phi eps)
  (bisect 
    (lambda (J) 
      ((principal-value pi)
       (- phi (map phi J (lambda (phip Jp) phip) list))))
    Jmin Jmax eps))


(define (find-invariant-curve map rn theta0 Jmin Jmax eps)
  (bisect (lambda (J) (which-way? rn theta0 J map))
          Jmin Jmax eps))

#|

(define (which-way? rn theta0 J0 map)
  (compare-streams
   (position-stream theta0
                    (orbit-stream map theta0 J0)
                    '())
   (position-stream theta0
                    (orbit-stream (circle-map rn) theta0 J0)
                    '())
   0))

(define (circle-map rotation-number)
  (let ((delta-theta (* :2pi rotation-number)))
    (lambda (theta y result fail)
      (result ((principal-value :2pi) (+ theta delta-theta))
              y))))

(define (orbit-stream the-map x y)
  (cons-stream (list x y)
               (the-map x y 
                       (lambda (nx ny)
                         (orbit-stream the-map nx ny))
                       (lambda () 'fail))))

(define (position-stream cut orbit list)
  (insert! ((principal-value cut) (car (head orbit)))
           list
           (lambda (nlist position)
             (cons-stream 
               position
               (position-stream cut (tail orbit) nlist)))))


(define (insert! x set cont)
  (cond ((null? set)
         (cont (list x) 1))
        ((< x (car set))
         (cont (cons x set) 0))
        (else
         (let lp ((i 1) (lst set))
           (cond ((null? (cdr lst))
                  (set-cdr! lst (cons x (cdr lst)))
                  (cont set i))
                 ((< x (cadr lst))
                  (set-cdr! lst (cons x (cdr lst)))
                  (cont set i))
                 (else
                  (lp (+ i 1) (cdr lst))))))))

(define (compare-streams s1 s2 count)
  (if (= (head s1) (head s2))
      (compare-streams (tail s1) (tail s2) (+ count 1))
      ((principal-range count) (- (head s2) (head s1)))))

(find-invariant-curve (standard-map 0.95)
                      (- 1 (/ 1 golden-mean))
                      0.0
                      2.0
                      2.2
                      1e-5)
;Value: 2.114462280273437

|#

(define (which-way? rotation-number x0 y0 map)
  (let ((pv (principal-value (+ x0 pi))))
    (let lp ((n 0) 
             (z x0) (zmin (- x0 2pi)) (zmax (+ x0 2pi))
             (x x0) (xmin (- x0 2pi)) (xmax (+ x0 2pi)) (y y0))
      (let ((nz (pv (+ z (* 2pi rotation-number)))))
        (map x y 
             (lambda (nx ny)
               (let ((nx (pv nx)))
                 (cond ((< x0 z zmax)
                        (if (< x0 x xmax)
                            (lp (+ n 1) nz zmin z nx xmin x ny)
                            (if (> x xmax) 1 -1)))
                       ((< zmin z x0)
                        (if (< xmin x x0)
                            (lp (+ n 1) nz z zmax nx x xmax ny)
                            (if (< x xmin) -1 1)))
                       (else 
                        (lp (+ n 1) nz zmin zmax nx xmin xmax ny)))))
             (lambda ()
               (error "Map failed" x y)))))))


;;;; Chapter 5

;;; Makes a canonical point transformation from a 
;;;  time-invariant coordinate transformation T(q)

(define (F->CT F)
  (define (CT H-state)
    (let ((t (time H-state))
          (q (coordinate H-state))
          (p (momentum H-state)))
      (->H-state t
                 (F H-state)
                 (* p (s:inverse 
		       (compatible-shape q)
		       (((partial 1) F) H-state)
		       (compatible-shape p)
		       )))))
  CT)

#|
;;; For display in book

(define ((F->CT F) H-state)
  (->H-state (time H-state)
             (F H-state)
             (* (momentum H-state)
                (invert (((partial 1) F) H-state)))))
|#
#|

(define ((H-central m V) state)
  (let ((x (coordinate state))
	(p (momentum state)))
    (+ (/ (square p) (* 2 m))
       (V (sqrt (square x))))))

(show-expression
 ((compose (H-central 'm (literal-function 'V))
           (F->CT p->r))
  (->H-state 't
             (coordinate-tuple 'r 'phi)
             (momentum-tuple 'p_r 'p_phi))))
(+ (V r)
   (/ (* 1/2 (expt p_r 2)) m)
   (/ (* 1/2 (expt p_phi 2)) (* m (expt r 2))))

|#


(define (canonical? C H Hprime)
  (- (compose (phase-space-derivative H) C)
     (* (D C) (phase-space-derivative Hprime))))

#|
(define (compositional-canonical? C H)
  (- (compose (phase-space-derivative H) C)
     (* (D C) 
        (phase-space-derivative (compose H C)))))
|#

(define (compositional-canonical? C H)
  (canonical? C H (compose H C)))

#|
(simplify
 ((compositional-canonical? (F->CT p->r)
			    (H-central 'm (literal-function 'V)))
  (->H-state 't
	     (coordinate-tuple 'r 'phi)
	     (momentum-tuple 'p_r 'p_phi))))
;Value: (up 0 (up 0 0) (down 0 0))
|#

#|
;;; The shape of (D (F->CT p->r))
(print-expression
 (* ((D (F->CT p->r))
     (->H-state 't
		(coordinate-tuple 'r 'theta)
		(momentum-tuple 'p_r 'p_theta)))
    (->H-state 'dt
	       (coordinate-tuple 'dr 'dtheta)
	       (momentum-tuple 'dp_r 'dp_theta))))
(up
 dt
 (up (+ (* -1 dtheta r (sin theta)) (* dr (cos theta)))
     (+ (* dtheta r (cos theta)) (* dr (sin theta))))
 (down
  (+ (* -1 dtheta p_r (sin theta))
     (* dp_r (cos theta))
     (/ (* -1 dtheta p_theta (cos theta)) r)
     (/ (* -1 dp_theta (sin theta)) r)
     (/ (* dr p_theta (sin theta)) (expt r 2)))
  (+ (* dtheta p_r (cos theta))
     (* dp_r (sin theta))
     (/ (* -1 dtheta p_theta (sin theta)) r)
     (/ (* dp_theta (cos theta)) r)
     (/ (* -1 dr p_theta (cos theta)) (expt r 2)))))
|#

(define (J-func DH)
  (->H-state 0
	     (ref DH 2)
	     (- (ref DH 1))))

(define (T-func DH)
  (->H-state 1 
	     (zero-like (ref DH 2))
	     (zero-like (ref DH 1))))

(define (J*-func s)
  (down 0
	(ref s 2)
	(- (ref s 1))))


#|
(print-expression
 ((D J-func)
  (->H-state 't
	     (coordinate-tuple 'x 'y)
	     (momentum-tuple 'p_x 'p_y))))
(down (up 0 (down 0 0) (up 0 0))
      (down (up 0 (down 0 0) (up -1 0)) (up 0 (down 0 0) (up 0 -1)))
      (up (up 0 (down 1 0) (up 0 0)) (up 0 (down 0 1) (up 0 0))))    
|#

(define (linear-function->multiplier F argument)
  ((D F) argument))

(define ((Phi A) v) (* A v))

(define ((Phi* A) w) (* w A))


(define ((time-independent-canonical? C) s)
  ((- J-func
      (compose (Phi ((D C) s)) 
               J-func
               (Phi* ((D C) s))))
   (compatible-shape s)))

(define (symplectic-two-form zeta1 zeta2)
  (- (* (momentum zeta2) (coordinate zeta1))
     (* (momentum zeta1) (coordinate zeta2))))

(define ((dual-canonical? C) s)
  (let ((DC* (flip-indices ((D C) s))))
    ((- J*-func
	(compose (Phi DC*) 
		 J*-func
		 (Phi* DC*)))
     (typical-object s))))



#|
(print-expression
 ((time-independent-canonical? (F->CT p->r))
  (->H-state 't
	     (coordinate-tuple 'r 'phi)
	     (momentum-tuple 'p_r 'p_phi))))
(up 0 (up 0 0) (down 0 0))


;;; but not all transforms are

(define (a-non-canonical-transform Istate)
  (let ((t (time Istate))
        (theta (coordinate Istate))
	(p (momentum Istate)))
    (let ((x (* p (sin theta)))
	  (p_x (* p (cos theta))))
      (->H-state t x p_x))))

(print-expression
 ((time-independent-canonical? a-non-canonical-transform)
  (->H-state 't 'theta 'p)))
(up 0 (+ (* -1 p x8102) x8102) (+ (* p x8101) (* -1 x8101)))

(print-expression
 ((time-independent-canonical? (polar-canonical 'alpha))
  (->H-state 't 'a 'I)))
(up 0 0 0)

(define (Cmix H-state)
  (let ((t (time H-state))
	(q (coordinate H-state))
	(p (momentum H-state)))
    (->H-state t
	       (coordinate-tuple (ref q 0) (- (ref p 1)))
	       (momentum-tuple   (ref p 0) (ref q 1)))))

(define a-state (->H-state 't 
			   (coordinate-tuple 'x 'y)
			   (momentum-tuple 'p_x 'p_y)))

(print-expression
 ((time-independent-canonical? Cmix)
  a-state))
(up 0 (up 0 0) (down 0 0))

(define (Cmix2 H-state)
  (let ((t (time H-state))
	(q (coordinate H-state))
	(p (momentum H-state)))
    (->H-state t
	       (flip-outer-index p)
	       (- (flip-outer-index q)))))

(print-expression
 ((time-independent-canonical? Cmix2)
  a-state))
(up 0 (up 0 0) (down 0 0))

(define ((F m0 m1) H-state)
  (let ((q (coordinate H-state)))
    (let ((x0 (ref q 0))
	  (x1 (ref q 1)))
      (coordinate-tuple (/ (+ (* m0 x0) (* m1 x1)) (+ m0 m1))
			(- x1 x0)))))

(define ((C m0 m1) state)
  (let ((x (coordinate state))
	(p (momentum state)))
    (let ((x0 (ref x 0))
	  (x1 (ref x 1))
	  (p0 (ref p 0))
	  (p1 (ref p 1)))
      (->H-state 
       (time state)
       (coordinate-tuple (/ (+ (* m0 x0) (* m1 x1)) (+ m0 m1))
			 (- x1 x0))
       (momentum-tuple (+ p0 p1)
		       (/ (- (* m0 p1) (* m1 p0))
			  (+ m0 m1)))))))

(define b-state
  (->H-state
   't
   (coordinate-tuple
    (coordinate-tuple 'x_1 'y_1)
    (coordinate-tuple 'x_2 'y_2))
   (momentum-tuple
    (momentum-tuple 'p_x_1 'p_y_1)
    (momentum-tuple 'p_x_2 'p_y_2))))

(pe (- ((F->CT (F 'm0 'm1)) b-state)
       ((C 'm0 'm1) b-state)))
(up 0 (up (up 0 0) (up 0 0)) (down (down 0 0) (down 0 0)))

(print-expression
 ((time-independent-canonical? (C 'm1 'm2)) b-state))
(up 0 (up (up 0 0) (up 0 0)) (down (down 0 0) (down 0 0)))
|#

#|
(define ((time-independent-canonical? C) s)
  (let ((s* (compatible-shape s)))
    (let ((J (linear-function->multiplier J-func s*)))
      (- J 
	 (* ((D C) s)
	    (* J
	       ((multiplicative-transpose s*) ((D C) s))))))))
|# 

(define ((multiplicative-transpose s) A)
  (linear-function->multiplier (transpose-function A) s))

(define ((transpose-function A) p) (* p A))

#|
(print-expression
 ((time-independent-canonical? (F->CT p->r))
  (->H-state 't
	     (coordinate-tuple 'r 'phi)
	     (momentum-tuple 'p_r 'p_phi))))
(up (up 0 (up 0 0) (down 0 0))
    (up (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0)))
    (down (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0))))


;;; but not all transforms are

(define (a-non-canonical-transform Istate)
  (let ((t (time Istate))
        (theta (coordinate Istate))
	(p (momentum Istate)))
    (let ((x (* p (sin theta)))
	  (p_x (* p (cos theta))))
      (->H-state t x p_x))))

(print-expression
 ((time-independent-canonical? a-non-canonical-transform)
  (->H-state 't 'theta 'p)))
(up (up 0 0 0) (up 0 0 (+ -1 p)) (up 0 (+ 1 (* -1 p)) 0))

(print-expression
 ((time-independent-canonical? (polar-canonical 'alpha))
  (->H-state 't 'a 'I)))
(up (up 0 0 0) (up 0 0 0) (up 0 0 0))

(define (Cmix H-state)
  (let ((t (time H-state))
	(q (coordinate H-state))
	(p (momentum H-state)))
    (->H-state t
	       (coordinate-tuple (ref q 0) (- (ref p 1)))
	       (momentum-tuple   (ref p 0) (ref q 1)))))

(define a-state (->H-state 't 
			   (coordinate-tuple 'x 'y)
			   (momentum-tuple 'p_x 'p_y)))

(print-expression
 ((time-independent-canonical? Cmix)
  a-state))
(up (up 0 (up 0 0) (down 0 0))
    (up (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0)))
    (down (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0))))

(define (Cmix2 H-state)
  (let ((t (time H-state))
	(q (coordinate H-state))
	(p (momentum H-state)))
    (->H-state t
	       (flip-outer-index p)
	       (- (flip-outer-index q)))))

(print-expression
 ((time-independent-canonical? Cmix2)
  a-state))
(up (up 0 (up 0 0) (down 0 0))
    (up (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0)))
    (down (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0))))

(define ((C m0 m1) state)
  (let ((x (coordinate state))
	(p (momentum state)))
    (let ((x0 (ref x 0))
	  (x1 (ref x 1))
	  (p0 (ref p 0))
	  (p1 (ref p 1)))
      (->H-state 
       (time state)
       (coordinate-tuple (/ (+ (* m0 x0) (* m1 x1)) (+ m0 m1))
			 (- x1 x0))
       (momentum-tuple (+ p0 p1)
		       (/ (- (* m0 p1) (* m1 p0))
			  (+ m0 m1)))))))

(define b-state
  (->H-state
   't
   (coordinate-tuple
    (coordinate-tuple 'x_1 'y_1)
    (coordinate-tuple 'x_2 'y_2))
   (momentum-tuple
    (momentum-tuple 'p_x_1 'p_y_1)
    (momentum-tuple 'p_x_2 'p_y_2))))

(print-expression
 ((time-independent-canonical? (C 'm1 'm2)) b-state))
(up
 (up 0 (up (up 0 0) (up 0 0)) (down (down 0 0) (down 0 0)))
 (up
  (up (up 0 (up (up 0 0) (up 0 0)) (down (down 0 0) (down 0 0)))
      (up 0 (up (up 0 0) (up 0 0)) (down (down 0 0) (down 0 0))))
  (up (up 0 (up (up 0 0) (up 0 0)) (down (down 0 0) (down 0 0)))
      (up 0 (up (up 0 0) (up 0 0)) (down (down 0 0) (down 0 0)))))
 (down
  (down (up 0 (up (up 0 0) (up 0 0)) (down (down 0 0) (down 0 0)))
        (up 0 (up (up 0 0) (up 0 0)) (down (down 0 0) (down 0 0))))
  (down (up 0 (up (up 0 0) (up 0 0)) (down (down 0 0) (down 0 0)))
        (up 0 (up (up 0 0) (up 0 0)) (down (down 0 0) (down 0 0))))))

|#

;;; in terms of matrices

(define ((symplectic? C) s)
  (let ((s* (compatible-shape s)))
    (let ((J (s->m s* ((D J-func) s*) s*))
	  (DC (s->m s* ((D C) s) s)))
      (- J (* DC J (m:transpose DC))))))

#|
;;; By the way

(define b-state
  (->H-state
   't
   (coordinate-tuple
    (coordinate-tuple 'x_1 'y_1)
    (coordinate-tuple 'x_2 'y_2))
   (momentum-tuple
    (momentum-tuple 'p_x_1 'p_y_1)
    (momentum-tuple 'p_x_2 'p_y_2))))

(let* ((s b-state)
       (s* (compatible-shape s))
       (A (typical-object (outer-product s s*))))
  (pe (- (s->m s ((D (Phi* A)) s*) s*)
	 (m:transpose (s->m s* ((D (Phi A)) s) s)))))
(matrix-by-rows (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0))

|#

#|
(print-expression
 ((symplectic? (F->CT p->r))
  (->H-state 't
	     (coordinate-tuple 'r 'phi)
	     (momentum-tuple 'p_r 'p_phi))))
(matrix-by-rows (list 0 0 0 0 0)
                (list 0 0 0 0 0)
                (list 0 0 0 0 0)
                (list 0 0 0 0 0)
                (list 0 0 0 0 0))


;;; but not all transforms are

(define (a-non-canonical-transform Istate)
  (let ((t (time Istate))
        (theta (coordinate Istate))
	(p (momentum Istate)))
    (let ((x (* p (sin theta)))
	  (p_x (* p (cos theta))))
      (->H-state t x p_x))))

(print-expression
 ((symplectic? a-non-canonical-transform)
  (->H-state 't 'theta 'p)))
(matrix-by-rows (list 0 0 0)
		(list 0 0 (+ 1 (* -1 p)))
		(list 0 (+ -1 p) 0))

(print-expression
 ((symplectic? (polar-canonical 'alpha))
  (->H-state 't 'a 'I)))
(matrix-by-rows (list 0 0 0)
		(list 0 0 0)
		(list 0 0 0))

(define (Cmix H-state)
  (let ((t (time H-state))
	(q (coordinate H-state))
	(p (momentum H-state)))
    (->H-state t
	       (coordinate-tuple (ref q 0) (- (ref p 1)))
	       (momentum-tuple   (ref p 0) (ref q 1)))))

(define a-state (->H-state 't 
			   (coordinate-tuple 'x 'y)
			   (momentum-tuple 'p_x 'p_y)))

(print-expression
 ((symplectic? Cmix)
  a-state))
(matrix-by-rows (list 0 0 0 0 0)
                (list 0 0 0 0 0)
                (list 0 0 0 0 0)
                (list 0 0 0 0 0)
                (list 0 0 0 0 0))

(define (Cmix2 H-state)
  (let ((t (time H-state))
	(q (coordinate H-state))
	(p (momentum H-state)))
    (->H-state t
	       (flip-outer-index p)
	       (- (flip-outer-index q)))))

(print-expression
 ((symplectic? Cmix2)
  a-state))
(matrix-by-rows (list 0 0 0 0 0)
                (list 0 0 0 0 0)
                (list 0 0 0 0 0)
                (list 0 0 0 0 0)
                (list 0 0 0 0 0))


(define ((C m0 m1) state)
  (let ((x (coordinate state))
	(p (momentum state)))
    (let ((x0 (ref x 0))
	  (x1 (ref x 1))
	  (p0 (ref p 0))
	  (p1 (ref p 1)))
      (->H-state 
       (time state)
       (coordinate-tuple (/ (+ (* m0 x0) (* m1 x1)) (+ m0 m1))
			 (- x1 x0))
       (momentum-tuple (+ p0 p1)
		       (/ (- (* m0 p1) (* m1 p0))
			  (+ m0 m1)))))))

(define b-state
  (->H-state
   't
   (coordinate-tuple
    (coordinate-tuple 'x_1 'y_1)
    (coordinate-tuple 'x_2 'y_2))
   (momentum-tuple
    (momentum-tuple 'p_x_1 'p_y_1)
    (momentum-tuple 'p_x_2 'p_y_2))))

(print-expression
 ((symplectic? (C 'm1 'm2)) b-state))
(matrix-by-rows (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0)
                (list 0 0 0 0 0 0 0 0 0))
|#

#|

(print-expression
 ((dual-canonical? (F->CT p->r))
  (->H-state 't
	     (coordinate-tuple 'r 'phi)
	     (momentum-tuple 'p_r 'p_phi))))
(down 0 (down 0 0) (up 0 0))

;;; but not all transforms are

(define (a-non-canonical-transform Istate)
  (let ((t (time Istate))
        (theta (coordinate Istate))
	(p (momentum Istate)))
    (let ((x (* p (sin theta)))
	  (p_x (* p (cos theta))))
      (->H-state t x p_x))))

(print-expression
 ((dual-canonical? a-non-canonical-transform)
  (->H-state 't 'theta 'p)))
(down 0 (+ (* -1 p x11190) x11190) (+ (* p x11189) (* -1 x11189)))


(print-expression
 ((dual-canonical? (polar-canonical 'alpha))
  (->H-state 't 'a 'I)))
(down 0 0 0)

(define (Cmix H-state)
  (let ((t (time H-state))
	(q (coordinate H-state))
	(p (momentum H-state)))
    (->H-state t
	       (coordinate-tuple (ref q 0) (- (ref p 1)))
	       (momentum-tuple   (ref p 0) (ref q 1)))))

(define a-state (->H-state 't 
			   (coordinate-tuple 'x 'y)
			   (momentum-tuple 'p_x 'p_y)))
(print-expression
 ((dual-canonical? Cmix)
  a-state))
(down 0 (down 0 0) (up 0 0))

(define (Cmix2 H-state)
  (let ((t (time H-state))
	(q (coordinate H-state))
	(p (momentum H-state)))
    (->H-state t
	       (flip-outer-index p)
	       (- (flip-outer-index q)))))

(print-expression
 ((dual-canonical? Cmix2)
  a-state))
(down 0 (down 0 0) (up 0 0))


(define ((C m0 m1) state)
  (let ((x (coordinate state))
	(p (momentum state)))
    (let ((x0 (ref x 0))
	  (x1 (ref x 1))
	  (p0 (ref p 0))
	  (p1 (ref p 1)))
      (->H-state 
       (time state)
       (coordinate-tuple (/ (+ (* m0 x0) (* m1 x1)) (+ m0 m1))
			 (- x1 x0))
       (momentum-tuple (+ p0 p1)
		       (/ (- (* m0 p1) (* m1 p0))
			  (+ m0 m1)))))))

(define b-state
  (->H-state
   't
   (coordinate-tuple
    (coordinate-tuple 'x_1 'y_1)
    (coordinate-tuple 'x_2 'y_2))
   (momentum-tuple
    (momentum-tuple 'p_x_1 'p_y_1)
    (momentum-tuple 'p_x_2 'p_y_2))))

(print-expression
 ((dual-canonical? (C 'm1 'm2)) b-state))
(down 0 (down (down 0 0) (down 0 0)) (up (up 0 0) (up 0 0)))

|#

#|
(define (T v)
  (* (down (up 'a 'c) (up 'b 'd)) v))

(pe (T (up 'x 'y)))
(up (+ (* a x) (* b y)) (+ (* c x) (* d y)))

(pe (* (* (down 'p_x 'p_y) ((D T) (up 'x 'y))) (up 'v_x 'v_y)))
(+ (* a p_x v_x) (* b p_x v_y) (* c p_y v_x) (* d p_y v_y))


(pe (* (down 'p_x 'p_y) (* ((D T) (up 'x 'y)) (up 'v_x 'v_y))))
(+ (* a p_x v_x) (* b p_x v_y) (* c p_y v_x) (* d p_y v_y))

(pe (* (* ((multiplicative-transpose (down 'p_x 'p_y)) ((D T) (up 'x 'y)))
	  (down 'p_x 'p_y))
       (up 'v_x 'v_y)))
(+ (* a p_x v_x) (* b p_x v_y) (* c p_y v_x) (* d p_y v_y))

;;; But strangely enough...
(pe (* (* (down 'p_x 'p_y)
	  ((multiplicative-transpose (down 'p_x 'p_y)) ((D T) (up 'x 'y))))
       (up 'v_x 'v_y)))
(+ (* a p_x v_x) (* b p_x v_y) (* c p_y v_x) (* d p_y v_y))
|#


#|
(define ((canonical-K? C K) s)
  (let ((DCs ((D C) s))
	(s* (compatible-shape s)))
    (- (T-func s*)
       (* DCs ((phase-space-derivative K) s)))))
|#

(define ((canonical-K? C K) s)
  (let ((s* (compatible-shape s)))
    (- (T-func s*)
       (+ (* ((D C) s) (J-func ((D K) s)))
	  (((partial 0) C) s)))))

(define ((F->K F) s)
  (* -1 (((partial 0) F) s) (momentum ((F->CT F) s))))

#|

(define ((rotating n) state)
  (let ((t (time state))
	(q (coordinate state)))
    (let ((x (ref q 0))
	  (y (ref q 1))
	  (z (ref q 2)))
      (coordinate-tuple (+ (* (cos (* n t)) x) (* (sin (* n t)) y))
			(- (* (cos (* n t)) y) (* (sin (* n t)) x))
			z))))

(define (C-rotating n) (F->CT (rotating n)))

(define ((K n) s)
  (let ((q (coordinate s))
	(p (momentum s)))
    (let ((x (ref q 0)) (y (ref q 1))
	  (px (ref p 0)) (py (ref p 1)))
      (* n (- (* x py) (* y px))))))

(define a-state 
  (->H-state 't 
	     (coordinate-tuple 'x 'y 'z)
	     (momentum-tuple 'p_x 'p_y 'p_z)))


(pe ((canonical-K? (C-rotating 'n) (K 'n)) a-state))
(up 0 (up 0 0 0) (down 0 0 0))

;;; or getting K directly from F
(pe ((canonical-K? (C-rotating 'n) (F->K (rotating 'n))) a-state))
(up 0 (up 0 0 0) (down 0 0 0))

(pe ((- (F->K (rotating 'n))
	(K 'n))
     a-state))
0

;;; Poisson brackets are compositional
(pe ((- (compose (Poisson-bracket Lx Ly) (C-rotating 'n))
	(Poisson-bracket (compose Lx (C-rotating 'n))
			 (compose Ly (C-rotating 'n))) )
     3-state))
0

;;; not all K's work

(define ((bad-K n) s)
  (- ((K n) s)))

(pe ((canonical-K? (C-rotating 'n) (bad-K 'n)) a-state))
(up
 0
 (up (+ (* 2 n x (sin (* n t))) (* -2 n y (cos (* n t))))
     (+ (* 2 n x (cos (* n t))) (* 2 n y (sin (* n t))))
     0)
 (down (+ (* 2 n p_x (sin (* n t))) (* -2 n p_y (cos (* n t))))
       (+ (* 2 n p_x (cos (* n t))) (* 2 n p_y (sin (* n t))))
       0))
|#

;;;----------------------------------------------------------------
;;; Poisson brackets in terms of J
#|

guaranteed to work only for scalar valued functions

(define ((PB f g) s)
  (* ((D f) s) (J-func ((D g) s))))

(define a-state 
  (->H-state 't 
	     (coordinate-tuple 'x 'y 'z)
	     (momentum-tuple 'p_x 'p_y 'p_z)))

(pe ((- (Poisson-bracket Lx Ly) Lz) a-state))
0
(pe ((- (PB Lx Ly) Lz) a-state))
0

(define ((PB f g) s)
  (let ((J ((D J-func) ((D g) s))))
    (* ((D f) s) (* J ((D g) s)))))

(define ((PB f g) s)
  (let ((J (linear-function->multiplier J-func ((D g) s))))
    (* ((D f) s) (* J ((D g) s)))))

(pe 
 (- ((Poisson-bracket (H-harmonic 'm 'k)
		      ((component 0) coordinate)) 
     a-state)
    ((PB (H-harmonic 'm 'k)
	 (compose (component 0) coordinate))
     a-state)
    ))
0

(pe 
 (- ((Poisson-bracket (H-harmonic 'm 'k) coordinate) 
     a-state)
    ((PB (H-harmonic 'm 'k) coordinate)
     a-state)
    ))
(up 0 0 0)

(pe ((PB momentum (H-harmonic 'm 'k))
     a-state))
(down (* -1 k x) (* -1 k y) (* -1 k z))

(pe ((PB coordinate (H-harmonic 'm 'k))
     a-state))
(up (/ p_x m) (/ p_y m) (/ p_z m))

|#

#|
;;;----------------------------------------------------------------
;;; generating functions

;;; identity
(define (F2-identity t q p-prime)
  (+ (* (ref q 0) (ref p-prime 0))
     (* (ref q 1) (ref p-prime 1))))

;;; flip q & p
(define (F1-flip t q q-prime)
  (+ (* (ref q 0) (ref q-prime 0))
     (* (ref q 1) (ref q-prime 1))))

(define ((F2->C-inv F2) H-state)
  (let ((t (time H-state))
	(q (coordinate H-state))
	(p (momentum H-state)))
    (let ((p-func
	   (lambda (p-prime)
	     (((partial 1) F2) t q p-prime))))
      (let ((p-prime ((linear-inverse p-func) p)))
	(let ((q-prime (((partial 2) F2) t q p-prime)))
	  (->H-state t q-prime p-prime))))))
    
;;; x = f(x') is linear
(define ((linear-inverse f) p)
  (let ((b (f (zero-like p)))
	(a ((D f) (zero-like p))))
    (* (- p b)
       (s:inverse (compatible-shape p) a p))))

(define a-state
  (->H-state 't 
	     (coordinate-tuple 'x 'y)
	     (momentum-tuple 'p_x 'p_y)))

(print-expression ((F2->C-inv F2-identity) a-state))
(up t (up x y) (down p_x p_y))

(define b-state
  (->H-state
   't
   (coordinate-tuple
    (coordinate-tuple 'x_1 'y_1)
    (coordinate-tuple 'x_2 'y_2))
   (momentum-tuple
    (momentum-tuple 'p_x_1 'p_y_1)
    (momentum-tuple 'p_x_2 'p_y_2))))

(pe ((F2->C-inv F2-identity) b-state))
(up t
    (up (up x_1 y_1) (up x_2 y_2))
    (down (down p_x_1 p_y_1) (down p_x_2 p_y_2)))
 
(print-expression
 ((time-independent-canonical? (F2->C-inv F2-identity))
  a-state))
(up 0 (up 0 0) (down 0 0))

(pe (((partial 1) F2-identity) 
     't 
     (coordinate b-state)
     (typical-object (momentum b-state))))
(down (down x14885 x14886) (down x14887 x14888))

|#
;;; end of J.scm 
;;; ----------------------------------------------------------------

(define (tz->tqp t z)
  (let* ((2n (vector-length z))
	 (n (quotient 2n 2)))
    (assert (even? 2n))
    (if (fix:= n 1)
	(->H-state t (vector-ref z 0) (vector-ref z 1))
	(->H-state t
		   (vector->column (subvector z 0 n))
		   (vector->row (subvector z n 2n))))))

(define (z->tqp z)
  (tz->tqp 'unknown-time z))

(define (tqp->z tqp)
  (let ((q (coordinate tqp))
	(p (momentum tqp)))
    (if (and (column? q) (row? p))
	(vector->up
	 (vector-append (column->vector q)
			(row->vector p)))
	(up q p))))

(define (tqp->tz tqp)
  (let ((t (time tqp))
	(q (coordinate tqp))
	(p (momentum tqp)))
    (if (and (column? q) (row? p))
	(up t
	    (vector->up
	     (vector-append (column->vector q)
			    (row->vector p))))
	(up t (up q p)))))


(define ((symplectic-transform? C) arg)
  (symplectic-matrix?
   (qp-submatrix
    (s->m (compatible-shape arg)
	  ((D C) arg)
	  arg))))

(define (qp-submatrix m)
  (m:submatrix m 1 (m:num-rows m) 1 (m:num-cols m)))


(define (symplectic-matrix? M)
  (let ((2n (m:dimension M)))
    (if (not (even? 2n))
	(error "Wrong type -- SYMPLECTIC-MATRIX?" M))
    (let ((J (symplectic-unit (quotient 2n 2))))
      (- J (* M J (m:transpose M))))))

(define (symplectic-unit n)
  (let ((2n (fix:* 2 n)))
    (m:generate 2n 2n
       (lambda (a b) 
	 (cond ((fix:= (fix:+ a n) b) 1)
	       ((fix:= (fix:+ b n) a) -1)
	       (else 0))))))

#|
;;; For example, point transforms are canonical

(print-expression
 ((symplectic-transform? (F->CT p->r))
  (->H-state 't
	     (coordinate-tuple 'r 'theta)
	     (momentum-tuple 'p_r 'p_theta))))
(matrix-by-rows (list 0 0 0 0) (list 0 0 0 0) (list 0 0 0 0) (list 0 0 0 0))
|#

#|
(print-expression
 ((symplectic-transform? a-non-canonical-transform)
  (->H-state 't 'theta 'p)))
(matrix-by-rows (list 0 (+ 1 (* -1 p))) (list (+ -1 p) 0))
|#


;;; One particularly useful canonical transform is the 
;;;  Poincare transform, which is good for simplifying 
;;;  oscillators.

(define ((polar-canonical alpha) Istate)
  (let ((t (time Istate))
        (theta (coordinate Istate))
        (I (momentum Istate)))
    (let ((x (* (sqrt (/ (* 2 I) alpha)) (sin theta)))
	  (p_x (* (sqrt (* 2 alpha I)) (cos theta))))
      (->H-state t x p_x))))

(define ((polar-canonical-inverse alpha) s)
  (let ((t (time s))
	(x (coordinate s))
	(p (momentum s)))
    (let ((I (/ (+ (* alpha (square x))
		   (/ (square p) alpha))
		2)))
      (let ((theta (atan (/ x (sqrt (/ (* 2 I) alpha)))
			 (/ p (sqrt (* 2 I alpha))))))
	(->H-state t theta I)))))

#|

(pe
 ((compose (polar-canonical-inverse 'alpha)
	   (polar-canonical 'alpha))
  (->H-state 't 'x 'p)))
(up t x p)

|#

#|
;;; It is clearly canonical.

(print-expression
 ((symplectic-transform? (polar-canonical 'alpha))
  (->H-state 't 'a 'I)))
(matrix-by-rows (list 0 0) (list 0 0))
|#

;;; time evolution transformations

(define ((shift-t delta-t) state)
  (->H-state 
   (+ (time state) delta-t)
   (coordinate state)
   (momentum state)))

(define ((C->Cp C) delta-t)
  (compose (C delta-t) (shift-t (- delta-t))))

(define ((H->Hp delta-t) H)
  (compose H (shift-t (- delta-t))))


;;; 

(define (commutator op1 op2)
  (- (* op1 op2) (* op2 op1)))

(define (anticommutator op1 op2)
  (+ (* op1 op2) (* op2 op1)))

;;; We define the Lie derivative of F, as a derivative-like operator,
;;;  relative to the given Hamiltonian-like function, H.

(define (Lie-derivative H)
  (make-operator
   (lambda (F)
     (Poisson-bracket F H))
   `(Lie-derivative ,H)))

;;; the flow derivative generalizes the Lie derivative to 
;;; allow for time dependent H and F ---
;;; computes the "time" derivative of F along the flow specified by H

(define (flow-derivative H)
  (make-operator 
   (lambda (F)
     (+ ((partial 0) F)
	(Poisson-bracket F H)))
   `(flow-derivative ,H)))


;;; The Lie transform is just the time-advance operator using the Lie
;;;  derivative.
#|
(define (((Lie-transform H delta-t) F) H-state)
  ((((make-operator 
      (exp (* delta-t (Lie-derivative H)))
      `(Lie-transform ,H ,delta-t))
     F)
    H-state)
   1))
|#

(define (Lie-transform H delta-t)
  (make-operator 
    (exp (* delta-t (Lie-derivative H)))
    `(Lie-transform ,H ,delta-t)))
  

;;; The generalization of Lie-transform to include time dependence.
#|
(define (((flow-transform H delta-t) F) H-state)
  ((((make-operator 
      (exp (* delta-t (flow-derivative H)))
      `(flow-transform ,H ,delta-t))
     F)
    H-state)
   1))
|#

(define (flow-transform H delta-t)
  (make-operator 
   (exp (* delta-t (flow-derivative H)))
   `(flow-transform ,H ,delta-t)))

#|
;;; The general solution for a trajectory is:
;;;
;;;  q(t,q0,p0) = A(q0,p0) cos (sqrt(k/m)*t + phi(q0,p0))
;;;
;;;  where A(q0,p0) = sqrt(2/k)*sqrt(p0^2/(2*m) + (k/2)*q0^2)
;;;                 = sqrt((2/k)*E0)
;;;
;;;  and   phi(q0,p0) = - atan((1/sqrt(k*m))*(p0/q0))
;;;
;;; Thus, with initial conditions q0, p0
;;;   we should get q(t) = q0*cos(sqrt(k/m)*t)+p0*sin(sqrt(k/m)*t)
;;;
;;; We can expand this as a Lie series:

(define ((H-harmonic m k) state)
  (let ((q (coordinate state))
	(p (momentum state)))
    (+ (/ (square p) (* 2 m))
       (* 1/2 k (square q)))))

;;; This works, but it takes forever! -- hung in deriv, not in simplify!

(series:for-each print-expression
 (((Lie-transform (H-harmonic 'm 'k) 'dt)
   state->q)
  (->H-state 0 'x_0 'p_0))
 6)
x_0
(/ (* dt p_0) m)
(/ (* -1/2 (expt dt 2) k x_0) m)
(/ (* -1/6 (expt dt 3) k p_0) (expt m 2))
(/ (* 1/24 (expt dt 4) (expt k 2) x_0) (expt m 2))
(/ (* 1/120 (expt dt 5) (expt k 2) p_0) (expt m 3))
;Value: ...

(series:for-each print-expression
 (((Lie-transform (H-harmonic 'm 'k) 'dt)
   momentum)
  (->H-state 0 'x_0 'p_0))
 6)
p_0
(* -1 dt k x_0)
(/ (* -1/2 (expt dt 2) k p_0) m)
(/ (* 1/6 (expt dt 3) (expt k 2) x_0) m)
(/ (* 1/24 (expt dt 4) (expt k 2) p_0) (expt m 2))
(/ (* -1/120 (expt dt 5) (expt k 3) x_0) (expt m 2))
;Value: ...

(series:for-each print-expression
 (((Lie-transform (H-harmonic 'm 'k) 'dt)
   (H-harmonic 'm 'k))
  (->H-state 0 'x_0 'p_0))
 6)
(/ (+ (* 1/2 k m (expt x_0 2)) (* 1/2 (expt p_0 2))) m)
0
0
0
0
0
;Value: ...

(series:for-each print-expression
 (((Lie-transform
    (H-central-polar 'm (literal-function 'U))
    'dt)
   state->q)
  (->H-state 0
	      (coordinate-tuple 'r_0 'phi_0)
	      (momentum-tuple 'p_r_0 'p_phi_0)))
 4)
(up r_0 phi_0)
(up (/ (* dt p_r_0) m) (/ (* dt p_phi_0) (* m (expt r_0 2))))
(up
 (+ (/ (* -1/2 (expt dt 2) ((D U) r_0)) m)
    (/ (* 1/2 (expt dt 2) (expt p_phi_0 2)) (* (expt m 2) (expt r_0 3))))
 (/ (* -1 (expt dt 2) p_phi_0 p_r_0) (* (expt m 2) (expt r_0 3))))
(up
 (+
  (/ (* -1/6 (expt dt 3) p_r_0 (((expt D 2) U) r_0)) (expt m 2))
  (/ (* -1/2 (expt dt 3) (expt p_phi_0 2) p_r_0) (* (expt m 3) (expt r_0 4))))
 (+ (/ (* 1/3 (expt dt 3) p_phi_0 ((D U) r_0)) (* (expt m 2) (expt r_0 3)))
    (/ (* -1/3 (expt dt 3) (expt p_phi_0 3)) (* (expt m 3) (expt r_0 6)))
    (/ (* (expt dt 3) p_phi_0 (expt p_r_0 2)) (* (expt m 3) (expt r_0 4)))))
;Value: ...


;;; I left this one that uses the Lagrangian because it appears to be 
;;; used for timings
(show-time
 (lambda ()
   (series:print
    (((Lie-transform
       (Lagrangian->Hamiltonian
	(L-central-polar 'm (lambda (r) (- (/ 'GM r)))))
       'dt)
      state->q)
     (->H-state 0
		 (coordinate-tuple 'r_0 'phi_0)
		 (momentum-tuple 'p_r_0 'p_phi_0)))
    4)))
(up r_0 phi_0)
(up (/ (* dt p_r_0) m) (/ (* dt p_phi_0) (* m (expt r_0 2))))
(up
 (+ (/ (* -1/2 GM (expt dt 2)) (* m (expt r_0 2)))
    (/ (* 1/2 (expt dt 2) (expt p_phi_0 2)) (* (expt m 2) (expt r_0 3))))
 (/ (* -1 (expt dt 2) p_phi_0 p_r_0) (* (expt m 2) (expt r_0 3))))
(up
 (+
  (/ (* 1/3 GM (expt dt 3) p_r_0) (* (expt m 2) (expt r_0 3)))
  (/ (* -1/2 (expt dt 3) (expt p_phi_0 2) p_r_0) (* (expt m 3) (expt r_0 4))))
 (+ (/ (* (expt dt 3) p_phi_0 (expt p_r_0 2)) (* (expt m 3) (expt r_0 4)))
    (/ (* 1/3 GM (expt dt 3) p_phi_0) (* (expt m 2) (expt r_0 5)))
    (/ (* -1/3 (expt dt 3) (expt p_phi_0 3)) (* (expt m 3) (expt r_0 6)))))

;;; HOD     
;;;         process time: 14590 (13610 RUN + 980 GC); real time: 14588
;;; PLANET003 600MHz PIII
;;;         process time: 19610 (17560 RUN + 2050 GC); real time: 19610
;;; HEIFETZ xeon 400MHz 512K
;;;         process time: 27380 (24250 RUN + 3130 GC); real time: 27385
;;; GEVURAH 300 MHz
;;;         process time: 36070 (33800 RUN + 2270 GC); real time: 36072
;;; MAHARAL 
;;;         process time: 56390 (50970 RUN + 5420 GC); real time: 56386
;;; ACTION1 200MHz Pentium Pro
;;;         process time: 55260 (49570 RUN + 5690 GC); real time: 55257
;;; PPA     200MHz Pentium Pro
;;;         process time: 58840 (56500 RUN + 2340 GC); real time: 59165
;;; ZOHAR   33MHz 486
;;;         process time: 463610 (443630 RUN + 19980 GC); real time: 485593
|#

;;;*************************
#|
(define r (literal-function 'r))
(define theta (literal-function 'theta))
(define phi (literal-function 'phi))

;;; Better to formulate this in terms of tensors.

(define ((T-curvilinear m g) state)     ; g is metric matrix
  (let ((qdot (velocity state)))
    (* 1/2 m (* qdot g qdot))))

(show-expression
 (((T-curvilinear 'm
		  (row (row 1 0)
		       (row 0 (square r))))
   (->local
    (lambda (t) t)
    (coordinate-tuple r theta)
    (velocity-tuple (D r) (D theta))))
  't))
(+ (* 1/2 (expt ((D r) t) 2) m)
   (* 1/2 (expt (r t) 2) (expt ((D theta) t) 2) m))
|#

