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


(define ((L-free-particle mass) local)
  (let ((v (velocity local)))
    (* 1/2 mass (square v))))

(define ((L-harmonic m k) local)
  (let ((q (coordinate local)) 
        (v (velocity local)))
    (- (* 1/2 m (square v))
       (* 1/2 k (square q)))))

(define ((L-uniform-acceleration m g) local)
  (let ((q (coordinate local))
        (v (velocity local)))
    (let ((y (ref q 1)))
      (- (* 1/2 m (square v)) (* m g y)))))

(define ((L-central-rectangular m V) local)
  (let ((q (coordinate local))
        (v (velocity local)))
    (- (* 1/2 m (square v))
       (V (sqrt (square q))))))

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
#|
(define (L-central-polar m V)
  (compose (L-central-rectangular m V)
	   (F->C p->r)))
|#

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

(define ((Lf m g) local)
  (let ((q (coordinate local))
        (v (velocity local)))
    (let ((h (ref q 1)))
      (- (* 1/2 m (square v)) (* m g h)))))

(define ((L-coupled-harmonic m k) state)
  (let ((q (coordinate state))
	(qdot (velocity state)))
    (- (* 1/2 qdot m qdot)
       (* 1/2 q k q))))


(define ((periodic-drive amplitude frequency phase) t)
  (* amplitude (cos (+ (* frequency t) phase))))

(define (L-periodically-driven-pendulum m l g a omega)
  (let ((ys (periodic-drive a omega 0)))
    (L-pend m l g ys)))



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
	    (* b thetadot (vector (cos theta) (sin theta))))
	   (pend-vel (+ rel-pend-vel (vector xdot 0)))
	   (Tpend (* 1/2 m2 (square pend-vel)))
	   (Tsupport (* 1/2 m1 (square xdot)))
	   (V (- (* m2 g b (cos theta)))))
      (+ Tpend Tsupport (- V)))))

(define ((L-pendulum g m l) state)
  (let ((theta (coordinate state))
	(thetadot (velocity state)))
    (+ (* 1/2 m (square (* l thetadot)))
       (* g m l (cos theta)))))

(define ((Rayleigh-dissipation k) state)
  (let ((qdot (velocity state)))
    (* qdot k qdot)))

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

(define ((H-central m V) state)
  (let ((x (coordinate state))
	(p (momentum state)))
    (+ (/ (square p) (* 2 m))
       (V (sqrt (square x))))))

(define ((H-central-polar m V) state)
  (let ((q (coordinate state))
        (p (momentum state)))
    (let ((r (ref q 0))
          (phi (ref q 1))
          (pr (ref p 0))
          (pphi (ref p 1)))
      (+ (/ (+ (square pr)
	       (square (/ pphi r)))
	    (* 2 m))
         (V r)))))


(define ((L-two-particle m1 m2) local)
  (let ((x (coordinate local))
	(v (velocity local))
	(V (literal-function 'V (-> (X (^ Real 2) (^ Real 2)) Real))))
    (let ((x1 (ref x 0)) (x2 (ref x 1))
          (v1 (ref v 0)) (v2 (ref v 1)))
      (- (+ (* 1/2 m1 (square v1))
	    (* 1/2 m2 (square v2)))
	 (V x1 x2)))))


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


(define ((H-harmonic m k) state)
  (let ((q (coordinate state))
	(p (momentum state)))
    (+ (/ (square p) (* 2 m))
       (* 1/2 k (square q)))))


(define ((T-curvilinear m g) state)     ; g is metric matrix
  (let ((qdot (velocity state)))
    (* 1/2 m (* qdot g qdot))))


(define ((H-rectangular m V) H-state)
  (let ((q (coordinate H-state))
        (p (momentum H-state)))
    (+ (/ (square p) (* 2 m))
       (V (ref q 0) (ref q 1)))))

(define ((L-rectangular m V) local)
  (let ((q (coordinate local))
        (qdot (velocity local)))
    (- (* 1/2 m (square qdot))
       (V (ref q 0) (ref q 1)))))