#| -*-Scheme-*-

Copyright (C) 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994,
    1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005,
    2006, 2007, 2008, 2009, 2010, 2011 Massachusetts Institute of
    Technology

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

;;;; Special Relativity

;;; To make :c evaluate to the symbol :c execute
;;; (symbolic-constants #f (list (get-constant-data ':c)))

;;; Lorentz transformation for beta = (vx/c, vy/c, vz/c)

(define ((4vector-boost beta) 4vector-prime)
  (let ((delta-ct-prime (4vector->ct 4vector-prime))
	(delta-x-prime (4vector->space 4vector-prime)))
    (let ((betasq (square beta)))
      (let ((bx (dot-product beta delta-x-prime))
	    (gamma (/ 1 (sqrt (- 1 betasq)))))
	(let ((alpha (/ (- gamma 1) betasq)))
	  (let ((delta-ct
		 (* gamma (+ delta-ct-prime bx)))
		(delta-x
		 (+ (* gamma beta delta-ct-prime)
		    delta-x-prime
		    (* alpha beta bx))))
	    (make-4vector delta-ct delta-x)))))))

(define ((4vector-boost-safe beta) 4vector-prime)
  (if (zero? beta)
      4vector-prime
      ((4vector-boost beta) 4vector-prime)))

(define (make-4vector ct space)
  (up ct (ref space 0) (ref space 1) (ref space 2)))

(define (4vector->ct v)
  (ref v 0))

(define (4vector->space v)
  (up (ref v 1) (ref v 2) (ref v 3)))

(define (flat-spacetime-interval 4vector)
  (sqrt (- (square (4vector->ct 4vector))
	   (square (4vector->space 4vector)))))

(define (4vector-interval 4vector)
  (sqrt (- (square (4vector->space 4vector))
	   (square (4vector->ct 4vector)))))

(define (4vector-metric 4v1 4v2)
  (- (dot-product (4vector->space 4v1)
		  (4vector->space 4v2))
     (* (4vector->ct 4v1)
	(4vector->ct 4v2))))

#|
;;; Test of booster, with :c symbolic with units:
(symbolic-constants #t (list (get-constant-data ':c)))

(define e-prime
  (make-4vector
   (* :c (& 'tp second))
   (up (& 'xp meter) (& 'yp meter) (& 'zp meter))))

(4vector-interval e-prime)
#|
(& (sqrt (+ (* -1 (expt :c 2) (expt tp 2))
	    (expt xp 2)
	    (expt yp 2)
	    (expt zp 2)))
   meter)
|#

;;; Revert to unitless but symbolic :c
(symbolic-constants #f (list (get-constant-data ':c)))

(define e-prime
  (make-4vector (* :c 'tp) (up 'xp 'yp 'zp)))

(4vector-interval e-prime)
#|
(sqrt (+ (* -1 (expt :c 2) (expt tp 2))
	 (expt xp 2)
	 (expt yp 2)
	 (expt zp 2)))
|#

	
;;; Boost by an arbitrary spatial velocity:
			     
(4vector-interval
 ((4vector-boost (up 'vx 'vy 'vz)) e-prime))
#|
(sqrt (+ (* -1 (expt :c 2) (expt tp 2))
	 (expt xp 2)
	 (expt yp 2)
	 (expt zp 2)))
|#

;;; So the boost of an arbitrary 4vector by an arbitrary
;;; velocity preserves the flat-spacetime length of the 4vector.
|#

(define (((boost SR-orthonormal-basis) 3velocity) 4vector-field)
  (define ((v SR-function) event)
    (* ((compose (4vector-boost 3velocity)
		 ((basis->1form-basis SR-orthonormal-basis) 4vector-field))
	event)
       (((basis->vector-basis SR-orthonormal-basis) SR-function)
	event)))
  (procedure->vector-field v))

(define (flat-spacetime-metric)
  (let ((b (coordinate-system->1form-basis spacetime-rect)))
    (let ((dct (ref b 0))
	  (dx (ref b 1))
	  (dy (ref b 2))
	  (dz (ref b 3)))
      (define (g u v)
	(+ (* -1 (dct u) (dct v))
	   (* (dx u) (dx v))
	   (* (dy u) (dy v))
	   (* (dz u) (dz v))))
      (declare-argument-types! g
	(list vector-field? vector-field?))
      g)))

(define (spacetime-orthonormal-basis)
  (let ((vs (coordinate-system->vector-basis spacetime-rect))
	(ws (coordinate-system->1form-basis spacetime-rect)))
    (let ((d/dct (ref vs 0))
	  (d/dx (ref vs 1))
	  (d/dy (ref vs 2))
	  (d/dz (ref vs 3))
	  (dct (ref ws 0))
	  (dx (ref ws 1))
	  (dy (ref ws 2))
	  (dz (ref ws 3)))
      (define SR-vector-basis
	(down d/dct d/dx d/dy d/dz))

      (define SR-1form-basis
	(up dct dx dy dz))

      (define SR-basis
	(make-basis SR-vector-basis
		    SR-1form-basis))
      SR-basis)))

#|
(let ((g (flat-spacetime-metric))
      (b (spacetime-orthonormal-basis))
      (v1 (literal-vector-field 'V1 spacetime-rect))
      (v2 (literal-vector-field 'V2 spacetime-rect))
      (u (up (/ 'v^x :c) (/ 'v^y :c) (/ 'v^z :c)))
      (event ((point spacetime-rect) (up 'ct 'x 'y 'z))))
  (let ((LAMBDA ((boost b) u)))
    ((- (g v1 v2) (g (LAMBDA v1) (LAMBDA v2)))
     event)))
#| 0 |#

((4vector-boost
  (up (/ 'v^x :c) (/ 'v^y :c) (/ 'v^z :c)))
 (make-4vector 'ctp (up 'xp 'yp 'zp)))
#| mess |#

(set! *divide-out-terms* #f)

((4vector-boost (up (/ 'v^x :c) 0 0))
 (make-4vector (* :c 'tp) (up 'xp 'yp 'zp)))
#|
(up
 (/ (+ (* (expt :c 2) tp) (* v^x xp))
    (sqrt (+ (expt :c 2) (* -1 (expt v^x 2)))))
 (/ (+ (* :c tp v^x) (* :c xp)) (sqrt (+ (expt :c 2) (* -1 (expt v^x 2)))))
 yp
 zp)
|#

(determinant
 ((D (4vector-boost (up (/ 'v^x :c) (/ 'v^y :c) (/ 'v^z :c))))
  (make-4vector (* :c 'tp) (up 'xp 'yp 'zp))))
#| 1 |#

(let ((4xp (make-4vector (* :c 'tp) (up 'xp 'yp 'zp))))
  (let ((4x
	 ((4vector-boost (up (/ 'v^x :c) (/ 'v^y :c) (/ 'v^z :c)))
	  4xp)))
    (- (4vector-metric 4x 4x) (4vector-metric 4xp 4xp))))
#| 0 |#


(* ((D (4vector-boost
	(up (/ (* -1'v^x) :c) (/ (* -1 'v^y) :c) (/ (* -1'v^z) :c))))
    (make-4vector (* :c 'tp) (up 'xp 'yp 'zp)))
   ((D (4vector-boost (up (/ 'v^x :c) (/ 'v^y :c) (/ 'v^z :c))))
    (make-4vector (* :c 'tp) (up 'xp 'yp 'zp))))
#| (down (up 1 0 0 0) (up 0 1 0 0) (up 0 0 1 0) (up 0 0 0 1)) |#
|#

;;; this one works for zero v/c ...
;;; direction is a unit 3-vector, v/c is the speed, a number.

(define ((4vector-boost2 direction v/c) 4vector-prime)
  (let ((delta-ct-prime (4vector->ct 4vector-prime))
	(delta-x-prime (4vector->space 4vector-prime)))
    (let ((betasq (square v/c)))
      (let ((bx (dot-product direction delta-x-prime))
	    (gamma (/ 1 (sqrt (- 1 betasq)))))
	(let ((alpha (- gamma 1)))
	  (let ((delta-ct
		 (* gamma (+ delta-ct-prime (* bx v/c))))
		(delta-x
		 (+ (* gamma v/c direction delta-ct-prime)
		    delta-x-prime
		    (* alpha direction bx))))
	    (make-4vector delta-ct delta-x)))))))

#|
(let ((beta (up (/ 'v^x :c) (/ 'v^y :c) (/ 'v^z :c))))
  (- ((4vector-boost2 (/ beta (sqrt (square beta)))
		      (sqrt (square beta)))
      (up 'u0 'u1 'u2 'u3))
     ((4vector-boost beta) (up 'u0 'u1 'u2 'u3))))
#| (up 0 0 0 0) |#

(let ((beta (up (/ 'v^x :c) (/ 'v^y :c) (/ 'v^z :c))))
  (- ((4vector-boost2 (up 1 0 0) 0) (up 'u0 'u1 'u2 'u3))
     (up 'u0 'u1 'u2 'u3)))
#|(up 0 0 0 0) |#
|#

;;; direction is a unit 3-vector, rapidity is a number.

(define ((4vector-boost-rapidity direction) rapidity)
  (4vector-boost2 direction (tanh rapidity)))

(define ((extended-rotation R) xi-p)
  (make-4vector
   (4vector->ct xi-p)
   (R (4vector->space xi-p))))

#|
((D (extended-rotation (rotate-x 'theta))) (up 0 0 0 0))
#|
(down (up 1 0 0 0)
      (up 0 1 0 0)
      (up 0 0 (cos theta) (sin theta))
      (up 0 0 (* -1 (sin theta)) (cos theta)))
|#
|#

#|
;;; check rotation-boost identity

(let ((beta (up 'bx 'by 'bz))
      (xi (make-4vector 'ct (up 'x 'y 'z)))
      (R (compose
	  (rotate-x 'theta)
	  (rotate-y 'phi)
	  (rotate-z 'psi)))
      (R-inverse
       (compose
	(rotate-z (- 'psi))
	(rotate-y (- 'phi))
	(rotate-x (- 'theta)))))
  (- ((4vector-boost beta) xi)
     ((extended-rotation R-inverse)
      ((4vector-boost (R beta))
       ((extended-rotation R) xi)))))
#| (up 0 0 0 0) |#
|#

#|
;;; direction specified by colatitude and longitude, to get a unit
;;; vector for direction.

((rotate-z 'phi) ((rotate-y 'theta) (up 0 0 1)))
#|
(up (* (sin theta) (cos phi)) (* (sin theta) (sin phi)) (cos theta))
|#

;;; B(beta) = Re^-1 circ B(R(beta)) circ Re
;;; B(R(beta)) = Re circ B(beta) circ Re^-1
;;; in terms of rapidity a boost in the z-direction is

((D ((4vector-boost-rapidity (up 0 0 1)) 'rapidity)) (up 0 0 0 0))
#|
(down
 (up (cosh rapidity) 0 0 (sinh rapidity))
 (up 0 1 0 0)
 (up 0 0 1 0)
 (up (sinh rapidity) 0 0 (cosh rapidity)))
|#
;;; cheat: cosh^2-sinh^2=1
|#


(define ((boost-z rapidity) 4vector)
  (* (down
      (up (cosh rapidity) 0 0 (sinh rapidity))
      (up 0 1 0 0)
      (up 0 0 1 0)
      (up (sinh rapidity) 0 0 (cosh rapidity)))
     4vector))

(define (((boost-z->4vector-boost theta phi) rapidity) 4vector)
  (let ((R (compose (rotate-z phi) (rotate-y theta)))
	(R^-1 (compose (rotate-y (- theta)) (rotate-z (- phi)))))
    ((extended-rotation R)
     ((boost-z rapidity)
      ((extended-rotation R^-1) 4vector)))))

#|
((D ((boost-z->4vector-boost ':pi/2 0) 'rapidity)) (up 0 0 0 0))
#|
(down
 (up (cosh rapidity) (sinh rapidity) 0 0)
 (up (sinh rapidity) (cosh rapidity) 0 0)
 (up 0 0 1 0)
 (up 0 0 0 1))
|#
; after some hand simplification...
|#

;;; derivation of integrated program, below
;;; ((D ((boost-z->4vector-boost 'theta 'phi) 'rapidity)) (up 0 0 0 0))
;;; alpha = (* (cos phi) (sin theta))
;;; beta  = (* (sin theta) (sin phi))
;;; gamma = (cos theta)

(define (((rapidity-and-direction->boost direction) rapidity) 4vector)
  (let ((alpha (ref direction 0))
	(beta (ref direction 1))
	(gamma (ref direction 2)))
    (* (down (up (cosh rapidity)
		 (* alpha (sinh rapidity))
		 (* beta (sinh rapidity))
		 (* gamma (sinh rapidity)))
	     (up (* alpha (sinh rapidity))
		 (+ (* (square alpha) (cosh rapidity))
		    (- 1 (square alpha)))
		 (+ (* alpha beta (cosh rapidity))
		    (* -1 alpha beta))
		 (+ (* alpha gamma (cosh rapidity))
		    (* -1 alpha gamma)))
	     (up (* beta (sinh rapidity))
		 (+ (* alpha beta (cosh rapidity))
		    (* -1 alpha beta))
		 (+ (* (square beta) (cosh rapidity))
		    (- 1 (square beta)))
		 (+ (* beta gamma (cosh rapidity))
		    (* -1 beta gamma)))
	     (up (* gamma (sinh rapidity))
		 (+ (* alpha gamma (cosh rapidity))
		    (* -1 alpha gamma))
		 (+ (* beta gamma (cosh rapidity))
		    (* -1 beta gamma))
		 (+ (* (square gamma) (cosh rapidity))
		    (- 1 (square gamma)))))
       4vector)))
#|
(- (((boost-z->4vector-boost 'theta 'phi) 'rapidity)
    (up 'u0 'u1 'u2 'u3))
   (((rapidity-and-direction->boost
      (up (* (cos 'phi) (sin 'theta))
	  (* (sin 'phi) (sin 'theta))
	  (cos 'theta)))
     'rapidity)
    (up 'u0 'u1 'u2 'u3)))
#| (up 0 0 0 0) |#
|#

;;;;         Special-relativity frames.
;;; A frame is defined by a Lorentz transformation from a
;;; background 4-space frame.  To keep us from going nuts, an
;;; SR frame has a name, which it uses to label coordinates in
;;; its frame.  The background frame is called "the-ether".

(define (make-SR-frame name from-frame boost-direction v/c coordinate-origin)
  (define (coordinates->event coords)
    (assert (SR-coordinates? coords))
    (assert (eq? me (SR-frame coords)))
    (if (eq? me the-ether)
	coords
	(make-SR-coordinates the-ether
	  ((from-frame 'coords->event)
	   (make-SR-coordinates from-frame
	     (+ ((4vector-boost2 boost-direction v/c) coords)
		coordinate-origin))))))
  (define (event->coordinates event)
    (assert (SR-coordinates? event))
    (assert (eq? the-ether (SR-frame event)))
    (if (eq? me the-ether)
	event
	(make-SR-coordinates me
	  ((4vector-boost2 boost-direction (- v/c))
	   (- ((from-frame 'event->coords) event)
	      coordinate-origin)))))
  (define (me m)
    (case m
      ((dimension) 4)
      ((coords->event) coordinates->event)
      ((event->coords) event->coordinates)
      ((type) SR)
      ((name) name)
      ((from-frame) from-frame)
      ((boost-components) boost-components)
      ((coordinate-origin) coordinate-origin)
      (else (error "Unknown message: SR" name m))))
  me)

;;; Implementation of the coordinates uses a put/get table.

(define (make-SR-coordinates frame 4vector)
  (assert (vector? 4vector))
  (assert (fix:= (vector-length 4vector) 4))
  (eq-put! 4vector 'type 'SR-coordinates)
  (eq-put! 4vector 'frame frame)
  4vector)

(define (SR-coordinates? coords)
  (eq? (eq-get coords 'type) 'SR-coordinates))

(define (SR-frame coords)
  (eq-get coords 'frame))

(define (SR-name coords)
  ((eq-get coords 'frame) 'name))

(define the-ether
  (make-SR-frame 'the-ether #f #f #f #f))

#|
;;; Velocity addition formula

(define A
  (make-SR-frame 'andy the-ether
		 (up 1 0 0) (/ 'va :c)
		 #(0 0 0 0)))
(define B
  (make-SR-frame 'beth A
		 (up 1 0 0) (/ 'vb :c)
		 #(0 0 0 0)))

(set! *divide-out-terms* #f)

(let ((foo ((the-ether 'event->coords)
	    ((B 'coords->event)
	     (make-SR-coordinates B
				  (up (* :c 'tau) 0 0 0))))))
  (/ (ref foo 1) (/ (ref foo 0) :c)))
#| 
(/ (+ (* (expt :c 2) va) (* (expt :c 2) vb))
   (+ (expt :c 2) (* va vb)))
|#
|#

(define (add-velocities v1 v2)
  (/ (+ v1 v2) (+ 1 (* v1 v2))))

#|
;;; Experiments

(define A
  (make-SR-frame 'andy the-ether (up 1 0 0) 'va (up 0 0 0 0)))
(define B
  (make-SR-frame 'beth A         (up 1 0 0) 'vb (up 0 0 0 0)))
(define C
  (make-SR-frame 'cris A         (up 1 0 0) 'vc (up 0 'xc 0 0)))

((B 'event->coords)
 ((A 'coords->event)
  (make-SR-coordinates A (up 'ta 'xa 0 0))))
#|
(up (/ (+ (* -1 vb xa) ta) (sqrt (+ 1 (* -1 (expt vb 2)))))
    (/ (+ (* -1 ta vb) xa) (sqrt (+ 1 (* -1 (expt vb 2)))))
    0
    0)
|#

((C 'event->coords)
 ((A 'coords->event)
  (make-SR-coordinates A (up 'ta 'xa 0 0))))
#|
(up (/ (+ (* -1 vc xa) (* vc xc) ta)
       (sqrt (+ 1 (* -1 (expt vc 2)))))
    (/ (+ (* -1 ta vc) xa (* -1 xc))
       (sqrt (+ 1 (* -1 (expt vc 2)))))
    0
    0)
|#
|#

#|
;;; Twin paradox: green Schutz p28

;;; We will assume Artemis is moving with unknown velocity
;;; starting at an arbitrary point, and Diana is traveling
;;; with respect to him: v_diana=0.96c

(define A
  (let* ((bb (/ (up 'vx 'vy 'vz) :c))
	 (v/c (euclidean-norm bb))
	 (d (/ bb v/c)))
    (make-SR-frame 'artemis the-ether
		   d v/c
		   (up (* :c 't0) 'x0 'y0 'z0))))


;;; In Artemis coordinates, the path of Artemis is:

(define (A-path ta)
  (make-SR-coordinates A (up (* :c ta) 0 0 0)))


;;; Artemis's positions in Artemis coordinates

(define A-origin-A-coords (A-path 0))

A-origin-A-coords
#| (up 0 0 0 0) |#

;;; After 25 years of Artemis time, Diana turns around.

(define A-turning-A-coords (A-path 25))

A-turning-A-coords
#| (up (* 25 :c) 0 0 0) |#


(define A-meeting-A-coords (A-path 50))

A-meeting-A-coords
#| (up (* 50 :c) 0 0 0) |#

(4vector-interval (- A-meeting-A-coords A-origin-A-coords))
#| (* +50i :c) |#

;;; Artemis ages 50 years.  (Imaginary, because time-like.)

;;; In Artemis coordinates, the path of Diana has two
;;; segments, outgoing and incoming:

(define (D-outgoing-A-coords ta)
  (make-SR-coordinates A (up (* :c ta) (* 96/100 :c ta) 0 0)))

(define (D-incoming-A-coords ta)
  (let ((t2 (- ta 25))
	(start-of-seg2 (D-outgoing-A-coords 25)))
    (make-SR-coordinates A
      (+ (up (* :c t2) (* -96/100 :c t2) 0 0)
	 start-of-seg2))))


;;; Diana's positions in Artemis coordinates

(define D-origin-A-coords
  (D-outgoing-A-coords 0))

;;; They start out together
D-origin-A-coords
#| (up 0 0 0 0) |#


(define D-turning-A-coords
  (D-outgoing-A-coords 25))

;;; After 25 Artemis years, Diana is 24 light-years away.
D-turning-A-coords
#| (up (* 25 :c) (* 24 :c) 0 0) |#


(define D-meeting-A-coords
  (D-incoming-A-coords 50))

;;; After 50 Artemis years, Diana is back with Artemis.
D-meeting-A-coords
#| (up (* 50 :c) 0 0 0) |#

;;; But in that time, Diana has aged only 14 years!

(+ (4vector-interval (- D-meeting-A-coords D-turning-A-coords))
   (4vector-interval (- D-turning-A-coords D-origin-A-coords)))
#| (* +14i :c) |#

;;; Now, what does this look like from Diana's point of view?

;;; Diana's frames for the two parts of her trip are boosts from
;;; Artemis's frame.

(define D1
  (make-SR-frame 'diana1 A (up 1 0 0) 96/100 A-origin-A-coords))

(define D2
  (make-SR-frame 'diana2 A (up 1 0 0) -96/100 D-turning-A-coords))


;;; In Diana's coordinates, during her outgoing trip Artemis is
;;; recessing symmetrically.

(define (A-outgoing-D1-coords td)
  (make-SR-coordinates D1 (up (* :c td) (* -96/100 :c td) 0 0)))

(define A-origin-D1-coords (A-outgoing-D1-coords 0))

A-origin-D1-coords
#| (up 0 0 0 0) |#

(define A-turning-D1-coords (A-outgoing-D1-coords 7))

A-turning-D1-coords
#| (up (* 7 :c) (* -168/25 :c) 0 0) |#

;;; In 7 years of Diana time, Artemis has aged only about 2 years.

(4vector-interval (- A-turning-D1-coords A-origin-D1-coords))
#| (* +49/25i :c) |#

(/ 49. 25)
#| 1.96 |#

;;; If she continued for 25 years, Artemis would age 7 years... 
;;; This is the symmetry of time dialation.

(4vector-interval (- (A-outgoing-D1-coords 25) A-origin-D1-coords))
#| (* +7i :c) |# 

;;; But Diana turns around at 7 years, getting a new viewpoint!

(define A-origin-D2-coords
  ((compose (D2 'event->coords) (D1 'coords->event))
   A-origin-D1-coords))

A-origin-D2-coords
#| (up (* -1201/7 :c) (* -1200/7 :c) 0 0) |#

(define A-turning-D2-coords
  ((compose (D2 'event->coords) (D1 'coords->event))
   A-turning-D1-coords))

A-turning-D2-coords
#| (up (* -1152/7 :c) (* -28824/175 :c) 0 0) |#

(4vector-interval (- A-turning-D2-coords A-origin-D2-coords))
#| (* +49/25i :c) |#		;The interval is invariant.

(define A-meeting-D2-coords
  ((compose (D2 'event->coords) (A 'coords->event))
   D-meeting-A-coords))

A-meeting-D2-coords
#| (up (* 7 :c) 0 0 0) |#
;;; Right place, 7 years from the origin of the D2 frame.

(4vector-interval (- A-meeting-D2-coords A-origin-D2-coords))
#| (* +50i :c) |# 
;;; Apparently, Artemis has aged for 50 years...

(4vector-interval (- A-meeting-D2-coords A-turning-D2-coords))
#| (* +1201/25i :c) |#
;;; All his ageing appears to have happened during Diana's return trip!

(/ 1201. 25)
#| 48.04 |#

;;; In the last episode we saw that Diana's perception of Artemis's
;;; ageing was asymmetrical: almost all of his ageing appears to have
;;; happened on the return journey.  But the incoming trip is just
;;; like the outgoing trip, so this makes no sense.  However the two
;;; measurements are not symmetrical because the start of the incoming
;;; trip was computed from the end of the outgoing trip, whereas the
;;; start of the outgoing trip was from the initial condition, where
;;; Artemis and Diana are coincident.  The only other time they are
;;; coincident is at the reunion.  What does the turning point look 
;;; like from there?

;;;  In terms of D2 time starting at the meeting point.  The path
;;;  taken by Artemis, until he rejoins Diana, is:

(define (A-incoming-D2-coords td)
  (make-SR-coordinates D2 (up (* :c (+ td 7)) (* 96/100 td :c) 0 0)))

;;; So, seven years before the reunion, Diana thought that Artemis was
;;; at the turning point.

(define A-turning-D2-coords-rev
  (A-incoming-D2-coords -7))

A-turning-D2-coords-rev
#| (up 0 (* -168/25 :c) 0 0) |#

;;; So Diana thinks that, over the last segment, Artemis only aged
;;; about 2 years, just like in the first segment.

(4vector-interval (- A-meeting-D2-coords A-turning-D2-coords-rev))
#| (* +49/25i :c) |#

;;; So the wierdness has to do with the different views of the turning
;;; point.

;;; *** Now look at this!

(4vector-interval
 (- A-turning-D2-coords-rev A-turning-D2-coords))
#| (* +1152/25i :c) |#

(+ 49/25 1152/25 49/25)
#| 50 |#

;;; So... this is the real place that the discrepancy is hiding.
|#

;;;; A differential geometry point-of-view

(define (register-SR-frame frame)
  (let ((name (frame 'name)))
    (attach-coordinate-system name 'origin R^n
      (lambda (manifold)
	(define (me m)
	  (case m
	    ((check-coordinates)
	     (lambda (coords)
	       (and (up? coords) (fix:= (s:dimension coords) 4))))
	    ((coords->point)
	     (lambda (coords)
	       (assert ((me 'check-coordinates) coords)
		       "Bad coordinates: SR-frame" (list name coords))
	       (make-manifold-point ((frame 'coords->event) coords)
				    manifold
				    me
				    coords)))
	    ((check-point)
	     (lambda (point)
	       (my-manifold-point? point manifold)))
	    ((point->coords)
	     (lambda (point)
	       (assert ((me 'check-point) point)
		       "Bad point: SR-frame" (list name point))
	       (get-coordinates point me
		 (lambda ()
		   (let ((prep (manifold-point-representation point)))
		     (assert (and (up? prep)
				  (fix:= (s:dimension prep)
					 (manifold 'embedding-dimension)))
			     "Bad point: SR-frame" (list name point))
		     ((frame 'event->coords) prep))))))
	    (else (error "Bad message: SR-frame" name m))))
	me))))

#|
(register-SR-frame the-ether)
(register-SR-frame A)
(register-SR-frame D1)
(register-SR-frame D2)

(define Artemis-frame
  (coordinate-system-at 'artemis 'origin spacetime))

(define Diana1-frame
  (coordinate-system-at 'diana1 'origin spacetime))

(define Diana2-frame
  (coordinate-system-at 'diana2 'origin spacetime))


(define-coordinates tau_a the-real-line)
(define-coordinates tau_d the-real-line)

(define dg-D-outgoing-A-coords
  (compose (point Artemis-frame) D-outgoing-A-coords (chart the-real-line)))

(define dg-D-incoming-A-coords
  (compose (point Artemis-frame) D-incoming-A-coords (chart the-real-line)))

;;; After 25 years of Artemis's time, Diana gets to the
;;; turning point event.

(define origin
  (dg-D-outgoing-A-coords ((the-real-line '->point) 0)))

((chart Artemis-frame) origin)
#| (up 0 0 0 0) |#

(define turning-event
  (dg-D-outgoing-A-coords ((the-real-line '->point) 25)))

((chart Artemis-frame) turning-event)
#| (up (* 25 :c) (* 24 :c) 0 0) |#

(SR-name ((chart Artemis-frame) turning-event))
#| artemis |#

((chart Diana1-frame) turning-event)
#| (up (* 7 :c) 0 0 0) |#

((chart Diana2-frame) turning-event)
#| (up (* 7 :c) 0 0 0) |#

;;; The tangent vector to D-outgoing-A-coords 
;;; (over the map D-outgoing-A-coords)

(define T1 ((differential dg-D-outgoing-A-coords) d/dtau_a))

;;; Note: the vector field transported is an Artemis vector

;;; The turning point should be the integral of the tangent over 25
;;; years.  Because D-outgoing-A-coords is a straight line, this can
;;; be evaluated at any time.

(- (- ((chart Artemis-frame) turning-event)
      ((chart Artemis-frame) origin))
   (((* T1 25) (chart Artemis-frame))
    ((point the-real-line) 't_a)))
#| (up 0 0 0 0) |#

;;; We can also examine the tangent vector parameterized by
;;; Diana's proper time.

(define dg-seg1_d
  (compose (point Artemis-frame) seg1_d (chart the-real-line)))

(define T1_d ((differential dg-seg1_d) d/dtau_d))

(((* T1_d 7) (chart Diana1-frame))
 ((point the-real-line) 't_d))
#| (up (* 7 :c) 0 0 0) |#

(((* T1_d 7) (chart Artemis-frame))
 ((point the-real-line) 't_d))
#| (up (* 25 :c) (* 24 :c) 0 0) |#

(4vector-interval
 (((* T1_d 7) (chart Artemis-frame))
  ((point the-real-line) 't_d)))
#| (* +7i :c) |#

(sqrt
 (((metric-over-map dg-seg1_d (flat-spacetime-metric))
   T1_d T1_d)
  ((point the-real-line) 't)))
#| (* +i :c) |#

;;; The 4-velocity of a time-like path, parameterized by its
;;; proper time, is always :c.

;;; Now, from Diana's point of view:

((chart Diana1-frame) turning-event)
#| (up (* 7 :c) 0 0 0) |#
;;; So, Diana has only aged 7 years.

;;; Indeed, Diana's proper time is the 4vector-interval of her
;;; world line, which is an integral of the metric measuring
;;; the tangent vector over her path.

(* (sqrt
     (- (((metric-over-map dg-D-outgoing-A-coords
			   (flat-spacetime-metric))
	  T1 T1)
	 ((the-real-line '->point) 't))))
    25)
#| (* 7 :c) |#
;;; measured in length!

;;; Now, on the second segment, Diana reverses course.  

(define T2 ((differential dg-D-incoming-A-coords) d/dtau_a))

(+ (* (sqrt
       (- (((metric-over-map dg-D-outgoing-A-coords
			     (flat-spacetime-metric))
	    T1 T1)
	   ((the-real-line '->point) 5))))
      25)
   (* (sqrt
       (- (((metric-over-map dg-D-incoming-A-coords
			     (flat-spacetime-metric))
	    T2 T2)
	   ((the-real-line '->point) 5))))
      25))
#| (* 14 :c) |#

;;; So Diana has only aged 14 years while Artemis has aged
;;; 50 years.  This is 1/2 of the computation.  We must also
;;; compute Artemis's proper aging from Diana's coordinates.

;;; Revert to numerical constant with units:

(numerical-constants #t (list (get-constant-data ':c)))
|#
