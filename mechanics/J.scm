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


(define (J-func DH)
  (->H-state (ref DH 0) 
	     (ref DH 2)
	     (- (ref DH 1))))

#|
(print-expression
 ((D J-func)
  (compatible-shape
   (->H-state 't
	      (coordinate-tuple 'x 'y)
	      (momentum-tuple 'p_x 'p_y)))))
(up (up 1 (up 0 0) (down 0 0))
    (up (up 0 (up 0 0) (down -1 0)) (up 0 (up 0 0) (down 0 -1)))
    (down (up 0 (up 1 0) (down 0 0)) (up 0 (up 0 1) (down 0 0))))    
|#

(define ((phase-space-derivative-2 H) s)
  (let ((J ((D J-func) (compatible-shape s))))
    (* J ((D H) s))))

#|
(print-expression
 ((phase-space-derivative-2 (H-harmonic 'm 'k))
  (->H-state 't
	     (coordinate-tuple 'x 'y)
	     (momentum-tuple 'p_x 'p_y))))
(up 0 (up (/ p_x m) (/ p_y m)) (down (* -1 k x) (* -1 k y)))
|#

(define (linear-function->multiplier F argument)
  ((D F) argument))

(define (compatible-shape s)
  (linear-function->multiplier (lambda (x) 0) s))

(define ((time-independent-c? C) s)
  (let ((s* (compatible-shape s)))
    (let ((J (linear-function->multiplier J-func s*)))
      (- J 
	 (* ((D C) s)
	    (* J
	       ((adjoint s*) ((D C) s))))))))

(define ((adjoint s) A)
  (linear-function->multiplier (adjoint-function A) s))

(define ((adjoint-function A) p) (* p A))

#|
(define (T v)
  (* (down (up 'a 'c) (up 'b 'd)) v))

(pe (T (up 'x 'y)))
(up (+ (* a x) (* b y)) (+ (* c x) (* d y)))

(pe (* (* (down 'p_x 'p_y) ((D T) (up 'x 'y))) (up 'v_x 'v_y)))
(+ (* a p_x v_x) (* b p_x v_y) (* c p_y v_x) (* d p_y v_y))


(pe (* (down 'p_x 'p_y) (* ((D T) (up 'x 'y)) (up 'v_x 'v_y))))
(+ (* a p_x v_x) (* b p_x v_y) (* c p_y v_x) (* d p_y v_y))

(pe (* (* ((adjoint (down 'p_x 'p_y)) ((D T) (up 'x 'y)))
	  (down 'p_x 'p_y))
       (up 'v_x 'v_y)))
(+ (* a p_x v_x) (* b p_x v_y) (* c p_y v_x) (* d p_y v_y))

;;; But strangely enough...
(pe (* (* (down 'p_x 'p_y)
	  ((adjoint (down 'p_x 'p_y)) ((D T) (up 'x 'y))))
       (up 'v_x 'v_y)))
(+ (* a p_x v_x) (* b p_x v_y) (* c p_y v_x) (* d p_y v_y))
|#


#|
(print-expression
 ((time-independent-c? (F->CT p->r))
  (->H-state 't
	     (coordinate-tuple 'r 'varphi)
	     (momentum-tuple 'p_r 'p_varphi))))
(up (up 0 (up 0 0) (down 0 0))
    (up (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0)))
    (down (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0))))

(print-expression
 ((time-independent-c? a-non-canonical-transform)
  (->H-state 't 'theta 'p)))
(up (up 0 0 0) (up 0 0 (+ -1 p)) (up 0 (+ 1 (* -1 p)) 0))

(print-expression
 ((time-independent-c? (theta-I->x-px 'alpha))
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
 ((time-independent-c? Cmix)
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
 ((time-independent-c? Cmix2)
  a-state))
(up (up 0 (up 0 0) (down 0 0))
    (up (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0)))
    (down (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0))))

|#

;;;----------------------------------------------------------------
;;; Poisson brackets in terms of J
#|

(define ((PB f g) s)
  (* ((D f) s) (J-func ((D g) s))))

(define (shape x) x)

(define ((PB f g) s)
  (let ((J ((D J-func) (shape ((D g) s)))))
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
(up 0 0)

(pe ((PB momentum (H-harmonic 'm 'k))
     a-state))
(down (* -1 k x) (* -1 k y))

(pe ((PB coordinate (H-harmonic 'm 'k))
     a-state))
(up (/ p_x m) (/ p_y m))


|#

;;;----------------------------------------------------------------
;;; generating functions

;;; identity
(define (F2-identity q p-prime t)
  (+ (* (ref q 0) (ref p-prime 0))
     (* (ref q 1) (ref p-prime 1))))

;;; flip q & p
(define (F1-flip q q-prime t)
  (+ (* (ref q 0) (ref q-prime 0))
     (* (ref q 1) (ref q-prime 1))))

(define ((F2->C-inv F2) H-state)
  (let ((t (time H-state))
	(q (coordinate H-state))
	(p (momentum H-state)))
    (let ((p-func
	   (lambda (p-prime)
	     (((partial 0) F2) q p-prime t))))
      (let ((p-prime ((linear-inverse p-func) p)))
	(let ((q-prime (((partial 1) F2) q p-prime t)))
	  (->H-state t q-prime p-prime))))))
    
;;; x = f(x') is linear
(define ((linear-inverse f) x)
  (let ((b (f (zero-like x)))
	(a ((D f) (zero-like x))))
    (/ (- x b) a)))  

#|
(print-expression ((F2->C-inv F2-identity) a-state))
(up t (up x y) (down p_x p_y))
 
(print-expression
 ((time-independent-c? (F2->C-inv F2-identity))
  a-state))
(up (up 0 (up 0 0) (down 0 0))
    (up (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0)))
    (down (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0))))
|#

