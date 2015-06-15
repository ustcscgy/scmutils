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

;;; H(t (q0 q1) [p0 p1]) ==> R
;;; (DH(t (q0 q1) [p0 p1]))*(dt (dq0 dq1) [dp0 dp1]) ==> R
;;; Thus (DH(t (q0 q1) [p0 p1])) = [D0H, [D10H D11H] (D20H D21H)]

;;; Hamilton's equations are
;;; D(identity (q0 q1) [p0 p1]) = J-func(DH(t (q0 q1) [p0 p1]))
;;; J-func[D0H, [D10H D11H] (D20H D21H)] = (1 (dq0 dq1) -[dp0 dp1])

;;; Thus J-func may be defined as follows

(define (J-func DH)
  (->H-state (ref DH 0) 
	     (ref DH 2)
	     (- (ref DH 1))))

#|
;;; J-func is a linear function.
;;; The tuple representation of J-func (as a multiplier) is thus:

(print-expression
 ((D J-func)
  (->H-state 't
	     (coordinate-tuple 'x 'y)
	     (momentum-tuple 'p_x 'p_y))))
(down (up 1 (down 0 0) (up 0 0))
      (down (up 0 (down 0 0) (up -1 0)) (up 0 (down 0 0) (up 0 -1)))
      (up (up 0 (down 1 0) (up 0 0)) (up 0 (down 0 1) (up 0 0))))    

;;; Suppose we name the tuple multiplier J.  (This is not a matrix.)
|#


;;; So we can implement the phase-space derivative using this multiplier.
;;;   Note that compatible-shape gives an object that when multiplied by
;;;   its argument produces a number.

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

;;; This is indeed of the form required to make a D(t,q,p).
|#

;;; To test whether a purported transformation is canonical.


;;; Functional method only. 

(define ((time-independent-c? C) s)
  (let ((s* (compatible-shape s)))
    (let ((J ((D J-func) s*)))
      (- J 
	 (* ((D C) s)
	    (* J
	       ((multiplicative-transpose s*) ((D C) s))))))))

(define ((multiplicative-transpose s) A)
  ((D (transpose-function A)) s))

(define ((transpose-function A) p) (* p A))

#|
;;; This idea apparently works.

(print-expression
 ((time-independent-c? (F->CT p->r))
  (->H-state 't
	     (coordinate-tuple 'r 'varphi)
	     (momentum-tuple 'p_r 'p_varphi))))
(up (up 0 (up 0 0) (down 0 0))
    (up (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0)))
    (down (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0))))



(define (a-non-canonical-transform Istate)
  (let ((t (time Istate))
        (theta (coordinate Istate))
	(p (momentum Istate)))
    (let ((x (* p (sin theta)))
	  (p_x (* p (cos theta))))
      (->H-state t x p_x))))

(print-expression
 ((time-independent-c? a-non-canonical-transform)
  (->H-state 't 'theta 'p)))
(up (up 0 0 0) (up 0 0 (+ -1 p)) (up 0 (+ 1 (* -1 p)) 0))
|#

#|
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

(define a-state
  (->H-state 't 
	     (coordinate-tuple 'x 'y)
	     (momentum-tuple 'p_x 'p_y)))

(print-expression ((time-independent-c? Cmix) a-state))
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

(print-expression ((time-independent-c? Cmix2) a-state))
(up (up 0 (up 0 0) (down 0 0))
    (up (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0)))
    (down (up 0 (up 0 0) (down 0 0)) (up 0 (up 0 0) (down 0 0))))

|#

#|
(pe (s:transpose (down 0 0 0) (down (up 'a 'b 'c) (up 'c 'd 'e)) (up 0 0)))
(down (up a c) (up b d) (up c e))
|#
#|
(define ((time-independent-c? C) s)
  (let ((s* (compatible-shape s)))
    (let ((J ((D J-func) s*)))
      (- J 
	 (* ((D C) s)
	    (* J
	       (flip-indices (s:transpose s* ((D C) s) s))))))))

|#

#|
;;; A matrix implementation of time-independent-c?

(define ((time-independent-c? C) s)
  (let ((s* (compatible-shape s)))
    (let ((J (s->m (compatible-shape (J-func s*))
		   ((D J-func) s*)
		   s*))
	  (DC (s->m s* ((D C) s) s)))
      (- J (* DC (* J (m:transpose DC)))))))

(print-expression
 ((time-independent-c? (F->CT p->r))
  (->H-state 't
	     (coordinate-tuple 'r 'varphi)
	     (momentum-tuple 'p_r 'p_varphi))))
(matrix-by-rows (list 0 0 0 0 0)
		(list 0 0 0 0 0)
		(list 0 0 0 0 0)
		(list 0 0 0 0 0)
		(list 0 0 0 0 0))

(print-expression
 ((time-independent-c? a-non-canonical-transform)
  (->H-state 't 'theta 'p)))
(matrix-by-rows (list 0 0 0) (list 0 0 (+ 1 (* -1 p))) (list 0 (+ -1 p) 0))

(print-expression
 ((time-independent-c? (theta-I->x-px 'alpha))
  (->H-state 't 'a 'I)))
(matrix-by-rows (list 0 0 0) (list 0 0 0) (list 0 0 0))
|#

#|
(print-expression ((time-independent-c? Cmix) a-state))
(matrix-by-rows (list 0 0 0 0 0)
		(list 0 0 0 0 0)
		(list 0 0 0 0 0)
		(list 0 0 0 0 0)
		(list 0 0 0 0 0))

(print-expression ((time-independent-c? Cmix2) a-state))
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

(print-expression ((time-independent-c? (C 'm_0 'm_1)) b-state))
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

;;; Poisson brackets, for time-independent scalar functions of phase space.

(define ((PB f g) s)
  (* ((D f) s) (J-func ((D g) s))))

#|
(pe 
 (- ((Poisson-bracket (H-harmonic 'm 'k)
		      ((component 0) coordinate)) 
     a-state)
    ((PB (H-harmonic 'm 'k)
	 (compose (component 0) coordinate))
     a-state)
    ))
0


;;; Not defined for structures, but it worked here, because 

(pe 
 (- ((Poisson-bracket (H-harmonic 'm 'k)
		      coordinate) 
     a-state)
    ((PB (H-harmonic 'm 'k)
	 coordinate)
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

#|
;;; Exploration of the 2-form

;;; Matrix version

(let* ((C (F->CT p->r))
       (s 
	(->H-state 't
		   (coordinate-tuple 'r 'varphi)
		   (momentum-tuple 'p_r 'p_varphi)))
       (s* (compatible-shape s))
       (J (s->m (compatible-shape (J-func s*))
		((D J-func) s*)
		s*))
       (DC (s->m s* ((D C) s) s))
       (ds1
	(s->m s*
	      (->H-state 't_1
			 (coordinate-tuple 'dr_1 'dvarphi_1)
			 (momentum-tuple 'dp_r_1 'dp_varphi_1))
	      1))
       (ds2
	(s->m s*
	      (->H-state 't_2
			 (coordinate-tuple 'dr_2 'dvarphi_2)
			 (momentum-tuple 'dp_r_2 'dp_varphi_2))
	      1)))
  (pe (- (* (m:transpose ds1) J ds2)
	 (* (m:transpose (* DC ds1)) J (* DC ds2)))))
(matrix-by-rows (list 0))

;;; ((D J-func) s*) is a 2-up version of the tensor

(let* ((C (F->CT p->r))
       (s 
	(->H-state 't
		   (coordinate-tuple 'r 'varphi)
		   (momentum-tuple 'p_r 'p_varphi)))
       (s* (compatible-shape s))
       (DC ((D C) s))
       (ds1
	(->H-state 't_1
		   (coordinate-tuple 'dr_1 'dvarphi_1)
		   (momentum-tuple 'dp_r_1 'dp_varphi_1)))
       (ds2
	(->H-state 't_2
		   (coordinate-tuple 'dr_2 'dvarphi_2)
		   (momentum-tuple 'dp_r_2 'dp_varphi_2))))
  (pe (- (* (flip-indices ds1) (J-func (flip-indices ds2)))
	 (* (flip-indices (* DC ds1))
	    (J-func (flip-indices (* DC ds2)))))))
0


(let* ((C (F->CT p->r))
       (s 
	(->H-state 't
		   (coordinate-tuple 'r 'varphi)
		   (momentum-tuple 'p_r 'p_varphi)))
       (s* (compatible-shape s))
       (DC ((D C) s))
       (J ((D J-func) s*))
       (ds1
	(->H-state 't_1
		   (coordinate-tuple 'dr_1 'dvarphi_1)
		   (momentum-tuple 'dp_r_1 'dp_varphi_1)))
       (ds2
	(->H-state 't_2
		   (coordinate-tuple 'dr_2 'dvarphi_2)
		   (momentum-tuple 'dp_r_2 'dp_varphi_2))))
  (pe (- (* (flip-indices ds1) (* J (flip-indices ds2)))
	 (* (flip-indices (* DC ds1))
	    (* J (flip-indices (* DC ds2)))))))
0

;;; We can make a 2-down version

(let* ((C (F->CT p->r))
       (s 
	(->H-state 't
		   (coordinate-tuple 'r 'varphi)
		   (momentum-tuple 'p_r 'p_varphi)))
       (s* (compatible-shape s))
       (DC ((D C) s))
       (J (flip-indices ((D J-func) s*)))
       (ds1
	(->H-state 't_1
		   (coordinate-tuple 'dr_1 'dvarphi_1)
		   (momentum-tuple 'dp_r_1 'dp_varphi_1)))
       (ds2
	(->H-state 't_2
		   (coordinate-tuple 'dr_2 'dvarphi_2)
		   (momentum-tuple 'dp_r_2 'dp_varphi_2))))
  (pe (- (* ds1 (* J ds2))
	 (* (* DC ds1) (* J (* DC ds2))))))
0


;;; This is NOT vacuous

(let* ((C (F->CT p->r))
       (s 
	(->H-state 't
		   (coordinate-tuple 'r 'varphi)
		   (momentum-tuple 'p_r 'p_varphi)))
       (s* (compatible-shape s))
       (DC ((D C) s))
       (J (flip-indices ((D J-func) s*)))
       (ds1
	(->H-state 't_1
		   (coordinate-tuple 'dr_1 'dvarphi_1)
		   (momentum-tuple 'dp_r_1 'dp_varphi_1)))
       (ds2
	(->H-state 't_2
		   (coordinate-tuple 'dr_2 'dvarphi_2)
		   (momentum-tuple 'dp_r_2 'dp_varphi_2))))
  (pe (* ds1 (* J ds2))))
(+ (* -1 dp_r_1 dr_2) (* dp_r_2 dr_1) (* -1 dp_varphi_1 dvarphi_2) (* dp_varphi_2 dvarphi_1) (* t_1 t_2))

;;; Now for something really hard ...

(define db-state
  (->H-state 'dt
	     (coordinate-tuple (coordinate-tuple 'dx_1 'dy_1)
			       (coordinate-tuple 'dx_2 'dy_2))
	     (momentum-tuple (momentum-tuple 'dp_x_1 'dp_y_1)
			     (momentum-tuple 'dp_x_2 'dp_y_2))))


(define cb-state
  (->H-state 'ct
	     (coordinate-tuple (coordinate-tuple 'cx_1 'cy_1)
			       (coordinate-tuple 'cx_2 'cy_2))
	     (momentum-tuple (momentum-tuple 'cp_x_1 'cp_y_1)
			     (momentum-tuple 'cp_x_2 'cp_y_2))))


(let* ((s b-state)
       (s* (compatible-shape s))
       (DC ((D (C 'm1 'm2)) s))
       (J (flip-indices ((D J-func) s*)))
       (ds1 db-state)
       (ds2 cb-state))
  (pe (- (* ds1 J ds2)
	 (* (* DC ds1) J (* DC ds2)))))
0
|#
