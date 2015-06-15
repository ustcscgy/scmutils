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


(define (ultra-flatten s)
  (if (structure? s)
      (apply append
	     (map ultra-flatten
		  (vector->list (s:->vector s))))
      (list s)))

#|

(ultra-flatten (up 1 2 'a (down 3 4) (up (down 'c 'd) 'e)))
;Value 16: (1 2 a 3 4 c d e)

;;; similar to s:fringe (reverse)

|#

(define (number-elements s)
  (if (structure? s)
      (apply +
	     (map number-elements
		  (vector->list (s:->vector s))))
      1))

#|

(number-elements (up 1 2 'a (down 3 4) (up (down 'c 'd) 'e)))
;Value: 8

|#

(define (ultra-unflatten shape list)
  (if (structure? shape)
      (let lp ((s '()) (i 0) (list list))
	(if (< i (s:length shape))
	    (lp (cons (ultra-unflatten (ref shape i) list) s)
		(+ i 1)
		(list-tail list (number-elements (ref shape i))))
	    (s:structure (s:same shape) (list->vector (reverse s)))
	    ))
      (car list)))

#|

(ultra-unflatten
      (up 'x 'x 'x (down 'x 'x) (up (down 'x 'x) 'x))
      (list 1 2 'a 3 4 'c 'd 'e))
;Value 30: #(1 2 a (*down* #(3 4)) #((*down* #(c d)) e))

(pe (- 
     (ultra-unflatten
      (up 'x 'x 'x (down 'x 'x) (up (down 'x 'x) 'x))
      (list 1 2 'a 3 4 'c 'd 'e))
     (up 1 2 'a (down 3 4) (up (down 'c 'd) 'e))))
(up 0 0 0 (down 0 0) (up (down 0 0) 0))

|#


#|
(define vs
  (velocity-tuple
   (velocity-tuple 'vx1 'vy1)
   (velocity-tuple 'vx2 'vy2)))

(define (L1 vs)
  (let ((v1 (ref vs 0))
	(v2 (ref vs 1)))
    (+ (* 1/2 'm1 (square v1))
       (* 1/2 'm2 (square v2)))))

(pe (((expt D 2) L1) vs))
(down (down (down (down m1 0) (down 0 0)) (down (down 0 m1) (down 0 0)))
      (down (down (down 0 0) (down m2 0)) (down (down 0 0) (down 0 m2))))
|#

(define (s->m ls ms rs)
  (assert (numerical-quantity? (* ls (* ms rs))))
  (let ((lv (vector->row-matrix (list->vector (ultra-flatten ls))))
	(rv (vector->column-matrix (list->vector (ultra-flatten rs)))))
    (let ((nrows (m:num-cols lv)) (ncols (m:num-rows rv)))
      (m:generate nrows ncols
		  (lambda (i j)
		    (* (ultra-unflatten
			ls
			(vector->list
			 (v:make-basis-unit nrows i)))
		       (* ms
			  (ultra-unflatten
			   rs
			   (vector->list
			    (v:make-basis-unit ncols j))))))))))
#|
(pe (s->m vs (((expt D 2) L1) vs) vs))
(matrix-by-rows (list m1 0 0 0)
                (list 0 m1 0 0)
                (list 0 0 m2 0)
                (list 0 0 0 m2))
|#


(define (m->s ls m rs)
  (let ((ncols (m:num-cols m))
	(col-shape (compatible-shape ls)))
    (ultra-unflatten (compatible-shape rs)
		     (let lp ((j 0))
		       (if (= j ncols)
			   '()
			   (let ((colj (m:nth-col m j)))
			     (cons (ultra-unflatten col-shape
						    (vector->list colj))
				   (lp (+ 1 j)))))))))

#|
(pe (m->s vs (s->m vs (((expt D 2) L1) vs) vs) vs))
(down (down (down (down m1 0) (down 0 0)) (down (down 0 m1) (down 0 0)))
      (down (down (down 0 0) (down m2 0)) (down (down 0 0) (down 0 m2))))
|#


#|

(define ((rot-z angle) vs)
  (define ((rotate angle) v)
    (let ((vx (ref v 0))
	  (vy (ref v 1)))
      (velocity-tuple (- (* (cos angle) vx)
			 (* (sin angle) vy))
		      (+ (* (sin angle) vx)
			 (* (cos angle) vy)))))
  (s:generate (s:length vs) (s:same vs)
	      (lambda (i) ((rotate angle) (ref vs i)))))


(pe ((rot-z 'a) vs))
(up
 (up (+ (* vx1 (cos a)) (* -1 vy1 (sin a)))
     (+ (* vx1 (sin a)) (* vy1 (cos a))))
 (up (+ (* vx2 (cos a)) (* -1 vy2 (sin a)))
     (+ (* vx2 (sin a)) (* vy2 (cos a)))))


(define R
  (linear-function->multiplier (rot-z 'a) vs))

(pe R)
(down
 (down (up (up (cos a) (sin a)) (up 0 0))
       (up (up (* -1 (sin a)) (cos a)) (up 0 0)))
 (down (up (up 0 0) (up (cos a) (sin a)))
       (up (up 0 0) (up (* -1 (sin a)) (cos a)))))

(pe (* R vs))
(up
 (up (+ (* vx1 (cos a)) (* -1 vy1 (sin a)))
     (+ (* vx1 (sin a)) (* vy1 (cos a))))
 (up (+ (* vx2 (cos a)) (* -1 vy2 (sin a)))
     (+ (* vx2 (sin a)) (* vy2 (cos a)))))


(define M
  (s->m (compatible-shape vs)
	R 
	vs))

(pe M)
(matrix-by-rows (list (cos a) (* -1 (sin a)) 0 0)
                (list (sin a) (cos a) 0 0)
                (list 0 0 (cos a) (* -1 (sin a)))
                (list 0 0 (sin a) (cos a)))

(pe (m->s (compatible-shape vs)
	  M
	  vs))
(down
 (down (up (up (cos a) (sin a)) (up 0 0))
       (up (up (* -1 (sin a)) (cos a)) (up 0 0)))
 (down (up (up 0 0) (up (cos a) (sin a)))
       (up (up 0 0) (up (* -1 (sin a)) (cos a)))))
|#

(define (new-transpose ls ms rs)
  (m->s ls
	(m:transpose
	 (s->m ls ms rs))
	rs))

(define (new-inverse ls ms rs)
  (m->s ls
	(m:invert
	 (s->m ls ms rs))
	rs))

#|

(pe (new-transpose



(define (structure-inverse ls ms rs)
  
...

|#