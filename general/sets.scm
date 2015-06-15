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

;;;; Sets -- Implementation as ordered lists of symbols

(declare (usual-integrations))

;;; The arguments determine the type of the elements.

(define (make-sets-package set-equal-elements? set-less-elements?)

  (define the-empty-set '())
  (define set-empty? null?)

  (define set-first car)
  (define set-rest cdr)

  (define (set-singleton? s)
    (if (null? s) false (null? (cdr s))))

  (define (set-singleton x)
    (list x))

  (define (set-adjoin x set)
    (cond ((null? set) (list x))
	  ((set-equal-elements? x (car set)) set)
	  ((set-less-elements? x (car set)) (cons x set))
	  (else (cons (car set) (set-adjoin x (cdr set))))))

  (define (set-remove x set)
    (cond ((null? set) '())
	  ((set-equal-elements? x (car set)) (cdr set))
	  ((set-less-elements? x (car set)) set)
	  (else (cons (car set) (set-remove x (cdr set))))))

  (define (set-element? x set)
    (cond ((null? set) false)
	  ((set-equal-elements? x (car set)) true)
	  ((set-less-elements? x (car set)) false)
	  (else (set-element? x (cdr set)))))

  (define (set-intersection set1 set2)
    (cond ((null? set1) '())
	  ((null? set2) '())
	  ((set-equal-elements? (car set1) (car set2))
	   (cons (car set1) (set-intersection (cdr set1) (cdr set2))))
	  ((set-less-elements? (car set1) (car set2))
	   (set-intersection (cdr set1) set2))
	  (else (set-intersection set1 (cdr set2)))))

  (define (set-union set1 set2)
    (cond ((null? set1) set2)
	  ((null? set2) set1)
	  ((set-equal-elements? (car set1) (car set2))
	   (cons (car set1) (set-union (cdr set1) (cdr set2))))
	  ((set-less-elements? (car set1) (car set2))
	   (cons (car set1) (set-union (cdr set1) set2)))
	  (else (cons (car set2) (set-union set1 (cdr set2))))))

  (define (set-difference set1 set2)
    (cond ((null? set2) set1)
	  ((null? set1) '())
	  ((set-equal-elements? (car set1) (car set2))
	   (set-difference (cdr set1) (cdr set2)))
	  ((set-less-elements? (car set2) (car set1))
	   (set-difference set1 (cdr set2)))
	  (else (cons (car set1) (set-difference (cdr set1) set2)))))

  (define (set-subset? s1 s2)
    (cond ((null? s1) true)
	  ((null? s2) false)
	  ((set-equal-elements? (car s1) (car s2))
	   (set-subset? (cdr s1) (cdr s2)))
	  ((set-less-elements? (car s1) (car s2)) false)
	  (else (set-subset? s1 (cdr s2)))))

  (define (list->set lst)
    (define (remove-duplicates lst)
      (cond ((null? lst) lst)
	    ((null? (cdr lst)) lst)
	    ((set-equal-elements? (car lst) (cadr lst))
	     (remove-duplicates (cdr lst)))
	    (else
	     (cons (car lst)
		   (remove-duplicates (cdr lst))))))
    (remove-duplicates (sort lst set-less-elements?)))

  (define (set->list set) set)

  (vector the-empty-set
	  set-empty?
	  set-singleton
	  set-singleton?
	  set-adjoin
	  set-remove
	  set-element?
	  set-intersection
	  set-union
	  set-difference
	  set-subset?
	  list->set
	  set->list))

(define (empty-set set-type) (vector-ref set-type 0))
(define (empty-set? set-type) (vector-ref set-type 1))
(define (singleton-set set-type) (vector-ref set-type 2))
(define (singleton-set? set-type) (vector-ref set-type 3))
(define (adjoin-set set-type) (vector-ref set-type 4))
(define (remove-set set-type) (vector-ref set-type 5))
(define (element-set? set-type) (vector-ref set-type 6))
(define (intersect-sets set-type) (vector-ref set-type 7))
(define (union-sets set-type) (vector-ref set-type 8))
(define (difference-sets set-type) (vector-ref set-type 9))
(define (subset-sets? set-type) (vector-ref set-type 10))
(define (list->set set-type) (vector-ref set-type 11))
(define (set->list set-type) (vector-ref set-type 12))

(define symbols (make-sets-package eq? symbol<?))
(define numbers (make-sets-package = <))

#|
;;; For example

((set->list symbols)
 ((union-sets symbols)
  ((list->set symbols) '(a c e))
  ((list->set symbols) '(d e f))))
;Value: (a c d e f)
|#