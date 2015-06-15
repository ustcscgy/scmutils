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

;;; variation of a path-dependent function

#|
(define ((((delta eta) f) q) t)
  ((D (lambda (h) ((f (+ q (* h eta))) t))) 0))
|#

(define (delta eta)
  (define (((delta-eta f) q) t)
    ((D (lambda (h) ((f (+ q (* h eta))) t))) 0))
  (make-operator delta-eta))

(print-expression
 ((((delta (literal-function 'eta))
    (lambda (path) path))
   (literal-function 'q))
  't))
(eta t)


(print-expression
 ((((delta (literal-function 'eta))
    (lambda (path) (D path)))
   (literal-function 'q))
  't))
((D eta) t)

(print-expression
 ((((delta (literal-function 'eta))
    (lambda (x) (- (* 1/2 'm (square (D x))) (* 1/2 'k (square x)))))
   (literal-function 'q))
  't))
(+ (* -1 k (q t) (eta t))
 (* m ((D q) t) ((D eta) t)))

#|
(print-expression
 ((((delta (literal-function 'eta))
    (compose (L-harmonic 'm 'k) Gamma))
   (literal-function 'q))
  't))
|#

(print-expression
 ((Gamma-bar
   (lambda (q)
     ((Gamma-bar 
       (lambda (eta)
	 ((((delta eta)
	    (lambda (x)
	      (- (* 1/2 'm (square (D x))) (* 1/2 'k (square x)))))
	   q)
	  't)))
     (->state 't 'eta 'Deta))))
  (->state 't 'q 'v)))
;;; fails now

(+ (* Deta eta m)
   (* -1 deta k q))

(define ((f q) t) ((literal-function 'f*) (q t)))
(define ((g q) t) ((literal-function 'g*) (q t)))


(pe (let ((q (literal-function 'q))
	  (eta (literal-function 'eta)))
      ((((delta eta) (+ f g)) q) 't)))
(+ (* ((D f*) (q t)) (eta t)) (* ((D g*) (q t)) (eta t)))

(pe (let ((q (literal-function 'q))
	  (eta (literal-function 'eta)))
      (- ((((delta eta) (+ f g)) q) 't)
	 (+ ((((delta eta) f) q) 't)
	    ((((delta eta) g) q) 't)))))
0


(pe (let ((q (literal-function 'q))
	  (eta (literal-function 'eta)))
      ((((delta eta) (* f g)) q) 't)))
(+ (* (g* (q t)) ((D f*) (q t)) (eta t))
   (* (f* (q t)) ((D g*) (q t)) (eta t)))

(pe (let ((q (literal-function 'q))
	  (eta (literal-function 'eta)))
      (- ((((delta eta) (* f g)) q) 't)
	 (+ (* ((((delta eta) f) q) 't) ((g q) 't))
	    (* ((f q) 't) ((((delta eta) g) q) 't))))))
0


(pe (let ((q (literal-function 'q))
	  (eta (literal-function 'eta)))
      ((((delta eta) (* 'c f)) q) 't)))
(* c ((D f*) (q t)) (eta t))


(pe (let ((q (literal-function 'q))
	  (eta (literal-function 'eta)))
      (- ((((delta eta) (* 'c f)) q) 't)
	 (* 'c ((((delta eta) f) q) 't)))))
0


(pe (let ((F (literal-function 'F))
	  (q (literal-function 'q))
	  (eta (literal-function 'eta)))
      ((((delta eta) 
	 (lambda (q) (compose F (f q))))
	 q) 't)))
(* ((D F) (f* (q t)))
   ((D f*) (q t)) (eta t))

(pe
 (let ((F (literal-function 'F))
       (q (literal-function 'q))
       (eta (literal-function 'eta)))
   (- ((((delta eta) 
	 (lambda (q) (compose F (f q))))
	q)
       't)
      ((* (compose (D F) (f q))
	  (((delta eta) f) q))
       't))))
0
	      

#|
(define ((f q) t)
  ((literal-function 'f*
		     (-> (X Real (^ Real 'integer) (^ Real 'integer))
			 Real))
   ((Gamma q) t)))

(define ((g q) t)
  ((literal-function 'g*
		     (-> (X Real (^ Real 'integer) (^ Real 'integer))
			 Real))
   ((Gamma q) t)))

;;; This should have worked, but the 3-argument thing should have been
;;; seen as 1 up-tuple argument.
|#

(define ((f q) t)
  ((literal-function 'f* (-> *vector* Real)) ((Gamma q) t)))

(define ((g q) t)
  ((literal-function 'g* (-> *vector* Real)) ((Gamma q) t)))


(pe (let ((F (literal-function 'F))
	  (q (literal-function 'q))
	  (eta (literal-function 'eta)))
      ((((delta eta) 
	 (lambda (q) (compose F (f q))))
	 q)
       't)))

(pe
 (let ((F (literal-function 'F))
       (q (literal-function 'q))
       (eta (literal-function 'eta)))
   (- ((((delta eta) 
	 (lambda (q) (compose F (f q))))
	q)
       't)
      ((* (compose (D F) (f q))
	  (((delta eta) f) q))
       't))))

#| D and delta commute


(define ((f q) t) ((literal-function 'f*) (q t)))
;Value: f

(define (delta eta)
  (define (((delta-eta f) q) t)
    ((D (lambda (h) ((f (+ q (* h eta))) t))) 0))
  (make-operator delta-eta))
;Value: delta

(pe ((((delta (literal-function 'eta)) f) (literal-function 'q)) 't))
(* ((D f*) (q t)) (eta t))
;Unspecified return value

(pe ((((delta (literal-function 'eta)) 
       (lambda (q) (lambda (t) ((D (f q)) t))))
      (literal-function 'q)) 't))
(+ (* (((expt D 2) f*) (q t)) (eta t) ((D q) t))
   (* ((D f*) (q t)) ((D eta) t)))
;Unspecified return value

(pe ((D (((delta (literal-function 'eta)) f)
      (literal-function 'q))) 't))
(+ (* (((expt D 2) f*) (q t)) ((D q) t) (eta t))
   (* ((D f*) (q t)) ((D eta) t)))
;Unspecified return value


|#

#|  look what we can do now!!


(define L (literal-function 'L (Lagrangian)))

(print-expression
 ((((delta (literal-function 'eta)) 
    (lambda (q) (compose L (Gamma q))))
   (literal-function 'q))
  't))
(+ (* (((partial 2) L) (up t (q t) ((D q) t))) ((D eta) t))
   (* (((partial 1) L) (up t (q t) ((D q) t))) (eta t)))
;Unspecified return value

(define L (literal-function 'L (Lagrangian 2)))

(print-expression
 ((((delta (literal-function 'eta)) 
    (lambda (q) (compose L (Gamma q))))
   (up (literal-function 'x) (literal-function 'y)))
  't))
;Bad structure -- MAKE-PARTIALS (1 0) (*diff* (() ...) (... ...))

;;; A red herring...

(print-expression
 ((((delta
     (up (literal-function 'eta)
	 (literal-function 'nu)))
    (lambda (q) (compose L (Gamma q))))
   (up (literal-function 'x) (literal-function 'y)))
  't))
(+
 (* (((partial 2 1) L) (up t (up (x t) (y t)) (up ((D x) t) ((D y) t))))
    ((D nu) t))
 (* (((partial 2 0) L) (up t (up (x t) (y t)) (up ((D x) t) ((D y) t))))
    ((D eta) t))
 (* (((partial 1 1) L) (up t (up (x t) (y t)) (up ((D x) t) ((D y) t))))
    (nu t))
 (* (((partial 1 0) L) (up t (up (x t) (y t)) (up ((D x) t) ((D y) t))))
    (eta t)))

|#

#|  can now do this ...

(define (F q)
  (compose 
   (literal-function 'f (-> (UP Real (UP* Real) (UP* Real)) Real))
   (Gamma q)))

(define (G q)
  (compose 
   (literal-function 'g (-> (UP Real (UP* Real) (UP* Real)) Real))
   (Gamma q)))

(print-expression
 (- ((((delta (literal-function 'eta)) (* F G))
      (literal-function 'q))
     't)
    (((+ (* ((delta (literal-function 'eta)) F) G)
	 (* F ((delta (literal-function 'eta)) G)))
      (literal-function 'q))
    't)))
0

|#
