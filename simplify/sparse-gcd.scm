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

;;;;             Sparse Multivariate Polynomial GCD 
;;;    a probabilistic method inspired by Richard Zippel's thesis
;;;      coded and debugged by Gerald Jay Sussman and Dan Zuras  
;;;                         June 1998                  

;;; This code differs from Zippel's in that it does not use modular
;;; arithmetic or do anything special with the Vandermonde matrices
;;; that arise in the problem.  This makes the idea of using sparse
;;; interpolation stand out in stark contrast without the confusing
;;; complications introduced by those optimizations.

(declare (usual-integrations))

;;; The purpose of sparse-gcd is to remove and replace the content, so
;;; as to present the sparse-multivariate-gcd program with primitive
;;; polynomials.  The sparse-content is defined to include the largest
;;; monomial factor, as well, thus lowering the degree as much as
;;; possible.

;;; This code is only called from pcf-fpf.scm, which ensures that the
;;; arguments to sparse-gcd are of the same arity.

(define (sparse-gcd u v)
  (cond ((null? u) v)
	((null? v) u)
	((equal? u v) u)
	((sparse-one? u) u)
	((sparse-one? v) v)
	(else
	 (let ((uc (sparse-content u))
	       (vc (sparse-content v)))
	   (let ((ans 
		  (if (sparse-one-term? uc)
		      (if (sparse-one-term? vc)
			  (sparse-multivariate-gcd u v)
			  (sparse-multivariate-gcd u (sparse-normalize v vc)))
		      (if (sparse-one-term? vc)
			  (sparse-multivariate-gcd (sparse-normalize u uc) v)
			  (let ((c (sparse-monomial-gcd uc vc)))
			    (if (sparse-one-term? c)
				(sparse-multivariate-gcd (sparse-normalize u uc)
							 (sparse-normalize v vc))
				(sparse-scale
				 (sparse-multivariate-gcd (sparse-normalize u uc)
							  (sparse-normalize v vc))
				 c)))))))
	     (sparse-abs ans))))))


;;; sparse-multivariate-gcd determines the maximum possible degree of
;;; the gcd in each variable -- the minimum of the max in each input
;;; polynomial.

(define (sparse-multivariate-gcd P Q)
  (let ((n (length (sparse-exponents (car P))))
	(dPs (map (lambda (l) (apply max l))
		  (list-transpose (map car P))))
	(dQs (map (lambda (l) (apply max l))
		  (list-transpose (map car Q)))))
    (assert (= n (length (sparse-exponents (car Q)))))
    (let ((ds (map min dPs dQs)))
      (sparse-multivariate-gcd-helper P Q n ds))))

(define (sparse-monomial-gcd m1 m2)
  (sparse-term (map min (sparse-exponents m1) (sparse-exponents m2))
	       (base/gcd (sparse-coefficient m1) (sparse-coefficient m2))))

(define (sparse-content poly)
  (let lp ((p (cdr poly)) (ans (car poly)))
    (cond ((null? p) ans)
	  ((sparse-one-term? ans) ans)
	  (else
	   (lp (cdr p)
	       (sparse-monomial-gcd (car p) ans))))))

(define (sparse-multivariate-gcd-helper P Q n ds)
  (reset-interpolation-args! ds
			     (apply max (map sparse-coefficient P))
			     (apply max (map sparse-coefficient Q)))
  (let restart ((time0 (runtime)))
    (if *sgcd-wallp* (pp 'restart))
    (let* ((rargs1 (make-interpolation-args (- n 1)))
	   (P1 (sparse-evaluate> P rargs1))
	   (Q1 (sparse-evaluate> Q rargs1))
	   (g1 (univariate-gcd P1 Q1)))
      (let stagelp ((k 1) (g g1) (rargs rargs1))	
	;; g has k vars interpolated to make it arity k.
	(if (= k n)
	    g
	    (let* ((skeleton (map sparse-exponents g))
		   (nterms (length skeleton))
		   (trial-arglists
		    (begin
		      (if *sgcd-wallp*
			  (pp `(sparse-gcd:
				(k ,k)
				(nterms ,nterms)
				(skeleton ,skeleton)
				(time ,(- (runtime) time0)))))
		      (generate-list nterms
			 (lambda (i)
			   (make-interpolation-args k)))))
		   (Pk (sparse-evaluate> P (cdr rargs)))
		   (Qk (sparse-evaluate> Q (cdr rargs)))
		   (Gks (map (lambda (arglist)
			       (univariate-gcd
				(sparse-evaluate< Pk arglist)
				(sparse-evaluate< Qk arglist)))
			     trial-arglists))
		   (GkSkels (map (lambda (Gk)
				   (map sparse-exponents Gk))
				 Gks)))

	      (if (not (all-equal? GkSkels))
		  (begin (if *sgcd-wallp*
			     (pp `(sparse-gcd: GkSkels-not-same ,GkSkels)))
			 (stagelp k g rargs))
		  (let ((xk+1s
			 (generate-list (+ (list-ref ds k) 1)
					interpolate-random)))
		    (lu-decompose
		     (matrix-by-row-list
		      (map (lambda (arguments)
			     (map (lambda (exponents)
				    (apply *
					   (map expt
						arguments
						exponents)))
				  skeleton))
			   trial-arglists))
		     (lambda (lu-matrix lu-permutation lu-sign)
		       (let ((coeffs
			      (map (lambda (xk+1)
				     (let ((values
					    (map (lambda (Gk)
						   (sparse-evaluate
						    Gk
						    (list xk+1)))
						 Gks)))
				       (vector->list
					(lu-backsubstitute
					 lu-matrix
					 lu-permutation
					 (list->vector values)))))
				   xk+1s)))
			 (let clp ((css (list-transpose coeffs)) (cps '()))
			   (if (null? css)
			       (let ((cps (reverse cps)))
				 (let ((new-g (expand-poly g cps)))
				   (if (and (sparse-divisible? Pk new-g)
					    (sparse-divisible? Qk new-g))
				       (stagelp (fix:+ k 1) new-g (cdr rargs))
				       (begin (if *sgcd-wallp*
						  (pp `(sparse-gcd: division)))
					      (restart time0)))))
			       (univariate-interpolate-values
				xk+1s (car css)
				(lambda (cp) (clp (cdr css) (cons cp cps)))
				(lambda ()
				  (if *sgcd-wallp*
				      (pp `(sparse-gcd: interpolation)))
				  (restart time0)))))))
		     (lambda (x)
		       (if *sgcd-wallp* (pp `(sparse-gcd: singular)))
		       (restart time0)))))))))))

(define *sgcd-wallp* #f)

(define *interpolate-primes-stream* '())

(define (reset-interpolation-args! ds max-c-p max-c-q)
  (set! *interpolate-primes-stream* prime-numbers-stream)
  'done)

(define (make-interpolation-args k)
  (let lp ((i 0) (s *interpolate-primes-stream*) (args '()))
    (if (= i k)
	(begin (set! *interpolate-primes-stream* s)
	       args)
	(lp (+ i 1) (tail s) (cons (head s) args)))))

#|
;;; This is trying to be a good boy, using the formula from Zippel for
;;; the mod prime, but I think it is not really necessary.  Timings 
;;; at the end are using this choice...

(define (reset-interpolation-args! ds max-c-p max-c-q)
  (first-prime-stream-exceeding!
   (max max-c-p max-c-p (apply max ds)))
  'done)

(define (first-prime-stream-exceeding! n)
  (let lp ((s prime-numbers-stream))
    (if (> (head s) n)
	(set! *interpolate-primes-stream* s)
	(lp (tail s)))))
|#
#|
;;; This randomly works pretty well, but... I don't understand why.

(define (reset-interpolation-args! ds max-c-p max-c-q)
  'done)

(define (make-interpolation-args k)
  (generate-list k interpolate-prime))

(define (interpolate-prime i)
  (stream-ref prime-numbers-stream (random 5000)))
|#
#|
;;; I tried relatively-prime stuff, and it doesn't really do the job!
(define *interpolation-args* '())

(define (reset-interpolation-args! ds max-c-p max-c-q)
  (set! *interpolation-args* '()))

(define (make-interpolation-args k)
  (let next ((i 0) (args '()))
    (if (= i k)
	args
	(let try-again ((trial (random *interpolate-size*)))
	  (if (for-all? *interpolation-args*
		(lambda (a) (= (base/gcd a trial) 1)))
	      (begin (set! *interpolation-args*
			   (cons trial *interpolation-args*))
		     (next (+ i 1) (cons trial args)))
	      (try-again (random *interpolate-size*)))))))
|#

(define (univariate-gcd u v)		;Euclid's Algorithm is OK here.
  (define (pgcd ppu ppv)
    (if *ugcd-wallp* (pp `((ppu: ,ppu) (ppv: ,ppv))))
    (cond ((null? ppv) ppu)		;v=0      => u
	  ((sparse-constant? ppv)	;deg(v)=0 => 1
	   univariate-one)
	  (else
	   (pgcd ppv
		 (univariate-primitive-part
		  (univariate-pseudo-remainder ppu ppv))))))
  (cond ((null? u) v)
	((null? v) u)
	((sparse-constant? u)
	 (univariate-constant (base/gcd (sparse-coefficient (car u)) (univariate-content v))))
	((sparse-constant? v)
	 (univariate-constant (base/gcd (univariate-content u) (sparse-coefficient (car v)))))
	(else
	 (let ((uc (univariate-content u))
	       (vc (univariate-content v)))
	   (let ((ans
		  (if (= uc 1)
		      (if (= vc 1)
			  (pgcd u v)
			  (pgcd u (univariate-normalize v vc)))
		      (if (= vc 1)
			  (pgcd (univariate-normalize u uc) v)
			  (let ((c (base/gcd uc vc)))
			    (if (= c 1)
				(pgcd (univariate-normalize u uc)
				      (univariate-normalize v vc))
				(univariate-scale
				 (pgcd (univariate-normalize u uc)
				       (univariate-normalize v vc))
				 c)))))))
	     (sparse-abs ans))))))

(define *ugcd-wallp* #f)

(define (univariate-content poly)
  (let lp ((p (cdr poly)) (ans (sparse-coefficient (car poly))))
    (cond ((null? p) ans)
	  ((= ans 1) 1)
	  (else
	   (lp (cdr p)
	       (base/gcd (sparse-coefficient (car p))
			 ans))))))

(define (univariate-primitive-part poly)
  (if (null? poly)
      '()
      (univariate-normalize poly (univariate-content poly))))

(define (univariate-pseudo-remainder u v)
  (let ((cvn (sparse-coefficient (car v)))            ;leading coefficient of v
	(n (car (sparse-exponents (car v)))))         ;degree v
    (let lp ((u u) )
      (if (null? u)
	  '()
	  (let ((cum (sparse-coefficient (car u)))    ;leading coefficient of u
		(m (car (sparse-exponents (car u))))) ;degree u
	    (if (< m n)
		u
		(lp (sparse-add
		     (univariate-scale u cvn)
		     (sparse-multiply-term
		      (sparse-term (list (- m n)) (- cum))
		      v)))))))))

(define (univariate-constant coeff)
  (list (sparse-term '(0) coeff)))

(define univariate-one
  (univariate-constant 1))

(define (univariate-scale p c)
  (map (lambda (term)
	 (sparse-term (sparse-exponents term)
		      (* c (sparse-coefficient term))))
       p))

(define (univariate-normalize p c)
  (map (lambda (term)
	 (sparse-term (sparse-exponents term)
		      (/ (sparse-coefficient term) c)))
       p))

#|
;;; Knuth's test
(univariate-gcd '(((8) . 1) ((6) . 1) ((4) . -3) ((3) . -3) ((2) . 8) ((1) . 2) ((0) . -5))
		'(((6) . 3) ((4) . 5) ((2) . -4) ((1) . -9) ((0) . 21)))

((ppu: (((8) . 1) ((6) . 1) ((4) . -3) ((3) . -3) ((2) . 8) ((1) . 2) ((0) . -5)))
 (ppv: (((6) . 3) ((4) . 5) ((2) . -4) ((1) . -9) ((0) . 21))))
((ppu: (((6) . 3) ((4) . 5) ((2) . -4) ((1) . -9) ((0) . 21)))
 (ppv: (((4) . -5) ((2) . 1) ((0) . -3))))
((ppu: (((4) . -5) ((2) . 1) ((0) . -3)))
 (ppv: (((2) . -13) ((1) . -25) ((0) . 49))))
((ppu: (((2) . -13) ((1) . -25) ((0) . 49)))
 (ppv: (((1) . -4663) ((0) . 6150))))
((ppu: (((1) . -4663) ((0) . 6150)))
 (ppv: (((0) . 1))))
;Value: (((0) . 1))
|#

#|
(define (gcd-test d f g)
  (let ((pd (fpf:expression-> d (lambda (p v) p)))
	(pf (fpf:expression-> f (lambda (p v) p)))
	(pg (fpf:expression-> g (lambda (p v) p))))
    (let ((pdf (fpf:* pd pf)) (pdg (fpf:* pd pg)))
      (let ((gcd?
	     (sparse-gcd
	      (fpf:->sparse pdf)
	      (fpf:->sparse pdg)))
	    (ans
	     (fpf:->sparse pd)))
	(if (equal? gcd? ans)
	    #t
	    (pp (list (list 'gcd? gcd?)
		      (list 'ans ans))))))))


(define d1
  '(+ (expt x1 2) x1 3))

(define f1
  '(+ (* 2 (expt x1 2)) (* 2 x1) 1))

(define g1
  '(+ (expt x1 2) (* 2 x1) 2))

;(show-time (lambda () (gcd-test d1 f1 g1)))
;process time: 10 (10 RUN + 0 GC); real time: 5
;Value: #t


(define d2
  '(+ (* 2 (expt x1 2) (expt x2 2))
      (* x1 x2)
      (* 2 x1)))

(define f2
  '(+ (expt x2 2)
      (* 2 (expt x1 2) x2)
      (expt x1 2)
      1))

(define g2
  '(+ (* (expt x1 2) (expt x2 2))
      (* (expt x1 2) x2)
      (* x1 x2)
      (expt x1 2)
      x1))

;(show-time (lambda () (gcd-test d2 f2 g2)))
;process time: 40 (40 RUN + 0 GC); real time: 34
;Value: #t

(define d3
  '(+ (* x2 x2 x3 x3)
      (* x2 x2 x3)
      (* 2 x1 x1 x2 x3)
      (* x1 x3)))

(define f3
  '(+ (* x3 x3)
      (* x2 x2 x3)
      (* x1 x1 x2 x3)
      (* x1 x3)
      (* x1 x1 x2 x2)))

(define g3
  '(+ (* x2 x3)
      (* 2 x1 x3)
      x3
      x1))

;(show-time (lambda () (gcd-test d3 f3 g3)))
;process time: 80 (80 RUN + 0 GC); real time: 83
;Value: #t


(define d4
  '(+ (* x1 x1 x4 x4)
      (* x2 x2 x3 x4)
      (* x1 x1 x2 x4)
      (* x2 x4)
      (* x1 x1 x2 x3)))

(define f4
  '(+ (* x1 x2 x3 x3 x4 x4)
      (* x1 x3 x3 x4 x4)
      (* x1 x4 x4)
      (* x4 x4)
      (* x1 x3 x4)))

(define g4
  '(+ (* x1 x3 x3 x4 x4)
      (* x3 x3 x4 x4)
      (* x4 x4)
      (* x1 x2 x2 x3 x4)
      (* x1 x2 x2)))

;(show-time (lambda () (gcd-test d4 f4 g4)))
;process time: 240 (240 RUN + 0 GC); real time: 238
;Value: #t

(define d5
  '(+ (* x1 x1 x1 x2 x2 x3 x3 x4 x5 x5)
      (* x1 x2 x2 x5 x5)
      (* x1 x1 x1 x3 x4 x4 x5)
      (* x1 x1 x1 x2 x3 x3 x4 x5)
      (* x1 x1 x2 x3 x3 x4 x4)))

(define f5
  '(+ (* x1 x2 x2 x5 x5)
      (* x1 x2 x3 x3 x4 x5)
      (* x1 x2 x3 x3 x4 x4)
      (* x1 x2 x2 x4 x4)
      1))

(define g5
  '(+ (* x1 x3 x3 x4 x5 x5)
      (* x2 x5 x5)
      (* x1 x2 x4 x5)
      (* x2 x5)
      (* x1 x2 x3 x4 x4)))

;(show-time (lambda () (gcd-test d5 f5 g5)))
;process time: 450 (450 RUN + 0 GC); real time: 453
;Value: #t


(define d6
  '(+ (* x1 x2 x4 x4 x5 x5 x6 x6)
      (* x1 x2 x2 x3 x3 x4 x5 x5 x6 x6)
      (* x1 x1 x3 x6 x6)
      (* x1 x1 x2 x3 x3 x4 x5 x5 x6)
      (* x1 x1 x3 x5 x6)))

(define f6
  '(+ (* x1 x1 x2 x4 x5 x5 x6 x6)
      (* x1 x3 x5 x5 x6 x6)
      (* x1 x2 x2 x6 x6)
      (* x1 x1 x2 x2 x3 x3 x5 x6)
      (* x1 x3 x3 x4 x5)))

(define g6
  '(+ (* x2 x2 x3 x3 x4 x5 x5 x6)
      (* x1 x4 x4 x5 x6)
      (* x2 x2 x3 x3 x4 x5 x6)
      (* x1 x2 x2 x3 x4 x4 x6)
      (* x1 x1 x3 x5 x5)))

;(show-time (lambda () (gcd-test d6 f6 g6)))
;process time: 460 (460 RUN + 0 GC); real time: 460
;Value: #t

(define d7
  '(+ (* x1 x2 x2 x4 x4 x6 x6 x7 x7)
      (* x1 x1 x3 x4 x6 x6 x7 x7)
      (* x3 x3 x4 x4 x7 x7)
      (* x1 x1 x2 x4 x4 x6)
      (* x3 x4 x5 x5)))

(define f7
  '(+ (* x1 x1 x2 x4 x4 x5 x6 x6 x7 x7)
      (* x1 x2 x3 x6 x7)
      (* x3 x4 x4 x5 x5 x7)
      (* x1 x1 x2 x3 x4 x4 x5 x6)))

(define g7
  '(+ (* x1 x3 x5 x6 x6 x7 x7)
      (* x2 x2 x3 x3 x4 x4 x5 x6 x7 x7)
      (* x4 x6 x7 x7)
      (* x1 x1 x2 x3 x5 x6 x7)
      (* x1 x1 x3 x3 x4 x5 x5)))


;(show-time (lambda () (gcd-test d7 f7 g7)))
;process time: 690 (690 RUN + 0 GC); real time: 685
;Value: #t


(define d8
  '(+ (* x2 x2 x4 x5 x6 x7 x8 x8)
      (* x1 x1 x2 x3 x3 x4 x4 x6 x6 x7 x7 x8)
      (* x1 x1 x3 x4 x4 x6 x6 x7 x7)
      (* x1 x1 x2 x2 x3 x3 x4 x5 x5 x6 x7 x7)
      (* x2 x2 x4 x6)))

(define f8
  '(+ (* x1 x1 x2 x2 x3 x4 x4 x5 x6 x6 x8 x8)
      (* x2 x5 x6 x6 x8 x8)
      (* x1 x1 x2 x2 x3 x3 x4 x4 x6 x6 x7 x7 x8)
      (* x1 x1 x3 x3 x4 x5 x5 x7 x7 x8)
      (* x1 x2 x2 x3 x3 x5 x5 x7)))

(define g8
  '(+ (* x1 x4 x4 x6 x6 x7 x8 x8)
      (* x1 x2 x2 x4 x4 x5 x5 x6 x6 x8)
      (* x1 x1 x2 x3 x4 x4 x6 x6 x8)
      (* x1 x1 x2 x2 x3 x3 x4 x5 x5 x8)
      (* x1 x2 x4 x4 x5 x5)))

;(show-time (lambda () (gcd-test d8 f8 g8)))
;process time: 900 (900 RUN + 0 GC); real time: 894
;Value: #t

(define d10
  '(+ (* x1 x2 x2 x4 x4 x8 x9 x9 x10 x10)
      (* x2 x2 x4 x5 x5 x6 x7 x9 x10 x10)
      (* x1 x1 x2 x3 x5 x5 x7 x7 x9 x9)
      (* x1 x3 x3 x4 x4 x7 x7 x9 x9)
      (* x1 x1 x3 x4 x7 x7 x8 x8)))

(define f10
  '(+ (* x1 x2 x3 x3 x4 x6 x7 x8 x9 x9 x10 x10)
      (* x2 x2 x3 x3 x4 x4 x6 x6 x9 x10 x10)
      (* x1 x2 x2 x3 x3 x4 x5 x6 x7 x8 x8 x9 x9 x10)
      (* x1 x1 x2 x4 x4 x5 x5 x8 x8 x9 x9 x10)
      (* x3 x4 x4 x5 x6 x7 x7 x9 x10)))

(define g10
  '(+ (* x1 x2 x2 x3 x3 x5 x5 x6 x6 x7 x8 x9 x9 x10 x10)
      (* x3 x8 x9 x9 x10 x10)
      (* x1 x2 x2 x3 x4 x5 x5 x6 x6 x8 x8 x9 x10)
      (* x1 x3 x6 x7 x8 x10)
      (* x4 x4 x5 x5 x6 x6 x7 x9 x9)))
 
;(show-time (lambda () (gcd-test d10 f10 g10)))
;process time: 1550 (1550 RUN + 0 GC); real time: 1553
;Value: #t

;;; These are a bit harder versions of the problem

(define d10a
  '(+ (* 2 x1 x2 x2 x4 x4 x8 x9 x9 x10 x10)
      (* 3 x2 x2 x4 x5 x5 x6 x7 x9 x10 x10)
      (* 4 x1 x1 x2 x3 x5 x5 x7 x7 x9 x9)
      (* 5 x1 x3 x3 x4 x4 x7 x7 x9 x9)
      (* 6 x1 x1 x3 x4 x7 x7 x8 x8)
      7))

(define f10a
  '(+ (* 8 x1 x2 x3 x3 x4 x6 x7 x8 x9 x9 x10 x10)
      (* 9 x2 x2 x3 x3 x4 x4 x6 x6 x9 x10 x10)
      (* 10 x1 x2 x2 x3 x3 x4 x5 x6 x7 x8 x8 x9 x9 x10)
      (* 11 x1 x1 x2 x4 x4 x5 x5 x8 x8 x9 x9 x10)
      (* 12 x3 x4 x4 x5 x6 x7 x7 x9 x10)
      13))

(define g10a
  '(+ (* 14 x1 x2 x2 x3 x3 x5 x5 x6 x6 x7 x8 x9 x9 x10 x10)
      (* 15 x3 x8 x9 x9 x10 x10)
      (* 16 x1 x2 x2 x3 x4 x5 x5 x6 x6 x8 x8 x9 x10)
      (* 17 x1 x3 x6 x7 x8 x10)
      (* 18 x4 x4 x5 x5 x6 x6 x7 x9 x9)
      19))

;(show-time (lambda () (gcd-test d10a f10a g10a)))
process time: 2540 (2540 RUN + 0 GC); real time: 2534
;Value: #t
|#
