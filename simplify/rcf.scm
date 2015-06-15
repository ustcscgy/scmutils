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

;;;; Rational Forms constructed over polynomials

(declare (usual-integrations))

;;; The following procedures are mapped for polynomial arithmetic.
;;; Other implementation of polynomial arithmetic may be constructed
;;; by modifying this map.

(define pcf:zero poly:zero)
(define pcf:one poly:one)
(define pcf:zero? poly:zero?)
(define pcf:one? poly:one?)
(define pcf:leading-base-coefficient poly:leading-base-coefficient)
(define pcf:normalize-by poly:normalize-by)

(define pcf:pcf? poly?)

(define pcf:= poly:=)
(define pcf:+ poly:+)
(define pcf:* poly:*)

#|
;;; Now set in pcf-fpf.scm to sparse stuff
(define (pcf:gcd x y)
  ((poly/heuristic-gcd poly/gcd-memoized) x y))
|#

(define pcf:quotient poly:quotient)
(define pcf:- poly:-)
(define pcf:negate poly:negate)
(define pcf:value poly:value)
(define pcf:derivative poly:derivative)
(define pcf:expt poly:expt)
(define pcf:new-variables poly:new-variables)
(define pcf:principal-reverse poly:principal-reverse)
(define pcf:arg-scale poly:arg-scale)
(define pcf:arg-shift poly:arg-shift)
(define pcf:degree poly:degree)
(define pcf:arity poly:arity)
(define pcf:extend poly:extend)

(define int:even? even?)

(define rcf-tag '*RCF*)

(define (ratform? object)
  (and (pair? object)
       (eq? (car object) rcf-tag)))

(define (make-ratform n d)
  (list rcf-tag n d))

(define ratform-numerator cadr)
(define ratform-denominator caddr)


(define rcf:zero pcf:zero)
(define rcf:one pcf:one)

(define (rcf:zero? q)
  (and (not (ratform? q)) (pcf:zero? q)))

(define (rcf:one? q)
  (and (not (ratform? q)) (pcf:one? q)))

(define (rcf:arity q)	 
  (define (check-same-arity p1 p2)
    (let ((a1 (pcf:arity p1)) (a2 (pcf:arity p2)))
      (cond ((fix:= a1 0) a2)
	    ((fix:= a2 0) a1)
	    ((fix:= a1 a2) a1)
	    (else
	     (error "Unequal arities in RCF" q)))))
  (cond ((ratform? q)
	 (check-same-arity (ratform-numerator q)
			   (ratform-denominator q)))
	((pcf:pcf? q) (pcf:arity q))
	(else (error "Wrong type -- RCF:ARITY" q))))

(define (make-rcf n d)
  (cond ((pcf:zero? d)
	 (error "Zero denominator -- MAKE-RCF" n d))
	((or (pcf:zero? n) (pcf:one? d)) n)
	((number? d) (pcf:* (/ 1 d) n))
	(else
	 (let ((dn (pcf:leading-base-coefficient d)))
	   (if (pcf:one? dn)
	       (make-ratform n d)
	       (make-ratform (pcf:normalize-by n dn)
			     (pcf:normalize-by d dn)))))))


(define (rcf:rcf? object)
  (or (ratform? object) (pcf:pcf? object)))

(define (rcf:pcf? object)
  (and (not (ratform? object)) (pcf:pcf? object)))

(define (rcf:= q r)
  (if (ratform? q)
      (if (ratform? r)
	  (and (pcf:= (ratform-numerator q) (ratform-numerator r))
	       (pcf:= (ratform-denominator q) (ratform-denominator r)))
	  (if (pcf:pcf? r)
	      #f
	      (error "Wrong type -- RCF:=" r)))
      (if (ratform? r)
	  (if (pcf:pcf? q)
	      #f
	      (error "Wrong type -- RCF:=" q))
	  (pcf:= q r))))


;;; The notation here is from Knuth (p. 291).
;;; In various places we take the gcd of two numbers and then call
;;; quotient to reduce those numbers.

(define (rcf:+ u/u* v/v*)
  (rcf:binary-operator u/u* v/v*
     pcf:+
     (lambda (u v v*)
       (if (pcf:zero? u)
	   v/v*
	   (make-rcf (pcf:+ (pcf:* u v*) v) v*)))
     (lambda (u u* v)
       (if (pcf:zero? v)
	   u/u*
	   (make-rcf (pcf:+ u (pcf:* u* v)) u*)))
     (lambda (u u* v v*)
       (if (pcf:= u* v*)
	   (let* ((n (pcf:+ u v)) (g (pcf:gcd u* n)))
	     (if (pcf:one? g)
		 (make-rcf n u*)
		 (make-rcf (pcf:quotient n g) (pcf:quotient u* g))))
	   (let ((d1 (pcf:gcd u* v*)))
	     (if (pcf:one? d1)
		 (make-rcf (pcf:+ (pcf:* u v*) (pcf:* u* v))
			   (pcf:* u* v*))
		 (let* ((u*/d1 (pcf:quotient u* d1))
			(t (pcf:+ (pcf:* u (pcf:quotient v* d1))
				  (pcf:* u*/d1 v))))
		   (if (pcf:zero? t)
		       rcf:zero
		       (let ((d2 (pcf:gcd t d1)))
			 (if (pcf:one? d2)
			     (make-rcf t (pcf:* u*/d1 v*))
			     (make-rcf
			      (pcf:quotient t d2)
			      (pcf:* u*/d1
				     (pcf:quotient v* d2)))))))))))))

(define (rcf:- u/u* v/v*)
  (rcf:binary-operator u/u* v/v*
     pcf:-
     (lambda (u v v*)
       (if (pcf:zero? u)
	   (make-ratform (pcf:negate v) v*)
	   (make-rcf (pcf:- (pcf:* u v*) v) v*)))
     (lambda (u u* v)
       (if (pcf:zero? v)
	   u/u*
	   (make-rcf (pcf:- u (pcf:* u* v)) u*)))
     (lambda (u u* v v*)
       (if (pcf:= u* v*)
	   (let* ((n (pcf:- u v)) (g (pcf:gcd u* n)))
	     (if (pcf:one? g)
		 (make-rcf n u*)
		 (make-rcf (pcf:quotient n g) (pcf:quotient u* g))))
	   (let ((d1 (pcf:gcd u* v*)))
	     (if (pcf:one? d1)
		 (make-rcf (pcf:- (pcf:* u v*) (pcf:* u* v))
			   (pcf:* u* v*))
		 (let* ((u*/d1 (pcf:quotient u* d1))
			(t (pcf:- (pcf:* u (pcf:quotient v* d1))
				  (pcf:* u*/d1 v))))
		   (if (pcf:zero? t)
		       rcf:zero
		       (let ((d2 (pcf:gcd t d1)))
			 (if (pcf:one? d2)
			     (make-rcf t (pcf:* u*/d1 v*))
			     (make-rcf
			      (pcf:quotient t d2)
			      (pcf:* u*/d1
				     (pcf:quotient v* d2)))))))))))))

(define (rcf:negate v/v*)
  (if (ratform? v/v*)
      (make-ratform (pcf:negate (ratform-numerator v/v*))
		    (ratform-denominator v/v*))
      (pcf:negate v/v*)))

(define (rcf:* u/u* v/v*)
  (rcf:binary-operator u/u* v/v*
    pcf:*
    (lambda (u v v*)
      (cond ((pcf:zero? u) rcf:zero)
	    ((pcf:one? u) v/v*)
	    (else
	     (let ((d (pcf:gcd u v*)))
	       (if (pcf:one? d)
		   (make-rcf (pcf:* u v) v*)
		   (make-rcf (pcf:* (pcf:quotient u d) v)
			     (pcf:quotient v* d)))))))
    (lambda (u u* v)
      (cond ((pcf:zero? v) rcf:zero)
	    ((pcf:one? v) u/u*)
	    (else
	     (let ((d (pcf:gcd u* v)))
	       (if (pcf:one? d)
		   (make-rcf (pcf:* u v) u*)
		   (make-rcf (pcf:* u (pcf:quotient v d))
			     (pcf:quotient u* d)))))))
    (lambda (u u* v v*)
      (let ((d1 (pcf:gcd u v*))
	    (d2 (pcf:gcd u* v)))
	(if (pcf:one? d1)
	    (if (pcf:one? d2)
		(make-rcf (pcf:* u v) (pcf:* u* v*))
		(make-rcf (pcf:* u (pcf:quotient v d2))
			  (pcf:* (pcf:quotient u* d2) v*)))
	    (if (pcf:one? d2)
		(make-rcf (pcf:* (pcf:quotient u d1) v)
			  (pcf:* u* (pcf:quotient v* d1)))
		(make-rcf (pcf:* (pcf:quotient u d1)
				 (pcf:quotient v d2))
			  (pcf:* (pcf:quotient u* d2)
				 (pcf:quotient v* d1)))))))))

(define (rcf:square q)
  (if (ratform? q)
      (make-ratform (let ((n (ratform-numerator q))) (pcf:* n n))
		    (let ((d (ratform-denominator q))) (pcf:* d d)))
      (pcf:* q q)))

(define (rcf:/ u/u* v/v*)
  (rcf:* u/u* (rcf:invert v/v*)))

(define (rcf:invert v/v*)
  (make-rcf (rcf:denominator v/v*)
	    (rcf:numerator v/v*)))

(define (rcf:gcd u/u* v/v*)
  (rcf:binary-operator u/u* v/v*
    pcf:gcd
    (lambda (u v v*)
      (cond ((pcf:zero? u) v/v*)
	    ((pcf:one? u) pcf:one)
	    (else (pcf:gcd u v))))
    (lambda (u u* v)
      (cond ((pcf:zero? v) u/u*)
	    ((pcf:one? v) pcf:one)
	    (else (pcf:gcd u v))))
    (lambda (u u* v v*)
      (let ((d1 (pcf:gcd u v))
	    (d2 (pcf:gcd u* v*)))
	(make-rcf d1 d2)))))

(define (rcf:binary-operator u/u* v/v* int*int int*rat rat*int rat*rat)
  (if (ratform? u/u*)
      (if (ratform? v/v*)
	  (rat*rat (ratform-numerator u/u*)
		   (ratform-denominator u/u*)
		   (ratform-numerator v/v*)
		   (ratform-denominator v/v*))
	  (rat*int (ratform-numerator u/u*)
		   (ratform-denominator u/u*)
		   v/v*))
      (if (ratform? v/v*)
	  (int*rat u/u*
		  (ratform-numerator v/v*)
		  (ratform-denominator v/v*))
	  (int*int u/u* v/v*))))

(define (rcf:numerator q)
  (cond ((ratform? q) (ratform-numerator q))
	((pcf:pcf? q) q)
	(else (error "Wrong type -- NUMERATOR" q))))

(define (rcf:denominator q)
  (cond ((ratform? q) (ratform-denominator q))
	((pcf:pcf? q) pcf:one)
	(else (error "Wrong type -- DENOMINATOR" q))))


(define (rcf:expt base exponent)
  (define (expt-iter x count answer)
    (if (fix:zero? count)
	answer
	(if (int:even? count)
	    (expt-iter (rcf:square x) (fix:quotient count 2) answer)
	    (expt-iter x (fix:-1+ count) (rcf:* x answer)))))
  (if (fix:negative? exponent)
      (rcf:invert (expt-iter base (int:negate exponent) rcf:one))
      (expt-iter base exponent rcf:one)))

(define (rcf:arg-scale r points)
  (if (ratform? r)
      (rcf:/ (apply pcf:arg-scale (ratform-numerator r) points)
	     (apply pcf:arg-scale (ratform-denominator r) points))
      (apply pcf:arg-scale r points)))

(define (rcf:arg-shift r points)
  (if (ratform? r)
      (rcf:/ (apply pcf:arg-shift (ratform-numerator r) points)
	     (apply pcf:arg-shift (ratform-denominator r) points))
      (apply pcf:arg-shift r points)))

(define (rcf:value r points)
  (if (ratform? r)
      (rcf:/ (apply pcf:value (ratform-numerator r) points)
	     (apply pcf:value (ratform-denominator r) points))
      (apply pcf:value r points)))

;;; The following only plugs r2 in for the principal indeterminate.

(define (rcf:compose r1 r2)
  (if (ratform? r2)
      (let ((nr1 (ratform-numerator r1))    (nr2 (ratform-numerator r2))
	    (dr1 (ratform-denominator r1))  (dr2 (ratform-denominator r2)))
	(let ((dn (pcf:degree nr1))
	      (dd (pcf:degree dr1))
	      (narity (fix:+ (pcf:arity dr1) 1)))
	  (let ((nnr1 (pcf:extend 1 (pcf:principal-reverse nr1)))
		(ndr1 (pcf:extend 1 (pcf:principal-reverse dr1))))
	    (let ((scales (list (cadr (pcf:new-variables narity)) 1)))
	      (let ((pn (pcf:value (pcf:principal-reverse
				      (pcf:arg-scale nnr1 scales))
				   nr2
				   dr2))
		    (pd (pcf:value (pcf:principal-reverse
				    (pcf:arg-scale ndr1 scales))
				   nr2
				   dr2)))
		(cond ((fix:> dn dd)
		       (rcf:/ pn (pcf:* (pcf:expt dr2 (fix:- dn dd)) pd)))
		      ((fix:< dn dd)
		       (rcf:/ (pcf:* (pcf:expt dr2 (fix:- dd dn)) pn) pd))
		      (else (rcf:/ pn pd))))))))
      (rcf:/ (pcf:value (ratform-numerator r1) r2)
	     (pcf:value (ratform-denominator r1) r2))))

(define (rcf:derivative r varnum)
  (if (ratform? r)
      (let ((u (ratform-numerator r)) (v (ratform-denominator r)))
	(rcf:/ (pcf:- (pcf:* (pcf:derivative u varnum) v)
		      (pcf:* u (pcf:derivative v varnum)))
	       (pcf:* v v)))
      (pcf:derivative r varnum)))

;;; I don't know if this stuff is ever important...GJS

(define (assoc-accumulation rat:op poly:op rat:identity)
  (define (operate rats)
    (cond ((null? rats) rat:identity)
	  ((null? (cdr rats)) (car rats))
	  ((ratform? (car rats))
	   (cond ((ratform? (cadr rats))
		  (operate (cons (rat:op (car rats) (cadr rats))
				 (cddr rats))))
		 ((null? (cddr rats)) (rat:op (car rats) (cadr rats)))
		 ((not (ratform? (caddr rats)))
		  (operate (cons (car rats)
				 (cons (poly:op (cadr rats) (caddr rats))
				       (cdddr rats)))))
		 (else
		  (operate (cons (rat:op (car rats) (cadr rats))
				 (cddr rats))))))
	  ((ratform? (cadr rats))
	   (operate (cons (rat:op (car rats) (cadr rats))
			  (cddr rats))))
	  (else
	   (operate (cons (poly:op (car rats) (cadr rats))
			  (cddr rats))))))
  (lambda rats (operate rats)))

(define +$rcf (assoc-accumulation rcf:+ pcf:+ rcf:zero))
(define *$rcf (assoc-accumulation rcf:* pcf:* rcf:one))


(define (assoc-inverse-accumulation rat:inv-op rat:op rat:invert poly:op rat:identity)
  (let ((direct-op (assoc-accumulation rat:op poly:op rat:identity)))
    (define (operate . rats)
      (cond ((null? rats) rat:identity)
	    ((null? (cdr rats)) (rat:invert (car rats)))
	    (else (rat:inv-op (car rats) (apply direct-op (cdr rats))))))
    operate))

(define -$rcf (assoc-inverse-accumulation rcf:- rcf:+ rcf:negate pcf:+ rcf:zero))
(define /$rcf (assoc-inverse-accumulation rcf:/ rcf:* rcf:invert pcf:* rcf:one))

;;; For simplifier

(define (rcf:->expression p vars)
  (if (poly? p)
      (poly:->expression p vars)
      (symb:/ (poly:->expression (ratform-numerator p) vars)
	      (poly:->expression (ratform-denominator p) vars))))

(define (rcf:expression-> expr cont #!optional less?)
  ;; cont = (lambda (ratp vars) ... )
  (let ((evars
	 (sort (list-difference (variables-in expr)
				rcf:operators-known)
		(if (default-object? less?) alphaless? less?))))
    (cont ((expression-walker
	    (pair-up evars
		     (pcf:new-variables (length evars))
		     rcf:operator-table))
	   expr)
	  evars)))

(define (rcf:->lambda p)
  (if (poly? p)
      (poly->function p)
      (let* ((n (rcf:arity p))
	     (vars (generate-list-of-symbols 'x n))
	     (exp (rcf:->expression p vars)))	  
	`(lambda ,vars ,exp))))


(define rcf:operator-table
  `((+        ,+$rcf)
    (-        ,-$rcf)
    (*        ,*$rcf)
    (/        ,/$rcf)
    (negate   ,rcf:negate)
    (invert   ,rcf:invert)
    (expt     ,rcf:expt)
    (square   ,rcf:square)
    (gcd      ,rcf:gcd)))

(define rcf:operators-known
  (map car rcf:operator-table))
