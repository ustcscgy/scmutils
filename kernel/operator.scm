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

;;;; Operators

(declare (usual-integrations))


(define (o:type o) operator-type-tag)

(define (o:type-predicate o) operator?)

(define (o:arity o) 
  (operator-arity o))


#|;;; In GENERIC.SCM

(define (make-operator p #!optional name subtype arity #!rest opts)
  (if (default-object? name) (set! name #f))
  (if (default-object? subtype) (set! subtype #f))
  (if (default-object? arity) (set! arity (procedure-arity p)))
  (make-apply-hook p `(,operator-type-tag ,subtype ,name ,arity ,@opts)))
|#


(define (operator-procedure op)
  (assert (operator? op))
  (apply-hook-procedure op))

(define (operator-subtype op)
  (assert (operator? op))
  (cadr (apply-hook-extra op)))

(define (operator-name op)
  (assert (operator? op))
  (caddr (apply-hook-extra op)))

(define (operator-arity op)
  (assert (operator? op))
  (cadddr (apply-hook-extra op)))

(define (operator-optionals op)
  (assert (operator? op))
  (cddddr (apply-hook-extra op)))


(define (operator-merge-subtypes op1 op2)
  (let ((t1 (operator-subtype op1))
	(t2 (operator-subtype op2)))
    (cond ((eq? t1 t2) t1)
	  ((not t1)    t2)
	  ((not t2)    t1)
	  (else
	   (error "Incompatible subtypes -- OPERATOR" t1 t2)))))

(define (operator-merge-arities op1 op2)
  (joint-arity (operator-arity op1) (operator-arity op2)))


(define (o:zero-like op)
  (assert (equal? (operator-arity op) *exactly-one*) "o:zero-like")
  (make-operator
   (lambda (f) (g:zero-like f))
   'zero
   (operator-subtype op)))

(define (o:one-like op)
  (assert (equal? (operator-arity op) *exactly-one*) "o:one-like")
  (make-operator g:identity 'identity (operator-subtype op)))

(define o:identity
  (make-operator g:identity 'identity))

(define (o:+ op1 op2)
  (make-operator
   (lambda fs
     (g:+ (apply op1 fs) (apply op2 fs)))
   `(+ ,(operator-name op1)
       ,(operator-name op2))
   (operator-merge-subtypes op1 op2)
   (operator-merge-arities op1 op2)))

(define (o:o+f op f)
  (make-operator
   (lambda (g)
     (g:+ (op g) (g:* f g)))
   `(+ ,(operator-name op) ,f)
   (operator-subtype op)))

(define (o:f+o f op)
  (make-operator
   (lambda (g)
     (g:+ (g:* f g) (op g)))
   `(+ ,f ,(operator-name op))
   (operator-subtype op)))


(define (o:- op1 op2)
  (make-operator
   (lambda fs
     (g:- (apply op1 fs) (apply op2 fs)))
   `(- ,(operator-name op1)
       ,(operator-name op2))
   (operator-merge-subtypes op1 op2)
   (operator-merge-arities op1 op2)))

(define (o:o-f op f)
  (make-operator
   (lambda (g)
     (g:- (op g) (g:* f g)))
   `(- ,(operator-name op) ,f)
   (operator-subtype op)))

(define (o:f-o f op)
  (make-operator
   (lambda (g)
     (g:- (g:* f g) (op g)))
   `(- ,f ,(operator-name op))
   (operator-subtype op)))

(define (o:negate op)
  (make-operator
   (lambda fs
     (g:negate (apply op fs)))
   `(- ,(operator-name op))
   (operator-subtype op)
   (operator-arity op)))

(define (o:* op1 op2)
  (let ((subtype
	 (operator-merge-subtypes op1 op2)))
    (if (procedure? subtype)
	(subtype op1 op2)
	(make-operator
	 (compose op1 op2)
	 `(* ,(operator-name op1)
	     ,(operator-name op2))
	 subtype
	 (operator-arity op2)))))

(define (o:f*o f op)
  (make-operator
   (lambda gs
     (g:* f (apply op gs)))
   `(* ,f ,(operator-name op))
   (operator-subtype op)
   (operator-arity op)))

(define (o:o*f op f)
  (make-operator
   (lambda gs
     (apply op (map (lambda (g) (g:* f g)) gs)))
   `(* ,(operator-name op) ,f)
   (operator-subtype op)
   (operator-arity op)))

(define (o:expt op n)
  (assert (equal? (operator-arity op) *exactly-one*) "o:expt")
  (make-operator
   (iterated op n o:identity)
   `(expt ,(operator-name op) ,n)
   (operator-subtype op)))

(define (o:exp op)
  (assert (equal? (operator-arity op) *exactly-one*) "o:exp")
  (make-operator
   (lambda (g)
     (lambda x
       (g:apply ((series:value exp-series (list op)) g) x)))
   `(exp ,(operator-name op))
   (operator-subtype op)))

(define (o:cos op)
  (assert (equal? (operator-arity op) *exactly-one*) "o:cos")
  (make-operator
   (lambda (g)
     (lambda x
       (g:apply ((series:value cos-series (list op)) g) x)))
   `(cos ,(operator-name op))
   (operator-subtype op)))

(define (o:sin op)
  (assert (equal? (operator-arity op) *exactly-one*) "o:sin")
  (make-operator
   (lambda (g)
     (lambda x
       (g:apply ((series:value sin-series (list op)) g) x)))
   `(sin ,(operator-name op))
   (operator-subtype op)))


;;; Optional order argument for exponentiation of operators.
;;; (((expn D 2) g) x)
;;;   = (((exp D)
;;;       (lambda (eps)
;;;        (((+ 1 (* (expt eps 2) D) (* 1/2 (expt eps 4) (expt D 2)) ...) g) x))
;;;      0)
;;; This is (exp (* (expt eps 2) D)) written as a power series in eps.

(define (expn op #!optional exponent)
  (assert (operator? op))
  (assert (equal? (operator-arity op) *exactly-one*) "o:expn")
  (if (default-object? exponent)
      (o:exp op)
      (make-operator
       (lambda (g)
	 (lambda x
	   (g:apply ((series:inflate (series:value exp-series (list op))
				      exponent)
		     g)
		    x)))
       `(exp ,(operator-name op))
       (operator-subtype op))))


(assign-operation 'type                o:type            operator?)
(assign-operation 'type-predicate      o:type-predicate  operator?)
(assign-operation 'arity               o:arity           operator?)

(assign-operation 'zero-like o:zero-like operator?)
(assign-operation 'one-like o:one-like operator?)
(assign-operation 'identity-like o:one-like operator?)

(assign-operation '+          o:+               operator? operator?) 
(assign-operation '+          o:o+f             operator? not-operator?) 
(assign-operation '+          o:f+o             not-operator? operator?) 

(assign-operation '-          o:-               operator? operator?)
(assign-operation '-          o:o-f             operator? not-operator?) 
(assign-operation '-          o:f-o             not-operator? operator?) 

(assign-operation '*          o:*               operator? operator?)
(assign-operation '*          o:o*f             operator? not-operator?) 
(assign-operation '*          o:f*o             not-operator? operator?)

(assign-operation 'negate     o:negate          operator?)
(assign-operation 'expt       o:expt            operator? exact-integer?)

(assign-operation 'exp                 o:exp             operator?)
(assign-operation 'sin                 o:sin             operator?)
(assign-operation 'cos                 o:cos             operator?)

