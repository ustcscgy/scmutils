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

(define (compile-and-run-sexp sexp #!optional environment declarations syntax-table)
  (scode-eval (if (default-object? declarations)
		  (compile-sexp sexp)
		  (compile-sexp sexp declarations))
	      (if (default-object? environment)
		  scmutils-base-environment
		  environment)))

(define (compile-sexp sexp #!optional declarations syntax-table)
  (let ((scode
	 (syntax&integrate (text/cselim sexp)
			   (if (default-object? declarations)
			       '((USUAL-INTEGRATIONS))
			       declarations)
			   (if (default-object? syntax-table)
			       system-global-syntax-table
			       syntax-table))))
    (if (or *do-not-compile*
	    (lexical-unreferenceable? system-global-environment 'COMPILE-PROCEDURE))
	scode
	(let* ((compiler-package (->environment '(compiler top-level)))
	       (compile-scode (access compile-scode compiler-package)))
	  (fluid-let (((access compiler:show-procedures? compiler-package)
		       false))
	    (compile-scode scode))))))


(define (compile-and-run-numerical sexp #!optional environment declarations syntax-table)
  (scode-eval (if (default-object? declarations)
		  (compile-sexp sexp)
		  (compile-numerical sexp declarations))       
	      (if (default-object? environment)
		  scmutils-base-environment
		  environment)))

(define (compile-numerical sexp #!optional declarations syntax-table)
  (let ((scode
	 (syntax&integrate (text/cselim (flonumize sexp))
			   (if (default-object? declarations)
			       *numerical-integrations*
			       declarations)
			   (if (default-object? syntax-table)
			       system-global-syntax-table
			       syntax-table))))
    (if (or *do-not-compile*
	    (lexical-unreferenceable? system-global-environment 'COMPILE-PROCEDURE))
	scode
	(let* ((compiler-package (->environment '(compiler top-level)))
	       (compile-scode (access compile-scode compiler-package)))
	  (fluid-let (((access compiler:show-procedures? compiler-package)
		       false))
	    (compile-scode scode))))))

(define (compile-generic->numerical generic-procedure #!optional simplify)
  (define (walk object)
    (cond ((pair? object)
	   (cons (walk (car object)) (walk (cdr object))))
	  ((primitive-procedure? object)
	   (primitive-procedure-name object))
	  (else object)))
  (let ((sexp (walk (unsyntax (procedure-lambda generic-procedure)))))
    (scode-eval (compile-numerical
		 (if (default-object? simplify)
		     sexp
		     (simplify sexp)))
		numerical-environment)))

(define *numerical-integrations*
  '((integrate-external "/usr/local/lib/mit-scheme/floint")))

(define *do-not-compile* false)

(define (flonumize exp)
  (cond ((list? exp)
	 (if (symbol? (car exp))
	     (case (car exp)
	       ((declare access quote) exp)
	       ((define define-integrable)
		`(,(car exp) ,(cadr exp) ,@(flonumize (cddr exp))))
	       ((lambda named-lambda)
		`(,(car exp) ,(cadr exp) ,@(map flonumize (cddr exp))))
	       ((let let* fluid-let)
		(if (symbol? (cadr exp))
		    `(,(car exp) ,(cadr exp)
				 ,(map (lambda (b)
					 (cons (car b)
					       (if (null? (cdr b))
						   '()
						   (list (flonumize (cadr b))))))
				       (caddr exp))
				 ,@(map flonumize (cdddr exp)))    
		    `(,(car exp)
		      ,(map (lambda (b)
			      (cons (car b)
				    (if (null? (cdr b))
					'()
					(list (flonumize (cadr b))))))
			    (cadr exp))
		      ,@(map flonumize (cddr exp)))))
	       ((vector-ref list-ref matrix-ref) exp)
	       ((make-vector)
		`(make-vector ,(cadr exp) ,(flonumize (caddr exp))))
	       ((vector-set!)
		`(vector-set! ,(cadr exp)
			      ,(caddr exp)
			      ,(flonumize (cadddr exp))))
	       ((matrix-set!)
		`(matrix-set! ,(cadr exp)
			      ,(caddr exp)
			      ,(cadddr exp)
			      ,(flonumize (list-ref exp 4))))
	       ((set!)
		`(set! ,(cadr exp) ,(flonumize (caddr exp))))
	       ((expt)
		(let ((base (flonumize (cadr exp)))
		      (e (caddr exp)))
		  (if (and (integer? e) (exact? e))
		      (case e
			((2) `(* ,base ,base))
			((3) `(* ,base ,base ,base))
			((4) `(* ,base ,base ,base ,base))
			(else `(expt ,base ,e))))))
	       (else
		(if (string-prefix? "fix:" (symbol->string (car exp)))
		    exp
		    (cons (car exp)
			  (map flonumize (cdr exp))))))
	     (map flonumize exp)))
	((number? exp) (exact->inexact exp))
	((vector? exp) ((vector-elementwise flonumize) exp))
	(else exp)))

#|
(define *numerical-integrations*
  '((usual-integrations + - * / < > = zero? positive? negative?
			sqrt abs exp sin cos #| atan |#)
    (reduce-operator
     (+ flo:+ (null-value 0. none) (group left))
     (- flo:- (null-value 0. single) (group left))
     (* flo:* (null-value 1. none) (group left))
     (/ flo:/ (null-value 1. single) (group left)))
    (integrate-primitive-procedures
     (< flonum-less?)
     (> flonum-greater?)
     (= flonum-equal?)
     (zero? flonum-zero?)
     (positive? flonum-positive?)
     (negative? flonum-negative?)
     (sqrt flonum-sqrt)
     (abs flonum-abs)
     (exp flonum-exp)
     (sin flonum-sin)
     (cos flonum-cos)
     #|(atan flonum-atan 1)		;
     (atan flonum-atan2 2)|#)))
|#

#|
Date: Thu, 30 Aug 90 07:06:55 edt
From: "Guillermo J. Rozas" <jinx@altdorf.ai.mit.edu>
Return-Path: <jinx@altdorf.ai.mit.edu>
To: hal@altdorf.ai.mit.edu, gjs@altdorf.ai.mit.edu, wisdom@altdorf.ai.mit.edu
Subject: Floating arithmetic
Reply-To: jinx@zurich.ai.mit.edu

To change the default arithmetic to floating point in a file,
use

(declare (usual-integrations + - * /)
	 (reduce-operator
	  (+ flo:+ (null-value 0. none) (group left))
	  (- flo:- (null-value 0. single) (group left))
	  (* flo:* (null-value 1. none) (group left))
	  (/ flo:/ (null-value 1. single) (group left))))

instead of 

(declare (usual-integrations))

You may need to add similar constructs for =, >, etc., and other
declarations for the irrational functions.

|#
