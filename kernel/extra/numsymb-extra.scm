#| -*-Scheme-*-

Copyright (C) 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994,
    1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005,
    2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013 Massachusetts
    Institute of Technology

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


(define (s:+ a1 a2)
  (cond ((sum? a1)
	 (cond ((sum? a2)
		(addup-args (append (operands a1) (operands a2)) '()))
	       ((difference? a2)
		(if (null? (cdr (operands a2)))
		    (addup-args (operands a1) (operands a2))
		    (addup-args (append (operands a1) (list (car (operands a2))))
				(cdr (operands a2)))))
	       (else (addup-args (append (operands a1) (list a2)) '()))))
	((difference? a1)
	 (if (null? (cdr (operands a1)))
	     (cond ((sum? a2) (addup-args (operands a2) (operands a1)))
		   ((difference? a2)
		    (if (null? (cdr (operands a2)))
			(addup-args '() (append (operands a1) (operands a2)))
			(addup-args (list (car (operands a2)))
				    (append (operands a1) (cdr (operands a2))))))
		   (else (addup-args (list a2) (operands a1))))
	     (cond ((sum? a2)
		    (addup-args (append (list (car (operands a1))) (operands a2))
				(cdr (operands a1))))
		   ((difference? a2)
		    (if (null? (cdr (operands a2)))
			(addup-args (list (car (operands a1)))
				    (append (cdr (operands a1)) (operands a2)))
			(addup-args (list (car (operands a1)) (car (operands a2)))
				    (append (cdr (operands a1)) (cdr (operands a2))))))
		   (else (addup-args (list (car (operands a1)) a2)
				     (cdr (operands a1)))))))
	(else
	 (cond ((sum? a2)
		(addup-args (append (list a1) (operands a2)) '()))
	       ((difference? a2)
		(if (null? (cdr (operands a2)))
		    (addup-args (list a1) (operands a2))
		    (addup-args (append (list a1) (list (car (operands a2))))
				(cdr (operands a2)))))
	       (else (addup-args (list a1 a2) '()))))))

	   
(define (addup-args pos neg)
  (define (make-answer sum pos neg)
    (if (zero? sum)
	(if (null? pos)
	    (if (null? neg)
		:zero
		(if (null? (cdr neg))
		    `(- ,(car neg))
		    `(- (+ ,@neg))))
	    (if (null? neg)
		(if (null? (cdr pos))
		    (car pos)
		    `(+ ,@pos))
		(if (null? (cdr pos))
		    (if (null? (cdr neg))
			`(- ,(car pos) ,(car neg))
			`(- ,(car pos) (+ ,@neg)))
		    (if (null? (cdr neg))
			`(- (+ ,@pos) ,(car neg))
			`(- (+ ,@pos) (+ ,@neg))))))
	(if (null? pos)
	    (if (null? neg)
		sum
		(if (null? (cdr neg))
		    `(- ,sum ,(car neg))
		    `(- ,sum (+ ,@neg))))
	    (if (null? neg)
		`(+ ,sum ,@pos)
		(if (null? (cdr neg))
		    `(- (+ ,sum ,@pos) ,(car neg))
		    `(- (+ ,sum ,@pos) (+ ,@neg)))))))
  (let plp ((p pos) (sum :zero) (respos '()))
    (cond ((null? p)
	   (let nlp ((n neg) (sum sum) (resneg '()))
	     (cond ((null? n)
		    (make-answer sum
				 (reverse respos)
				 (reverse resneg)))
		   ((number? (car n))
		    (nlp (cdr n) (- sum (car n)) resneg))
		   (else
		    (nlp (cdr n) sum (cons (car n) resneg))))))
	  ((number? (car p))
	   (plp (cdr p) (+ sum (car p)) respos))
	  (else
	   (plp (cdr p) sum (cons (car p) respos))))))


(define (s:* a1 a2)
  (cond ((product? a1)
	 (cond ((product? a2)
		(mulup-args (append (operands a1) (operands a2)) '()))
	       ((quotient? a2)
		(if (null? (cdr (operands a2)))
		    (mulup-args (operands a1) (operands a2))
		    (mulup-args (append (operands a1) (list (car (operands a2))))
				(cdr (operands a2)))))
	       (else (mulup-args (append (operands a1) (list a2)) '()))))
	((quotient? a1)
	 (if (null? (cdr (operands a1)))
	     (cond ((product? a2) (mulup-args (operands a2) (operands a1)))
		   ((quotient? a2)
		    (if (null? (cdr (operands a2)))
			(mulup-args '() (append (operands a1) (operands a2)))
			(mulup-args (list (car (operands a2)))
				    (append (operands a1) (cdr (operands a2))))))
		   (else (mulup-args (list a2) (operands a1))))
	     (cond ((product? a2)
		    (mulup-args (append (list (car (operands a1))) (operands a2))
				(cdr (operands a1))))
		   ((quotient? a2)
		    (if (null? (cdr (operands a2)))
			(mulup-args (list (car (operands a1)))
				    (append (cdr (operands a1)) (operands a2)))
			(mulup-args (list (car (operands a1)) (car (operands a2)))
				    (append (cdr (operands a1)) (cdr (operands a2))))))
		   (else (mulup-args (list (car (operands a1)) a2)
				     (cdr (operands a1)))))))
	(else
	 (cond ((product? a2)
		(mulup-args (append (list a1) (operands a2)) '()))
	       ((quotient? a2)
		(if (null? (cdr (operands a2)))
		    (mulup-args (list a1) (operands a2))
		    (mulup-args (append (list a1) (list (car (operands a2))))
				(cdr (operands a2)))))
	       (else (mulup-args (list a1 a2) '()))))))


(define (mulup-args pos neg)
  (define (make-answer factor pos neg)
    (if (zero? factor)
	factor
	(if (one? factor)
	    (if (null? pos)
		(if (null? neg)
		    :one
		    (if (null? (cdr neg))
			`(/ ,:one ,(car neg))
			`(/ ,:one (* ,@neg))))
		(if (null? neg)
		    (if (null? (cdr pos))
			(car pos)
			`(* ,@pos))
		    (if (null? (cdr pos))
			(if (null? (cdr neg))
			    `(/ ,(car pos) ,(car neg))
			    `(/ ,(car pos) (* ,@neg)))
			(if (null? (cdr neg))
			    `(/ (* ,@pos) ,(car neg))
			    `(/ (* ,@pos) (* ,@neg))))))
	    (if (null? pos)
		(if (null? neg)
		    factor
		    (if (null? (cdr neg))
			`(/ ,factor ,(car neg))
			`(/ ,factor (* ,@neg))))
		(if (null? neg)
		    `(* ,factor ,@pos)
		    (if (null? (cdr neg))			
			`(/ (* ,factor ,@pos) ,(car neg))
			`(/ (* ,factor ,@pos) (* ,@neg))))))))
  (let plp ((p pos) (factor :one) (respos '()))
    (cond ((null? p)
	   (let nlp ((n neg) (factor factor) (resneg '()))
	     (cond ((null? n)
		    (make-answer factor
				 (reverse respos)
				 (reverse resneg)))
		   ((number? (car n))
		    (nlp (cdr n) (/ factor (car n)) resneg))
		   (else
		    (nlp (cdr n) factor (cons (car n) resneg))))))
	  ((number? (car p))
	   (plp (cdr p) (* factor (car p)) respos))
	  (else
	   (plp (cdr p) factor (cons (car p) respos))))))

	   
(define (s:- a1 a2)
  (cond ((sum? a1)
	 (cond ((sum? a2)
		(addup-args (operands a1) (operands a2)))
	       ((difference? a2)
		(if (null? (cdr (operands a2)))
		    (addup-args (append (operands a1) (operands a2)) '())
		    (addup-args (append (operands a1) (cdr (operands a2)))
				(list (car (operands a2))))))
	       (else (addup-args (operands a1) (list a2)))))
	((difference? a1)
	 (if (null? (cdr (operands a1)))
	     (cond ((sum? a2) (addup-args '() (append (operands a1) (operands a2))))
		   ((difference? a2)
		    (if (null? (cdr (operands a2)))
			(addup-args (operands a2) (operands a1))
			(addup-args (cdr (operands a2))
				    (append (operands a1)
					    (list (car (operands a2)))))))
		   (else (addup-args '() (append (operands a1) (list a2)))))
	     (cond ((sum? a2)
		    (addup-args (list (car (operands a1)))
				(append (cdr (operands a1)) (operands a2))))
		   ((difference? a2)
		    (if (null? (cdr (operands a2)))
			(addup-args (append (list (car (operands a1))) (operands a2))
				    (cdr (operands a1)))
			(addup-args (cons (car (operands a1)) (cdr (operands a2)))
				    (append (cdr (operands a1))
					    (list (car (operands a2)))))))
		   (else (addup-args (list (car (operands a1)))
				     (append (cdr (operands a1)) (list a2)))))))
	(else
	 (cond ((sum? a2)
		(addup-args (list a1) (operands a2)))
	       ((difference? a2)
		(if (null? (cdr (operands a2)))
		    (addup-args (append (list a1) (operands a2)) '())
		    (addup-args (append (list a1) (cdr (operands a2)))
				(list (car (operands a2))))))
	       (else (addup-args (list a1) (list a2)))))))


(define (s:/ a1 a2)
  (cond ((product? a1)
	 (cond ((product? a2)
		(mulup-args (operands a1) (operands a2)))
	       ((quotient? a2)
		(if (null? (cdr (operands a2)))
		    (mulup-args (append (operands a1) (operands a2)) '())
		    (mulup-args (append (operands a1) (cdr (operands a2)))
				(list (car (operands a2))))))
	       (else (mulup-args (operands a1) (list a2)))))
	((quotient? a1)
	 (if (null? (cdr (operands a1)))
	     (cond ((product? a2) (mulup-args '() (append (operands a1) (operands a2))))
		   ((quotient? a2)
		    (if (null? (cdr (operands a2)))
			(mulup-args (operands a2) (operands a1))
			(mulup-args (cdr (operands a2))
				    (append (operands a1)
					    (list (car (operands a2)))))))
		   (else (mulup-args '() (append (operands a1) (list a2)))))
	     (cond ((product? a2)
		    (mulup-args (list (car (operands a1)))
				(append (cdr (operands a1)) (operands a2))))
		   ((quotient? a2)
		    (if (null? (cdr (operands a2)))
			(mulup-args (append (list (car (operands a1))) (operands a2))
				    (cdr (operands a1)))
			(mulup-args (cons (car (operands a1)) (cdr (operands a2)))
				    (append (cdr (operands a1))
					    (list (car (operands a2)))))))
		   (else (mulup-args (list (car (operands a1)))
				     (append (cdr (operands a1)) (list a2)))))))
	(else
	 (cond ((product? a2)
		(mulup-args (list a1) (operands a2)))
	       ((quotient? a2)
		(if (null? (cdr (operands a2)))
		    (mulup-args (append (list a1) (operands a2)) '())
		    (mulup-args (append (list a1) (cdr (operands a2)))
				(list (car (operands a2))))))
	       (else (mulup-args (list a1) (list a2)))))))
