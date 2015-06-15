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

;;;; Memoizers

(declare (usual-integrations))

;;; A general linear-time memoizer for functions.

;;;(define *auditing-memoizers* #f)
(define *auditing-memoizers* #t)

(define (linear-memoize f #!optional max-table-size ass)
  (let ((max-table-size		 ;set to -1 for no limit
	 (if (default-object? max-table-size)
	     12
	     max-table-size))
	(ass (if (default-object? ass) assoc ass))
	(table '())
	(size 0)
	(memo-hits 0)
	(memo-misses 0))
    (define (info)
      (list memo-hits memo-misses table))
    (define (reset)
      (set! memo-hits 0)
      (set! memo-misses 0)
      (set! size 0)
      (set! table '()))
    (set! *memoizers*
	  (cons (list f max-table-size info reset)
		*memoizers*))
    (let ((sm1 (fix:- max-table-size 1)))
      (define (memo-f . x)
	(let ((seen (ass x table)))
	  (if seen
	      (begin (if *auditing-memoizers*
			 (set! memo-hits (int:+ memo-hits 1)))
		     (cadr seen))
	      (let ((ans (apply f x)))
		(if *auditing-memoizers*
		    (set! memo-misses (int:+ memo-misses 1)))
		(cond ((fix:= size max-table-size)
		       (set! table
			     (cons (list x ans)
				   (list-head table sm1))))
		      (else ;CPH: Interrupt hole here!
		       (set! table
			     (cons (list x ans) table))
		       (set! size (fix:+ size 1))))
		ans))))
      memo-f)))

;;; If auditing, information is stored here:

(define *memoizers* '())

(define (show-memoizer-statistics)
  (pp (map (lambda (memoizer)
	   (cons (car memoizer) (list-head ((caddr memoizer)) 2)))
	 *memoizers*)))

(define (clear-memoizer-tables)
  (for-each (lambda (m) ((cadddr m))) *memoizers*))


;;; One appropriate lookup checks the arguments for EQV?uality. 

(define (ass-same-args args table)
  (let lp ((entries table))
    (cond ((null? entries) #f)
	  ((eqv-args? args (caar entries)) (car entries))
	  (else (lp (cdr entries))))))

(define (eqv-args? args1 args2)
  (cond ((null? args1)
	 (cond ((null? args2) #t)
	       (else #f)))
	((null? args2) #f)
	((eqv? (car args1) (car args2))
	 (eqv-args? (cdr args1) (cdr args2)))
	(else #f)))

(define *not-seen* (list '*not-seen*))

(define (hash-memoize f)
  (let ((table) (memo-hits) (memo-misses))
    (define (info)
      (list memo-hits memo-misses table))
    (define (reset)
      (set! memo-hits 0)
      (set! memo-misses 0)
      (set! table
	    ((weak-hash-table/constructor equal-hash-mod equal? #t))))
    (reset)
    (set! *memoizers*
	  (cons (list f -1 info reset)
		*memoizers*))
    (define (memo-f . x)
	(let ((seen (hash-table/get table x *not-seen*)))
	  (if (not (eq? seen *not-seen*))
	      (begin (if *auditing-memoizers*
			 (set! memo-hits (int:+ memo-hits 1)))
		     seen)
	      (let ((ans (apply f x)))
		(if *auditing-memoizers*
		    (set! memo-misses (int:+ memo-misses 1)))
		(hash-table/put! table x ans)
		ans))))
    memo-f))


