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

;;;; Property Tables

(declare (usual-integrations))

;;; Properties are n-dimensional sparse tables implemented as 
;;; nests of association lists.

;;; For any given sequence of keys, there can be both a value
;;; and a subtable.  A table is a list of a value and some entries.
;;; An entry is a pair, whose CAR is a key and whose CDR is a
;;; the subtable for that key.

(define (make-table table-name assoc)
  (let ((local-table (list *no-value*)))

    (define (lookup keys)
      (define (loop keys table)
	(if (null? keys) (car table)
	    (let ((entry (assoc (car keys) (cdr table))))
	      (if entry
		  (loop (cdr keys) (cdr entry))
		  *no-value*))))
      (loop keys local-table))

    (define (smash! keys value)
      (define (loop keys table)
	(if (null? keys) (set-car! table value)
	    (let ((entry (assoc (car keys) (cdr table))))
	      (if entry
		  (loop (cdr keys) (cdr entry))
		  (set-cdr! table
			    (cons (cons (car keys)
					(make-subtable (cdr keys) value))
				  (cdr table)))))))
      (loop keys local-table)
      local-table)

    (define (make-subtable keys value)
      (if (null? keys) (list value)
	  (list *no-value*
		(cons (car keys)
		      (make-subtable (cdr keys) value)))))

    (define (accumulator! increment-procedure initial-value keys value)
      (define (loop keys table)
	(if (null? keys)
	    (if (eq? (car table) *no-value*)
		(set-car! table (increment-procedure value initial-value))
		(set-car! table (increment-procedure value (car table))))
	    (let ((entry (assoc (car keys) (cdr table))))
	      (if entry
		  (loop (cdr keys) (cdr entry))
		  (set-cdr! table
			    (cons (cons (car keys)
					(make-subtable (cdr keys)
						       (increment-procedure value
									    initial-value)))
				  (cdr table)))))))
      (loop keys local-table)
      local-table)

    (define (remove! keys) (smash! keys *no-value*))

    (vector table-name lookup smash! accumulator! remove!)))


(define *no-value* (list '*no-value*))

(define (no-value? value)
  (eq? value *no-value*))


(define (get table . keys)
  ((vector-ref table 1) keys))

(define ((getter table) . keys)
  ((vector-ref table 1) keys))


(define (put! table value . keys)
  ((vector-ref table 2) keys value)
  'done)

(define ((putter! table) value . keys)
  ((vector-ref table 2) keys value)
  'done)


(define (get-with-default table default . keys)
  (let ((v ((vector-ref table 1) keys)))
    (if (eq? v *no-value*)
	default
	v)))

(define ((getter-with-default table default) . keys)
  (let ((v ((vector-ref table 1) keys)))
    (if (eq? v *no-value*)
	default
	v)))


(define (get-with-check table . keys)
  (let ((v ((vector-ref table 1) keys)))
    (if (eq? v *no-value*)
	(error "can't find value in table"
	       (list table keys))
	v)))

(define ((getter-with-check table) . keys)
  (let ((v ((vector-ref table 1) keys)))
    (if (eq? v *no-value*)
	(error "can't find value in table"
	       (list table keys))
	v)))


(define (add-to-list! object table . keys)
  ((vector-ref table 3) cons '() keys object)
  'done)

(define (adjoin-to-list! object table . keys)
  ((vector-ref table 3) list-adjoin '() keys object)
  'done)

(define (store! object table . keys)
  ((vector-ref table 2) keys object)
  'done)
