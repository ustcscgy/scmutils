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

(declare (usual-integrations))

;;; Generates a list of all permutations of a list of distinct elements.

(define (permutations lst)
  (if (null? lst)
      '(())
      (apply append
	     (map (lambda (item)
		    (map (lambda (perm) (cons item perm))
			 (permutations (delete item lst))))
		  lst))))


;;; Returns the number of interchanges required to generate the
;;;  permuted list from the original list.

(define (list-interchanges permuted-list original-list)
  (let lp1 ((plist permuted-list) (count 0))
    (if (null? plist)
	count
	(let ((bigger (cdr (member (car plist) original-list)))
	      (rest (cdr plist)))
	  (let lp2 ((l rest) (increment 0))
	    (if (null? l)
		(lp1 rest
		     (fix:+ count increment))
		(lp2 (cdr l)
		     (if (not (member (car l) bigger))
			 (fix:+ increment 1)
			 increment))))))))

(define (split-permutations original-list list-of-permutations cont)
  ;; cont = (lambda (even-permutations odd-permutations) ... )
  (let lp ((perms list-of-permutations) (evens '()) (odds '()))
    (if (null? perms)
	(cont evens odds)
	(let ((sig (list-interchanges (car perms) original-list)))
	  (if (even? sig)
	      (lp (cdr perms) (cons (car perms) evens) odds)
	      (lp (cdr perms) evens (cons (car perms) odds)))))))


;;; Returns the number of interchanges required to generate the
;;;  permuted list of numbers from an ordered list.

(define (permutation-interchanges permuted-list)
  (let lp1 ((plist permuted-list) (count 0))
    (if (null? plist)
	count
	(let ((first (car plist))
	      (rest (cdr plist)))
	  (let lp2 ((l rest) (increment 0))
	    (if (null? l)
		(lp1 rest
		     (fix:+ count increment))
		(lp2 (cdr l)
		     (if (int:> (car l) first)
			 increment
			 (fix:+ increment 1)))))))))


;;; Given a permutation (represented as a list of numbers),
;;;  and a list to be permuted, construct the list so permuted.

(define (permute permutation lst)
  (map (lambda (p)
	 (list-ref lst p))
       permutation))


;;; Sometimes we want to permute some of the elements of a list, as follows:
;;; (subpermute '((1 . 4) (4 . 2) (2 . 3) (3 . 1)) '(a b c d e))
;;; ;Value 6: (a e d b c)

(define (subpermute the-map lst)
  (let* ((n (length lst)))
    (let lp ((i 0) (source lst) (answer '()))
      (if (fix:= i n)
	  (reverse answer)
	  (let ((entry (assoc i the-map)))
	    (if (not entry)
		(lp (fix:+ i 1)
		    (cdr source)
		    (cons (car source) answer))
		(lp (fix:+ i 1)
		    (cdr source)
		    (cons (list-ref lst (cdr entry)) answer))))))))


(define (factorial n)
  (if (int:< n 2)
      1
      (int:* n (factorial (int:- n 1)))))

(define number-of-permutations factorial)

(define (number-of-combinations  n k)
  (int:quotient (factorial n)
		(int:* (factorial (int:- n k))
		       (factorial k))))
