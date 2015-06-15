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

;;;;         Match and Substitution Language Interpreter

;;; Needs rule-memoize, rule-environment

(declare (usual-integrations))

;;;   This is the infamous 6.001 rule interpreter, originally written
;;; by GJS for a lecture in the faculty course held at MIT in the
;;; summer of 1983, and subsequently used and tweaked from time to
;;; time.  This subsystem has been a serious pain in the ass, because
;;; of its expressive limitations, but I have not had the guts to
;;; seriously improve it since its first appearance. -- GJS

(define (rule-simplifier the-rules)
  (define simplify-exp
    (rule-memoize
     (lambda (exp)
       (let ((result
	      (try-rules (if (list? exp)
			     (map simplify-exp exp)
			     exp)
			 the-rules)))
	 (cond ((equal? result exp) exp)
	       (else (simplify-exp result)))))))
  simplify-exp)

(define (try-rules exp the-rules)
  (define (scan rules)
    (if (null? rules)
	exp
	(match (rule-pattern (car rules))
	       exp
	       the-empty-dictionary
	       (lambda ()
		 (scan (cdr rules)))
	       (lambda (dictionary fail)
		 (if (check-predicate (rule-predicate (car rules))
				      dictionary)
		     (execute (rule-skeleton (car rules))
			      dictionary)
		     (fail))))))
  (scan the-rules))

(define (match pat dat dict fail succeed)
  (cond ((eq? pat dat)
	 (succeed dict fail))
	((number? pat)
	 (if (number? dat)
	     (if (= pat dat)
		 (succeed dict fail)
		 (fail))
	     (fail)))
	((arbitrary-element? pat)
	 (element-match pat dat dict fail succeed))
	((and (pair? pat) (arbitrary-segment? (car pat)))
	 (segment-match pat dat dict fail succeed))
	((and (pair? pat) (pair? dat))
	 (match (car pat) (car dat) dict fail
		(lambda (dict fail)
		  (match (cdr pat) (cdr dat) dict fail succeed))))
	(else (fail))))


(define (element-match pat dat dict fail succeed)
  (let ((vname (match-var-name pat))
	(p (match-var-restriction pat)))
    (if vname
	(if (symbol? vname)
	    (let ((v (match-lookup vname dict)))
	      (if v
		  (if (element-var? v)
		      (if (and (equal? (element-in v) dat)
			       (if p (p dat) #t))
			  (succeed dict fail)
			  (fail))
		      (error "Not an element variable" pat))
		  (if (if p (p dat) #t)
		      (succeed (extend-dict-element vname dat dict)
			       fail)
		      (fail))))
	    (match-call vname dat dict fail succeed))
	(if (if p (p dat) #t)
	    (succeed dict fail)
	    (fail)))))

(define (segment-match pat dat dict fail succeed)
  (let ((vname (match-var-name (car pat)))
	(p (match-var-restriction (car pat))))
    (if vname
	(if (symbol? vname)
	    (let ((v (match-lookup vname dict)))
	      (if v
		  (if (segment-var? v)
		      (let ((end (segment-end v)))
			(if (and (if p (p dat end) #t)
				 (let scan ((vptr (segment-beg v))
					    (dptr dat))
				   (cond ((eq? vptr end) #t)
					 ((not (pair? dptr)) #f)
					 ((equal? (car vptr) (car dptr))
					  (scan (cdr vptr) (cdr dptr)))
					 (else #f))))
			    (match (cdr pat) end dict fail succeed)
			    (fail)))
		      (error "Not a segment variable" pat))
		  (let try-seg ((end dat))
		    (if (if p (p dat end) #t)
			(match (cdr pat)
			       end
			       (extend-dict-segment vname dat end dict)
			       (lambda ()
				 (if (pair? end) (try-seg (cdr end)) (fail)))
			       succeed)
			(try-next)))))
	    (match-call vname dat dict fail succeed))
	(let try-seg ((end dat))
	  (if (if p (p dat end) #t)
	      (match (cdr pat)
		     end
		     dict
		     (lambda ()
		       (if (pair? end) (try-seg (cdr end)) (fail)))
		     succeed)
	      (try-next))))))

;;; Evaluation in matching and in instantiation

(define (check-predicate predicate-expression dictionary)
  (if (not (eq? predicate-expression 'none))
      (evaluate-element predicate-expression dictionary)
      #t))


(define (execute skeleton dictionary)
  (if (evaluation? skeleton)
      (evaluate
       `(let ((*dictionary* ',dictionary))
	  ,(evaluation-expression skeleton)))
      (instantiate skeleton dictionary)))

(define (instantiate skeleton dictionary)
  (cond ((not (pair? skeleton)) skeleton)
	((element-value? skeleton)
	 (evaluate-element (element-expression skeleton) dictionary))
	((segment-value? (car skeleton))
	 (evaluate-segment (segment-expression (car skeleton))
			   dictionary
			   (instantiate (cdr skeleton) dictionary)))
	(else
	 (cons (instantiate (car skeleton) dictionary)
	       (instantiate (cdr skeleton) dictionary)))))


(define (evaluate-element expression dictionary)
  (evaluate (substitute-in expression dictionary)))

(define (evaluate-segment expression dictionary tail)
  (append (evaluate (substitute-in expression dictionary))
	  tail))

(define (substitute-in expression dictionary)
  (let walk ((e expression))
    (if (pair? e)
	(cons (walk (car e)) (walk (cdr e)))
	(let ((v (match-lookup e dictionary)))
	  (if v
	      (list 'quote (element-from v))
	      e)))))

;;; Syntax of the expressions being manipulated:

;;; Rule syntax

(define (rule-pattern rule) (car rule))
(define (rule-predicate rule) (cadr rule))
(define (rule-skeleton rule) (caddr rule))


;;; Pattern matching syntax.

(define (arbitrary-element? pat)
  (and (pair? pat) (eq? (car pat) '?)))
(define (arbitrary-segment? pat)
  (and (pair? pat) (eq? (car pat) '??)))

(define (match-var-name pat)
  (if (null? (cdr pat)) #f (cadr pat)))

(define (match-var-restriction pat)
  (cond ((null? (cdr pat)) #f)
	((null? (cddr pat)) #f)
	(else (evaluate (caddr pat)))))


;;; Skeleton instantation syntax.

(define (evaluation? skel)
  (and (pair? skel) (eq? (car skel) 'evaluate)))
(define (evaluation-expression skel) (cadr skel))


(define (element-value? skel)
  (and (pair? skel) (eq? (car skel) ':)))
(define (element-expression skel) (cadr skel))

(define (segment-value? skel)
  (and (pair? skel) (eq? (car skel) '::)))
(define (segment-expression skel) (cadr skel))

(define (evaluate exp)
  (eval exp rule-environment))

;;; Dictionaries

(define the-empty-dictionary '())

(define match-lookup assq)

(define (segment-var? vcell)
  (not (null? (cddr vcell))))
(define segment-beg cadr)
(define segment-end caddr)

(define (empty-segment? seg)
  (eq? (segment-end seg)
       (segment-beg seg)))


(define (element-var? vcell)
  (null? (cddr vcell)))
(define element-in cadr)

(define (extend-dict-segment name beg end dict)
  (cons (list name beg end) dict))
(define (extend-dict-element name dat dict)
  (cons (list name dat) dict))

(define (element-from entry)
  (if (fix:= (length entry) 2)
      (element-in entry)
      (segment->list (segment-beg entry) (segment-end entry))))

(define (segment->list beg end)
  (define (append-seg point)
    (if (eq? point end)
	'()
	(cons (car point)
	      (append-seg (cdr point)))))
  (append-seg beg))


(define ((match-assign-element name dictionary value)
	 fail succeed)
  (succeed (extend-dict-element name value dictionary)
	   fail))
  
(define (match-get-value vname dict)
  (let ((v (match-lookup vname dict)))
    (if v
	(element-from v)
	(error "Unbound match variable" vname))))

(define (match-call call dat dict fail succeed)
  ((evaluate (car call)) dat dict fail succeed))

