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

;;;;                  UTILS.SCM
;;; A few utilities
;;; 9/15/89 gjs -- added FORALL, EXISTS; moved ACCUMULATION, INVERSE-ACCUMULATION here.
;;; 7/16/89 (mh) correcting bug in DEFAULT-LOOKUP
;;; 9/22/89 (gjs) reduce->a-reduce

(declare (usual-integrations))

(define (assert p #!optional error-comment irritant)
  (if (not p)
      (begin
	(if (not (default-object? irritant))
	    (pp irritant))
	(error (if (default-object? error-comment)
		   "Failed assertion"
		   error-comment)))))


(define (forall p? l)
  (let loop ((l l))
    (cond ((null? l) true)
	  ((p? (car l)) (loop (cdr l)))
	  (else false))))

(define (exists p? l)
  (let loop ((l l))
    (cond ((null? l) false)
	  ((p? (car l)) true)
	  (else (loop (cdr l))))))


(define (&or disjuncts)
  (cond ((null? disjuncts) false)
	((car disjuncts) true)
	(else (&or (cdr disjuncts)))))

(define (*or . disjuncts) (&or disjuncts))


(define (&and conjuncts)
  (cond ((null? conjuncts) true)
	((car conjuncts) (&and (cdr conjuncts)))
	(else false)))

(define (*and . conjuncts) (&and conjuncts))


(define (conjunction predicate1 predicate2)
  (lambda (x)
    (and (predicate1 x) (predicate2 x))))

(define (disjunction predicate1 predicate2)
  (lambda (x)
    (or (predicate1 x) (predicate2 x))))

(define (negation predicate)
  (lambda (x) (not (predicate x))))

(define (implication antecedent consequent)
  (lambda (x) (or (not (antecedent x)) (consequent x))))

(define (do-up low hi proc)
  ;; execute PROC for values beginning at LOW up to HI (exclusive)
  (if (fix:< low hi)
      (begin (proc low)
	     (do-up (fix:+ low 1) hi proc))))

(define (do-down hi low proc)
  ;; execute PROC for values beginning at HI down to LOW (exclusive)
  (if (fix:< low hi)
      (begin (proc hi)
	     (do-down (fix:- hi 1) low proc))))


;;; List utilities

(define (count-elements p? l)
  (let loop ((count 0) (l l))
    (cond ((null? l) count)
          ((p? (car l)) (loop (fix:+ count 1) (cdr l)))
          (else (loop count (cdr l))))))

(define (find-first pred lst)
  (cond ((null? lst) #f)
	((pred (car lst)) (car lst))
	(else (find-first pred (cdr lst)))))

(define (countsymbols exp)
  (cond ((pair? exp)
	 (fix:+ (countsymbols (car exp))
		(countsymbols (cdr exp))))
	((symbol? exp) 1)
	(else 0)))

(define (butlast l)
  (if (null? (cdr l)) 
      '()
      (cons (car l)
            (butlast (cdr l)))))

(define (last l)
  (car (last-pair l)))

(define (list-transpose l)
  (apply map list l))

(define (list-index-of x lst)
  (cond ((null? lst)
	 (error "Not in list -- LIST-INDEX-OF" x))
	((equal? x (car lst)) 0)
	(else (fix:+ (list-index-of x (cdr lst)) 1))))

(define (delete-nth n list)
  (if (fix:= n 0)
      (cdr list)
      (cons (car list)
	    (delete-nth (fix:- n 1) (cdr list)))))

(define ((list:elementwise proc) . lists)
  (apply map proc lists))

;;; Sets represented as unsorted lists

(define (list-adjoin item list)
  (if (member item list)
      list
      (cons item list)))

(define (list-union l1 l2)
  (cond ((null? l1) l2)
	((member (car l1) l2)
	 (list-union (cdr l1) l2))
	(else (cons (car l1)
		    (list-union (cdr l1) l2)))))

(define (list-intersection l1 l2)
  (cond ((null? l1) '())
	((member (car l1) l2)
	 (cons (car l1)
	       (list-intersection (cdr l1) l2)))
	(else (list-intersection (cdr l1) l2))))

(define (list-difference l1 l2)
  (cond ((null? l1) '())
	((member (car l1) l2)
	 (list-difference (cdr l1) l2))
	(else
	 (cons (car l1)
	       (list-difference (cdr l1) l2)))))

(define (duplications? lst)
  (cond ((null? lst) false)
	((member (car lst) (cdr lst)) true)
	(else (duplications? (cdr lst)))))

(define (remove-duplicates list)
  (if (null? list)
      '()
      (let ((rest (remove-duplicates (cdr list))))
        (if (member (car list) rest)
            rest
            (cons (car list) rest)))))

(define (subset? s1 s2)
  (if (null? s1)
      true
      (and (member (car s1) s2)
	   (subset? (cdr s1) s2))))

(define (same-set? s1 s2)
  (and (subset? s1 s2)
       (subset? s2 s1)))

;;; MAP-DISTINCT-PAIRS APPLYs a function, F, to every distinct pair of
;;; values chosen from the list, M, producing a list of the results.

(define (map-distinct-pairs f lst)
  (map (lambda (p) (apply f p))
       (distinct-pairs lst)))

(define (distinct-pairs lst)
  (if (null? lst)
      '()
      (let ((f (car lst))
	    (r (distinct-pairs (cdr lst))))
	(let loop ((left (cdr lst)))
	  (if (null? left)
	      r
	      (cons (list f (car left))
		    (loop (cdr left))))))))

(define (for-each-distinct-pair proc list)
  (if list
      (let loop ((first (car list)) (rest (cdr list)))
	(for-each (lambda (other-element)
		    (proc first other-element))
		  rest)
	(if rest (loop (car rest) (cdr rest))))))

;;; Elementary table utilities implemented in ALISTs

(define (lookup key table)
  (let ((val (assq key table)))
    (if val
	(cadr val)
	(error "key not in table -- LOOKUP" key))))

(define (rlookup key table)
  (cond ((null? table) false)
	((eq? key (cadar table)) (car table))
	(else (rlookup key (cdr table)))))

(define (rassq key table)
  (cond ((null? table) false)
	((eq? key (cdar table)) (car table))
	(else (rassq key (cdr table)))))

(define (rassoc key table)
  (cond ((null? table) false)
	((equal? key (cdar table)) (car table))
	(else (rassoc key (cdr table)))))

(define (disassoc key alist)
  (cond ((null? alist) '())
	((equal? key (caar alist))
	 (cdr alist))
	(else
	 (cons (car alist)
	       (disassoc key (cdr alist))))))


;;; Elementary table utility implemented as PLISTs

(define (default-lookup name default list)
  (let ((L (memq name list)))
    (if L (cadr L) default)))

(define (table-of is? keys values)
  (define (lookup key)
    (let next ((ks keys) (vs values))
      (cond ((null? ks)
             (error "Key not in table" key))
            ((is? key (car ks)) (car vs))
            (else (next (cdr ks) (cdr vs))))))
  lookup)

(define make-pairwise-test
  (lambda (pred)
    (lambda args
      (define (loop x y rem)
	(and (pred x y)
	     (or (null? rem)
		 (loop y (car rem) (cdr rem)))))
      (if (or (null? args) (null? (cdr args)))
	  (error "Pred needs 2 args" pred args)
	  (loop (car args) (cadr args) (cddr args))))))

(define (all-equal? lst)
  (define (lp lst)
    (if (null? (cdr lst))
	#t
	(and (equal? (car lst) (cadr lst))
	     (lp (cdr lst)))))
  (if (null? lst)
      #t
      (lp lst)))



(define accumulation
  (lambda (operation identity)
    (lambda rest
      (define (loop accum rem)
	(if (null? rem)
	    accum
	    (loop (operation accum (car rem)) (cdr rem))))
      (cond ((null? rest) identity)
	    ((null? (cdr rest)) (car rest))
	    (else (operation (car rest) (loop (cadr rest) (cddr rest))))))))

(define inverse-accumulation
  (lambda (operation1 operation2 invop identity)
    (lambda rest
      (define (loop accum rem)
	(if (null? rem)
	    accum
	    (loop (operation2 accum (car rem)) (cdr rem))))
      (cond ((null? rest) identity)
	    ((null? (cdr rest)) (invop (car rest)))
	    ((null? (cddr rest)) (operation1 (car rest) (cadr rest)))
	    (else (operation1 (car rest) (loop (cadr rest) (cddr rest))))))))


(define (left-circular-shift l)
  (if (or (null? l) (null? (cdr l)))
      l
      (append (cdr l) (list (car l)))))

(define (right-circular-shift l)
  (if (or (null? l) (null? (cdr l)))
      l
      (let ((r (reverse l)))
        (cons (car r) (reverse! (cdr r))))))

;;; Functional operators

;;; Arity is important to special case.

(define *at-least-zero* '(0 . #f))
(define *exactly-zero* '(0 . 0))
(define *at-least-one* '(1 . #f))
(define *exactly-one* '(1 . 1))
(define *at-least-two* '(2 . #f))
(define *exactly-two* '(2 . 2))
(define *at-least-three* '(3 . #f))
(define *exactly-three* '(3 . 3))
(define *one-or-two*    '(1 . 2))

(define (exactly-n? arity)
  (fix:= (car arity) (cdr arity)))

(define (any-number? arity)
  (and (null? (cdr arity))
       (fix:= (car arity) 0)))

(define (compose . fs)
  (compose-n fs))

(define (compose-n fs)
  (define (lp fs)
    (cond ((null? (cdr fs)) (car fs))
	  (else (compose-2 (car fs) (lp (cdr fs))))))
  (cond ((null? fs) identity)
	((null? (cdr fs)) (car fs))
	(else				;compose-bin preserves arity
	 (compose-bin (lp (butlast fs))
		      (car (last-pair fs))))))

(define (identity x) x)

(define (compose-2 f g)
  (cond ((pair? g)
	 (lambda x
	   (apply f
		  (map (lambda (gi)
			 (apply gi x))
		       g))))
	(else
	 (lambda x
	   (f (apply g x))))))

(define (compose-bin f g)
  (cond ((pair? g)
	 (let ((a
		(a-reduce joint-arity
			  (map procedure-arity g))))
	   (cond ((equal? a *at-least-zero*)
		  (lambda x
		    (apply f
			   (map
			    (lambda (gi)
			      (apply gi x))
			    g))))
		 ((equal? a *exactly-zero*)
		  (lambda ()
		    (apply f
			   (map (lambda (gi)
				  (gi))
				g))))
		 ((equal? a *at-least-one*)
		  (lambda (x . y)
		    (apply f
			   (map (lambda (gi)
				  (apply gi x y))
				g))))
		 ((equal? a *exactly-one*)
		  (lambda (x)
		    (apply f
			   (map (lambda (gi)
				  (gi x))
				g))))

		 ((equal? a *at-least-two*)
		  (lambda (x y . z)
		    (apply f
			   (map (lambda (gi)
				  (apply gi x y z))
				g))))
		 ((equal? a *exactly-two*)
		  (lambda (x y)
		    (apply f
			   (map (lambda (gi)
				  (gi x y))
				g))))

		 ((equal? a *at-least-three*)
		  (lambda (u x y . z)
		    (apply f
			   (map (lambda (gi)
				  (apply gi u x y z))
				g))))
		 ((equal? a *exactly-three*)
		  (lambda (x y z)
		    (apply f
			   (map (lambda (gi)
				  (gi x y z))
				g))))
		 ((equal? a *one-or-two*)
		  (lambda (x #!optional y)
		    (if (default-object? y)
			(apply f
			       (map (lambda (gi)
				      (gi x))
				    g))
			(apply f
			       (map (lambda (gi)
				      (gi x y))
				    g)))))
		 (else
		  (lambda x
		    (apply f
			   (map
			    (lambda (gi)
			      (apply gi x))
			    g)))))))
	(else
	 (let ((a (procedure-arity g)))
	   (cond ((equal? a *at-least-zero*)
		  (lambda x
		    (f (apply g x))))
		 ((equal? a *exactly-zero*)
		  (lambda ()
		    (f (g))))
		 ((equal? a *at-least-one*)
		  (lambda (x . y)
		    (f (apply g x y))))
		 ((equal? a *exactly-one*)
		  (lambda (x)
		    (f (g x))))
		 ((equal? a *at-least-two*)
		  (lambda (x y . z)
		    (f (apply g x y z))))
		 ((equal? a *exactly-two*)
		  (lambda (x y)
		    (f (g x y))))
		 ((equal? a *at-least-three*)
		  (lambda (u x y . z)
		    (f (apply g u x y z))))
		 ((equal? a *exactly-three*)
		  (lambda (x y z)
		    (f (g x y z))))
		 ((equal? a *one-or-two*)
		  (lambda (x #!optional y)
		    (if (default-object? y)
			(f (g x))
			(f (g x y)))))
		 (else
		  (lambda x
		    (f (apply g x)))))))))

(define (any? . args) #t)
(define (none? . args) #f)

(define ((constant x) . y) x)

(define (joint-arity a1 a2)
  (if (and a1 a2)
      (let ((amin (max (car a1) (car a2)))
	    (amax
	     (let ((a1max (cdr a1)) (a2max (cdr a2)))
	       (if a1max
		   (if a2max
		       (min a1max a2max)
		       a1max)
		   a2max))))
	(if (and amax (fix:< amax amin))
	    #f
	    (cons amin amax)))
      #f))

(define (a-reduce f l)
  (define (loop l)
     (if (null? (cdr l))
         (car l)
         (loop (cons (f (car l) (cadr l)) (cddr l)))))
  (if (null? l)
      (error "Reduce no elements")
      (loop l)))

(define (filter pred l)
  (let lp ((l l))
    (cond ((null? l) '())
          ((pred (car l)) (cons (car l) (lp (cdr l))))
          (else (lp (cdr l))))))

(define (make-map f)		; very neat, e.g. ((make-map -) '(3 2) '(1 1)) = '(2 1)
  (lambda x (apply map (cons f x))))

(define ((bracket . fl) . x)
  (map (lambda (f) (apply f x))
       fl))

(define ((apply-to-all f) x)
  (map f x))

(define (((nary-combine fnary) . fs) . xs)
  (apply fnary
	 (map (lambda (f) (apply f xs))
	      fs)))

(define (((binary-combine fbin) f1 f2) . xs)
  (fbin (apply f1 xs) (apply f2 xs)))

(define (((unary-combine funary) f) . xs)
  (funary (apply f xs)))

(define (iterated f n #!optional id)
  (if (fix:< n 0)
      (error "I don't know how to invert -- ITERATED" f n)
      (let ((ident (if (default-object? id) identity id)))
	(if (fix:= n 0)
	    ident
	    (let lp ((n n))
	      (if (fix:= n 1)
		  f
		  (compose-2 f (lp (fix:- n 1)))))))))


;;; Generalization of fixed point stuff

(define (iterate-until-stable f done? x0)
  (let lp ((x x0))
    (let ((nx (f x)))
      (if (done? nx x)
	  nx
	  (lp nx)))))


;;; given a function f of a variable number of arguments, return a new 
;;; function that accepts a single vector argument and calls f with the 
;;; vector elements as arguments

(define (make-function-of-vector f)
  (lambda (v)
    (apply f (vector->list v))))

;;; given a function of a single vector argument, return a new function 
;;; that takes multiple arguments, being the vector elements

(define (make-function-of-arguments f)
  (lambda args
    (f (list->vector args))))


;;; The following procedure came from SCHEME 6.1.2 RUNTIME

(define alphaless?
  (let ()
    (define (stringify object)
      (cond ((symbol? object) (symbol->string object))
	    ((string? object) object)
	    (else (error "ALPHALESS?: Wrong type argument" object))))

    (named-lambda (alphaless? x y)
      (string<? (stringify x) (stringify y)))))


(define (concatenate-names-maker concat-string)
  (define (cn strings)
    (cond ((null? strings) "")
	  ((null? (cdr strings)) (car strings))
	  (else
	   (a-reduce (lambda (n1 n2)
		       (string-append n1 concat-string n2))
		     strings))))
  (define (concatenate-names . names)
    (cond ((null? names) the-null-symbol)
	  ((null? (cdr names)) (car names))
	  (else
	   (string->symbol
	    (cn (map symbol->string
		     (filter (lambda (x)
			       (not (eq? x the-null-symbol)))
			     names)))))))
  concatenate-names)


(define the-null-symbol (string->symbol ""))

(define concatenate-names (concatenate-names-maker "."))

;;; Special property of MIT CScheme

(define (print-depth #!optional newval)
  (if (default-object? newval) (set! newval #F))
  (if (or (not newval)
	  (and (integer? newval)
	       (positive? newval)))
      (set! *unparser-list-depth-limit* newval)
      (error "PRINT-DEPTH: Wrong type argument" newval)))

(define (print-breadth #!optional newval)
  (if (default-object? newval) (set! newval #F))
  (if (or (not newval)
	  (and (integer? newval)
	       (positive? newval)))
      (set! *unparser-list-breadth-limit* newval)
      (error "PRINT-BREADTH: Wrong type argument" newval)))


;;;for printing things out

(define (wallp-pp p? . objs)
  (if p? (for-each pp objs)))

(define (pp-it x)
  (pp x)
  x)

(define (watch-it wallp message)
  (lambda (e)
    (if wallp
	(begin (newline)
	       (display message)
	       (pp e)))
    e))

(define (cpp x #!optional port)
  (pp x
      (if (default-object? port)
	  (current-output-port)
	  port)
      ;; as code
      true))
