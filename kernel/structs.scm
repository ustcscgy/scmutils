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

;;;;                      Structures

(declare (usual-integrations))

;;; Structures are primitive tensor-like objects.  They are
;;; represented as recursive combinations of row vectors and column
;;; vectors, useful for dealing with derivatives of things with
;;; structured inputs and outputs.


#| in TYPES.SCM

(define (column? x)
  ;;(and (pair? x) (eq? (car x) column-type-tag))
  (vector? x))

(define (row? x)
  (and (pair? x)
       (eq? (car x) row-type-tag)))

(define (structure? x)
  (or (column? x) (row? x)))


(define (abstract-structure? x)
  (or (abstract-column? x) (abstract-row? x)))

|#

(define (s:type v)
  (cond ((column? v) column-type-tag)
	((row? v) row-type-tag)
	((abstract-column? v) column-type-tag)
	((abstract-row? v) row-type-tag)
	(else
	 (error "Bad structure -- S:TYPE" v))))

(define (sc:type-predicate v) column-quantity?)
(define (sr:type-predicate v) row-quantity?)


(define (vector->column v)
  ;;(list column-type-tag v)
  v)

(define (vector->row v)
  (list row-type-tag v))


(define (vector->up v)
  ;;(list column-type-tag v)
  v)

(define (vector->down v)
  (list row-type-tag v))


(define (s:structure up/down v)
  (case up/down
    ((up contravariant vector column)
     (vector->column v))
    ((down covariant covector row)
     (vector->row v))
    (else
     (error "Bad up/down spec -- S:STRUCTURE"
	    up/down v))))


(define (column->vector v)
  ;;(cadr v)
  v)

(define (row->vector v)
  (cadr v))


(define (up->vector v)
  ;;(cadr v)
  v)

(define (down->vector v)
  (cadr v))

(define (s:->vector v)
  (cond ((column? v) (column->vector v))
	((row? v) (row->vector v))
	(else
	 (error "Bad structure -- S:->VECTOR" v))))


(define (column . args)
  (vector->column (list->vector args)))

(define (row . args)
  (vector->row (list->vector args)))

(define (up . args)
  (vector->column (list->vector args)))

(define (down . args)
  (vector->row (list->vector args)))

#|
;;; Defined in types.scm
(define (up? v) (column? v))
(define (down? v) (row? v))
|#

(define (s:opposite v)
  (cond ((column? v) 'down)
	((row? v) 'up)
	(else
	 (error "Bad structure -- S:OPPOSITE" v))))

(define (s:same v)
  (cond ((column? v) 'up)
	((row? v) 'down)
	(else
	 (error "Bad structure -- S:SAME" v))))

(define (s:length v)
  (if (structure? v)
      (vector-length (s:->vector v))
      1))


(define (s:ref v i)
  (if (structure? v)
      (vector-ref (s:->vector v) i)
      (if (fix:= i 0)
	  v
	  (error "Bad structure -- S:REF" v i))))

(define (s:with-substituted-coord v i xi)
  (if (structure? v)
      (s:structure (s:same v)
		   (vector-with-substituted-coord
		    (s:->vector v)
		    i xi))
      (if (fix:= i 0)
	  xi
	  (error "Bad structure -- S:WITH-SUBSTITUTED-COORD" v i xi))))

(define (s:subst struct newval . chain)
  (s:subst-internal struct newval chain))

(define (s:subst-internal struct newval chain)
  (let lp ((chain chain) (struct struct))
    (if (null? chain)
        newval
        (s:generate (s:length struct)
                    (s:same struct)
                    (lambda (i)
                      (if (fix:= i (car chain))
                          (lp (cdr chain) (s:ref struct i))
                          (s:ref struct i)))))))

(define (s:generate n up/down proc)
  (s:structure up/down (v:generate n proc)))

(define (s:forall p s)
  (let ((n (s:length s)))
    (let lp ((i 1) (ans (p (s:ref s 0))))
      (cond ((fix:= i n) ans)
	    ((not ans) ans)
	    (else
	     (lp (fix:+ i 1) (p (s:ref s i))))))))  

(define (s:select . selectors)
  (let lp ((selectors selectors)
	   (ans g:identity)) 
    (if (null? selectors)
	ans
	(lp (cdr selectors)
	    (compose (lambda (s)
		       (s:ref s (car selectors)))
		     ans)))))


(define (s:map-chain proc s)
  (define (walk s rev-chain)
    (if (structure? s)
	(s:generate (s:length s)
		    (s:same s)
		    (lambda (i)
		      (walk (s:ref s i)
			    (cons i rev-chain))))
	(proc s (reverse rev-chain))))
  (walk s '()))

#|
(pe (s:map-chain (up 'a 'b 'c) cons))
(up (a 0) (b 1) (c 2))

(pe (s:map-chain (up 'a (down 'b 'c) 'd) cons))
(up (a 0) (down (b 1 0) (c 1 1)) (d 2))
|#

;;; S:FRINGE recursively traverses a structure, making up a list of
;;; the terminal elements.

(define (s:fringe s)
  (define (walk s ans)
    (if (structure? s)
	(let ((n (s:length s)))
	  (let lp ((i 0) (ans ans))
	    (if (fix:= i n)
		ans
		(lp (fix:+ i 1)
		    (walk (s:ref s i) ans)))))
	(cons s ans)))
  (walk s '()))

(define (s:foreach proc s)
  (define (walk s)
    (if (structure? s)
	(let ((n (s:length s)))
	  (let lp ((i 0))
	    (if (fix:= i n)
		'done
		(begin (walk (s:ref s i))
		       (lp (fix:+ i 1))))))
	(proc s)))
  (walk s))

;;; The following mappers only make sense if, when there is more than
;;; one structure they are all isomorphic.

(define (s:map/r proc . structures)
  (s:map/r/l proc structures))

(define (s:map/r/l proc structures)
  (s:map/l (lambda elements
	     (if (structure? (car elements))
		 (s:map/r/l proc elements)
		 (apply proc elements)))
	   structures))

(define (s:map proc . structures)
  (s:map/l proc structures))

(define (s:map/l proc structures)
  (if (structure? (car structures))
      (s:generate (s:length (car structures))
		  (s:same (car structures))
		  (lambda (i)
		    (apply proc
			   (map (lambda (s)
				  (s:ref s i))
				structures))))
      (apply proc structures)))

#|
(define (s:map/l proc structures)
  (s:generate (s:length (car structures))
	      (s:same (car structures))
	      (lambda (i)
		(apply proc
		       (map (lambda (s)
			      (s:ref s i))
			    structures)))))
|#

(define ((s:elementwise proc) . structures)
  (s:map/l proc structures))

(define structure:elementwise s:elementwise)

;;; Is there a part of thing that the predicate is true of?

(define (rexists pred thing)
  (let tlp ((thing thing))
    (cond ((pred thing) #t)
	  ((vector? thing)
	   (let ((n (vector-length thing)))
	     (let lp ((i 0))
	       (cond ((fix:= i n) #f)
		     ((tlp (vector-ref thing i)) #t)
		     (else (lp (fix:+ i 1)))))))
	  ((structure? thing)
	   (let ((n (s:length thing)))
	     (let lp ((i 0))
	       (cond ((fix:= i n) #f)
		     ((tlp (s:ref thing i)) #t)
		     (else (lp (fix:+ i 1)))))))
	  ((matrix? thing)
	   (tlp (matrix->array thing)))
	  ((pair? thing)
	   (cond ((memq (car thing) type-tags)
		  (let ((v (get-property thing 'expression)))
		    (if (not v)
			#f
			(tlp v))))
		 ((list? thing)
		  (there-exists? thing tlp))
		 (else (tlp (cdr thing)))))
	  (else #f))))

(define (s:arity v) (v:arity (s:->vector v)))

(define (s:inexact? v) (v:inexact? (s:->vector v)))

(define (s:zero? v) (v:zero? (s:->vector v)))

(define (s:unary v:op v)
  (s:structure (s:same v) (v:op (s:->vector v))))

(define (s:zero-like v) (s:unary v:zero-like v))

(define (s:negate v) (s:unary v:negate v))

(define (s:magnitude v) (complex-norm (s:->vector v)))
(define (s:abs v) (euclidean-norm (s:->vector v)))
(define (s:conjugate v) (s:unary v:conjugate v))

(define (structure=structure v1 v2)
  (and (eq? (s:same v1) (s:same v2))
       (fix:= (s:length v1) (s:length v2))
       (vector=vector (s:->vector v1)
		      (s:->vector v2))))

(define (s:binary v:op v1 v2)
  (let ((up/down (s:same v1)))
    (assert (eq? up/down (s:same v2)))
    (assert (fix:= (s:length v1) (s:length v2)))
    (s:structure up/down
		 (v:op (s:->vector v1)
		       (s:->vector v2)))))

(define (structure+structure v1 v2)
  (s:binary vector+vector v1 v2))

(define (structure-structure v1 v2)
  (s:binary vector-vector v1 v2))

#|
(define (s:multiply v1 v2)
  (if (s:compatible-for-contraction? v1 v2)
      (v:dot-product (s:->vector v1) (s:->vector v2))
      (begin
	(if (not *allowing-incompatible-multiplication*)
	    (bkpt "Incompatible multiplication" v1 v2))
	(s:generate (s:length v2) (s:same v2)
		    (lambda (i)
		      (g:* v1 (s:ref v2 i)))))))
|#

;;; Want to allow matrix multiply too...

(define (s:multiply v1 v2)
  (cond ((s:compatible-for-contraction? v1 v2)
	 (v:dot-product (s:->vector v1) (s:->vector v2)))
	((or *allowing-incompatible-multiplication*
	     (and (or (and (down? v1) (down? v2))
		      (and (up? v1) (up? v2)))
		  (s:forall (lambda (c)
			      (s:compatible-for-contraction? v1 c))
			    v2)))
	 (s:generate (s:length v2) (s:same v2)
		     (lambda (i)
		       (g:* v1 (s:ref v2 i)))))
	(else
	 (bkpt "Incompatible multiplication" v1 v2))))

(define *allowing-incompatible-multiplication* #t)



(define (s:compatible-for-contraction? v1 v2)
  (or (and (down? v1) (up? v2)
	   (s:compatable-elements? v1 v2))
      (and (up? v1) (down? v2)
	   (s:compatable-elements? v1 v2))))

(define (s:compatable-elements? v1 v2)
  (let ((n (s:length v1)))
    (and (fix:= n (s:length v2))
	 (let lp ((i 0))
	   (cond ((fix:= i n) #t)
		 ((or (not (structure? (s:ref v1 i)))
		      (not (structure? (s:ref v2 i))))
		  (lp (fix:+ i 1)))
		 ((s:compatible-for-contraction? (s:ref v1 i)
						 (s:ref v2 i))
		  (lp (fix:+ i 1)))
		 (else #f))))))

#|
;;; down-down rule is designed to effect matrix multiply.

(print-expression
 (matrix*matrix
  (matrix-by-rows '(a b) '(c d))
  (matrix-by-rows '(e f) '(g h))))
(matrix-by-rows (list (+ (* a e) (* b g)) (+ (* a f) (* b h)))
                (list (+ (* c e) (* d g)) (+ (* c f) (* d h))))

(print-expression
 (s:multiply
  (down (up 'a 'c) (up 'b 'd))
  (down (up 'e 'g) (up 'f 'h))))
(down (up (+ (* a e) (* b g)) (+ (* c e) (* d g)))
      (up (+ (* a f) (* b h)) (+ (* c f) (* d h))))

(print-expression
 (s:multiply
  (down (up 'a 'c) (up 'b 'd))
  (up 'dx 'dy)))
(up (+ (* a dx) (* b dy)) (+ (* c dx) (* d dy)))

|#


;;; Given two structures their outer product makes a structure
#|
(define (s:outer-product s1 s2)
  (let lp ((s s2))
    (if (structure? s)
	(s:generate (s:length s) (s:same s)
		    (lambda (i) (lp (s:ref s i))))
	(g:* s1 s))))
#|
(pe (s:outer-product (up 'a0 'a1)
		   (down 'b0 'b1 'b2)))
(down (up (* a0 b0) (* a1 b0))
      (up (* a0 b1) (* a1 b1))
      (up (* a0 b2) (* a1 b2)))
|#
|#

(define (s:outer-product struct2 struct1)
  (s:map/r (lambda (s1)
	     (s:map/r (lambda (s2)
			(g:* s1 s2))
		      struct2))
	   struct1))
#|
(pe (s:outer-product (up 'a0 'a1)
		   (down 'b0 'b1 'b2)))
(down (up (* a0 b0) (* a1 b0))
      (up (* a0 b1) (* a1 b1))
      (up (* a0 b2) (* a1 b2)))

;;; cf.

(pe (m:outer-product (column-matrix 'a0 'a1)
		     (row-matrix 'b0 'b1 'b2)))
(matrix-by-rows (list (* a0 b0) (* a0 b1) (* a0 b2))
		(list (* a1 b0) (* a1 b1) (* a1 b2)))



(let* ((s (up 'dt (up 'dx 'dy) (down 'dpx 'dpy)))
       (s* (compatible-shape s)))
  (* s* (s:outer-product s s*) s))
#|
(+ (* (expt dpx 2) (expt x29083 2))
   (* 2 dpx dpy x29083 x29084)
   (* 2 dpx dt x29080 x29083)
   (* 2 dpx dx x29081 x29083)
   (* 2 dpx dy x29082 x29083)
   (* (expt dpy 2) (expt x29084 2))
   (* 2 dpy dt x29080 x29084)
   (* 2 dpy dx x29081 x29084)
   (* 2 dpy dy x29082 x29084)
   (* (expt dt 2) (expt x29080 2))
   (* 2 dt dx x29080 x29081)
   (* 2 dt dy x29080 x29082)
   (* (expt dx 2) (expt x29081 2))
   (* 2 dx dy x29081 x29082)
   (* (expt dy 2) (expt x29082 2)))
|#
|#

(define (structure:expt t n)
  (cond ((fix:= n 1) t)
	((fix:> n 1) (g:* t (structure:expt t (fix:- n 1))))
	(else (error "Cannot: " `(expt ,t ,n)))))


(define (scalar*structure s v)
  (s:structure (s:same v)
	       (scalar*vector s (s:->vector v))))

(define (structure*scalar v s)
  (s:structure (s:same v)
	       (vector*scalar (s:->vector v) s)))


(define (structure/scalar v s)
  (s:structure (s:same v)
	       (vector/scalar (s:->vector v) s)))


(define (s:square v)
  (let ((vv (s:->vector v)))
    (let ((n (vector-length vv)))
      (if (fix:= n 0)
	  :zero
	  (let lp ((i 1) (sum (g:square (vector-ref vv 0))))
	    (if (fix:= i n)
		sum
		(lp (fix:+ i 1)
		    (g:+ sum (g:square (vector-ref vv i))))))))))

(define (s:dot-product v1 v2)
  (if (not (and (eq? (s:same v1) (s:same v2)) (= (s:dimension v1) (s:dimension v2))))
      (error "Incompatible structures -- S:DOT-PRODUCT" v1 v2))
  (apply g:+ (map g:* (s:fringe v1) (s:fringe v2))))


(define (s:partial-derivative v varspecs)
  (s:structure (s:same v)
	       (v:partial-derivative (s:->vector v) varspecs)))

(define (s:apply v args)
  (s:structure (s:same v)
	       (v:apply (s:->vector v) args)))

(assign-operation 'type                s:type            structure?)
(assign-operation 'type-predicate      sc:type-predicate column?)
(assign-operation 'type-predicate      sr:type-predicate row?)
(assign-operation 'arity               s:arity           structure?)
(assign-operation 'inexact?            s:inexact?        structure?)

(assign-operation 'zero?               s:zero?           structure?)
(assign-operation 'zero-like           s:zero-like       structure?)
(assign-operation 'negate              s:negate          structure?)

(assign-operation 'magnitude           s:magnitude       structure?)
(assign-operation 'abs                 s:abs             structure?)
(assign-operation 'conjugate           s:conjugate       structure?)


(assign-operation '=           structure=structure     structure? structure?)
(assign-operation '+           structure+structure     structure? structure?)
(assign-operation '-           structure-structure     structure? structure?)

(assign-operation '*           s:multiply              structure?    structure?)

(assign-operation 'square      s:square                structure?)

(assign-operation 'dot-product s:dot-product           structure?    structure?)
(assign-operation 'outer-product s:outer-product       structure?    structure?)

(assign-operation '*           scalar*structure        scalar?    structure?)
(assign-operation '*           structure*scalar        structure? scalar?)

(assign-operation '/           structure/scalar        structure? scalar?)

(assign-operation 'expt        structure:expt          structure? exact-integer?)


#| ;;; Should be subsumed by deriv:pd in deriv.scm.

(assign-operation 'partial-derivative
		  s:partial-derivative
		  structure? any?)
|#
(assign-operation 'apply       s:apply           structure? any?)

;;; Abstract structures generalize structural quantities.

(define (literal-column symbol)
  (make-literal column-type-tag symbol))

(define (literal-row symbol)
  (make-literal abstract-row-type-tag symbol))

(define (literal-up-tuple symbol)
  (make-literal column-type-tag symbol))

(define (literal-down-tuple symbol)
  (make-literal abstract-row-type-tag symbol))

(define (as:arity v)
  ;; Default is vector of numbers.
  (get-property v 'arity *at-least-zero*))

(define (ac:zero-like v)
  (let ((z (literal-column (list 'zero-like v))))
    (add-property! z 'zero #t)
    z))

(define (ar:zero-like v)
  (let ((z (literal-row (list 'zero-like v))))
    (add-property! z 'zero #t)
    z))

(define (make-column-combination operator #!optional reverse?)
  (if (default-object? reverse?)
      (lambda operands 
	(make-combination column-type-tag
			  operator operands))
      (lambda operands 
	(make-combination column-type-tag
			  operator (reverse operands)))))

(define (make-row-combination operator #!optional reverse?)
  (if (default-object? reverse?)
      (lambda operands 
	(make-combination abstract-row-type-tag
			  operator operands))
      (lambda operands 
	(make-combination abstract-row-type-tag
			  operator (reverse operands)))))

(assign-operation 'type           s:type             abstract-structure?)
(assign-operation 'type-predicate sc:type-predicate  abstract-column?)
(assign-operation 'type-predicate sr:type-predicate  abstract-row?)
(assign-operation 'arity          as:arity           abstract-structure?)
(assign-operation 'inexact? (has-property? 'inexact) abstract-structure?)

(assign-operation 'zero?    (has-property? 'zero)    abstract-structure?)
(assign-operation 'zero-like      ac:zero-like       abstract-column?)
(assign-operation 'zero-like      ar:zero-like       abstract-row?)

(assign-operation
   'negate     (make-column-combination 'negate)     abstract-column?)
(assign-operation
   'negate     (make-row-combination 'negate)        abstract-row?)

(assign-operation
   'magnitude  (make-column-combination 'magnitude)  abstract-column?)
(assign-operation
   'magnitude  (make-row-combination 'magnitude)     abstract-row?)

(assign-operation
   'abs        (make-column-combination 'abs)        abstract-column?)
(assign-operation
   'abs        (make-row-combination 'abs)           abstract-row?)

(assign-operation
   'conjugate  (make-column-combination 'conjugate)  abstract-column?)
(assign-operation
   'conjugate  (make-row-combination 'conjugate)     abstract-row?)


;(assign-operation '= structure=structure abstract-structure? abstract-structure?)

(assign-operation
   '+  (make-vector-combination '+) abstract-column? abstract-column?)
(assign-operation
   '+  (make-column-combination '+) column?  abstract-column?)
(assign-operation
   '+  (make-column-combination '+ 'r) abstract-column? column?)

(assign-operation
   '+  (make-row-combination '+)    row?   abstract-row?)
(assign-operation
   '+  (make-row-combination '+ 'r) abstract-row? row?)

(assign-operation
   '-  (make-vector-combination '-) abstract-column? abstract-column?)
(assign-operation
   '-  (make-column-combination '-) column?  abstract-column?)
(assign-operation
   '-  (make-column-combination '-) abstract-column? column?)

(assign-operation
   '-  (make-row-combination '-)    row?   abstract-row?)
(assign-operation
   '-  (make-row-combination '-)    abstract-row? row?)

(assign-operation
   '*  (make-numerical-combination '*)    abstract-structure? abstract-structure?)
(assign-operation
   '*  (make-numerical-combination '*)    structure?          abstract-structure?)
(assign-operation
   '*  (make-numerical-combination '*)    abstract-structure? structure?)

(assign-operation
   '*  (make-column-combination '*)    scalar? abstract-column?)
(assign-operation
   '*  (make-column-combination '* 'r) abstract-column? scalar?)

(assign-operation
   '*  (make-row-combination '*)       scalar?    abstract-row?)
(assign-operation
   '*  (make-row-combination '* 'r)    abstract-row?    scalar?)

		     
(assign-operation
   '/  (make-column-combination '/)    abstract-column? scalar?)

(assign-operation
   '/  (make-row-combination '/)       abstract-row?    scalar?)

(assign-operation 'partial-derivative
		  (make-column-combination 'partial-derivative)
		  abstract-column? any?)

(assign-operation 'partial-derivative
		  (make-row-combination 'partial-derivative)
		  abstract-row? any?)

;;; Sometimes a structure may be usefully represented as a matrix, and
;;; sometimes a matrix should be seen as a structure.  S:CANONICALIZE
;;; attempts to make things matrices when they should be.  A matrix is
;;; needed when we have something that is really a row of columns.

(define (s:canonicalize s)
  (if (and (row? s) (column? (s:ref s 0)))
      (let ((nrows (s:length (s:ref s 0))))
	(if (s:forall (lambda (c)
			(and (column? c)
			     (fix:= (s:length c) nrows)))
		      s)
	    (m:generate nrows (s:length s)
			(lambda (i j)
			  (s:ref (s:ref s j) i)))
	    s))
      s))


;;; An argument list really wants to be represented as an (up) vector
;;; of arguments.  Also, any matrix in the argument list wants to be
;;; converted to a row of columns.

(define (list->up-structure lst)
  (vector->column
   (list->vector (map matrix->structure lst))))

(define (matrix->structure mat)
  (cond ((row? mat) mat)
	((column? mat) mat)
	((matrix? mat)
	 (s:generate (m:num-cols mat) 'row
		     (lambda (j)
		       (s:generate (m:num-rows mat) 'column
				   (lambda (i)
				     (matrix-ref mat i j))))))
	(else mat)))

(define (submatrix s lowrow hirow+1 lowcol hicol+1)
  (cond ((structure? s)
	 (m:submatrix (structure->matrix s) lowrow hirow+1 lowcol hicol+1))
	((matrix? s)
	 (m:submatrix s lowrow hirow+1 lowcol hicol+1))
	(else (error "Wrong type submatrix" s))))

#|
(define (up-structure->list s)
  (map s:canonicalize
       (vector->list s)))
|#

(define (up-structure->list s)
  (vector->list (up->vector s)))

;;; Sometimes a 2-tensor must be viewed as a matrix for some purpose,
;;; for example to invert it.  The following are the required coercions 
;;; between tensor structures and matrices.

;;; Convention for A^m_n: rightmost index, n, is length of outermost
;;; structure.

;;; a row of n rows each m long -> n rows X m columns

(define (A_mn->Mnm s)			
  (if (and (row? s) (row? (s:ref s 0)))
      (let ((nrows (s:length s))
	    (ncols (s:length (s:ref s 0))))
	(if (s:forall (lambda (r)
			(and (row? r)
			     (fix:= (s:length r) ncols)))
		      s)
	    (m:generate nrows ncols
			(lambda (i j)
			  (s:ref (s:ref s i) j)))
	    (error "Not A_mn -- A_mn->Mnm" s)))
      (error "Not A_mn -- A_mn->Mnm" s)))
#|
(A_mn->Mnm (row (row 'a 'b) (row 'c 'd) (row 'e 'f)))
;Value: (*matrix* (3 . 2) #(#(a b) #(c d) #(e f)))
|#

(define (Mnm->A_mn mat)
  (assert (matrix? mat) "Not a matrix -- Mnm->A_mn" mat)
  (s:generate (m:num-rows mat) 'row
	      (lambda (i)
		(s:generate (m:num-cols mat) 'row
			    (lambda (j)
			      (matrix-ref mat i j))))))

#|
((compose mnm->a_mn a_mn->mnm) (row (row 'a 'b) (row 'c 'd) (row 'e 'f)))
#|
(down (down a b) (down c d) (down e f))
|#
|#

;;; a col of n cols each m long -> m rows X n columns

(define (A^mn->Mmn s)
  (if (and (column? s) (column? (s:ref s 0)))
      (let ((ncols (s:length s))
	    (nrows (s:length (s:ref s 0))))
	(if (s:forall (lambda (r)
			(and (column? r)
			     (fix:= (s:length r) nrows)))
		      s)
	    (m:generate nrows ncols
			(lambda (i j)
			  (s:ref (s:ref s j) i)))
	    (error "Not A^mn -- A^mn->Mmn" s)))
      (error "Not A^mn -- A^mn->Mmn" s)))
#|
(A^mn->Mmn (column (column 'a 'b) (column 'c 'd) (column 'e 'f)))
;Value: (*matrix* (2 . 3) #(#(a c e) #(b d f)))
|#

(define (Mmn->A^mn mat)
  (assert (matrix? mat) "Not a matrix -- Mmn->A^mn" mat)
  (s:generate (m:num-cols mat) 'column
	      (lambda (j)
		(s:generate (m:num-rows mat) 'column
			    (lambda (i)
			      (matrix-ref mat i j))))))
#|
((compose mmn->a^mn a^mn->mmn) (column (column 'a 'b) (column 'c 'd) (column 'e 'f)))
#|
(up (up a b) (up c d) (up e f))
|#
|#

;;; a row of n cols each m long -> m rows X n columns

(define (A^m_n->Mmn s)
  (if (and (row? s) (column? (s:ref s 0)))
      (let ((nrows (s:length (s:ref s 0))))
	(if (s:forall (lambda (c)
			(and (column? c)
			     (fix:= (s:length c) nrows)))
		      s)
	    (m:generate nrows (s:length s)
			(lambda (i j)
			  (s:ref (s:ref s j) i)))
	    (error "Not A^m_n -- A^m_n->Mmn" s)))
      (error "Not A^m_n -- A^m_n->Mmn" s)))
#|
(A^m_n->Mmn (row (column 'a 'b) (column 'c 'd) (column 'e 'f)))
;Value: (*matrix* (2 . 3) #(#(a c e) #(b d f)))
|#

(define (Mmn->A^m_n mat)
  (assert (matrix? mat) "Not a matrix -- Mmn->A^m_n" mat)
  (s:generate (m:num-cols mat) 'row
	      (lambda (j)
		(s:generate (m:num-rows mat) 'column
			    (lambda (i)
			      (matrix-ref mat i j))))))
#|
((compose mmn->A^m_n A^m_n->mmn) (row (column 'a 'b) (column 'c 'd) (column 'e 'f)))
#|
(down (up a b) (up c d) (up e f))
|#
|#

;;; a col of n rows each m long -> n rows X m columns

(define (A_m^n->Mnm s)
  (if (and (column? s) (row? (s:ref s 0)))
      (let ((ncols (s:length (s:ref s 0))))
	(if (s:forall (lambda (c)
			(and (row? c)
			     (fix:= (s:length c) ncols)))
		      s)
	    (m:generate (s:length s) ncols
			(lambda (i j)
			  (s:ref (s:ref s i) j)))
	    (error "Not A_m^n -- A_m^n->Mmn" s)))
      (error "Not A_m^n -- A_m^n->Mmn" s)))
#|
(A_m^n->Mnm (column (row 'a 'b) (row 'c 'd) (row 'e 'f)))
;Value: (*matrix* (3 . 2) #(#(a b) #(c d) #(e f)))
|#

(define (Mnm->A_m^n mat)
  (assert (matrix? mat) "Not a matrix -- Mnm->A_m^n" mat)
  (s:generate (m:num-rows mat) 'column
	      (lambda (i)
		(s:generate (m:num-cols mat) 'row
			    (lambda (j)
			      (matrix-ref mat i j))))))
#|
((compose mnm->A_m^n A_m^n->mnm) (column (row 'a 'b) (row 'c 'd) (row 'e 'f)))
#|
(up (down a b) (down c d) (down e f))
|#
|#

;;; A few lonely tensor operations here -- this will expand later.

(define (2-down? s)
  (and (row? s) (row? (s:ref s 0))))

(define (2-up? s)
  (and (column? s) (column? (s:ref s 0))))

(define (up-of-downs? s)
  (and (column? s) (row? (s:ref s 0))))

(define (down-of-ups? s)
  (and (row? s) (column? (s:ref s 0))))

(define (2-tensor? s)
  (or (2-down? s) (2-up? s) (up-of-downs? s) (down-of-ups? s)))

(define (single-layer-down? s)
  (and (down? s)
       (not (there-exists? (vector->list (s:->vector s)) structure?))))

(define (single-layer-up? s)
  (and (up? s)
       (not (there-exists? (vector->list (s:->vector s)) structure?))))


(define (structure->matrix s)
  (cond ((2-down? s) (A_mn->Mnm s))
	((2-up? s) (A^mn->Mmn s))
	((up-of-downs? s) (A_m^n->Mnm s))
	((down-of-ups? s) (A^m_n->Mmn s))
	(else (error "structure->matrix" s))))

(define (s:invert s)
  (cond ((2-down? s)
	 (Mmn->A^mn (m:invert (m:transpose (A_mn->Mnm s)))))
	((2-up? s)
	 (Mnm->A_mn (m:invert (m:transpose (A^mn->Mmn s)))))
	((up-of-downs? s)
	 (Mnm->A_m^n (m:invert (A_m^n->Mnm s))))
	((down-of-ups? s)
	 (Mmn->A^m_n (m:invert (A^m_n->Mmn s))))
	(else (error "s:invert" s))))

#|
;;;; Test by equation solving.  All answers should be <e f>.

;;; down of downs
(let ((a (down (down 'a 'b) (down 'c 'd)))
      (b (up 'e 'f)))
  (* a (* (s:invert a) b)))
#|
(up e f)
|#

;;; up of ups
(let ((a (up (up 'a 'b) (up 'c 'd)))
      (b (down 'e 'f)))
  (* a (* (s:invert a) b)))
#|
(down e f)
|#

;;; up of downs
(let ((a (up (down 'a 'b) (down 'c 'd)))
      (b (down 'e 'f)))
  (* a (* (s:invert a) b)))
#|
(down e f)
|#

;;; down of ups
(let ((a (down (up 'a 'b) (up 'c 'd)))
      (b (up 'e 'f)))
  (* a (* (s:invert a) b)))
#|
(up e f)
|#
|#

(define (scalar/tensor x s)
  (g:* x (s:invert s)))

(define (s:solve-rs rv s)
  (g:* rv (s:invert s)))

(define (s:solve-cs cv s)
  (g:* (s:invert s) cv))

(define (s:solve-du-cs cv s)
  (g:* (s:invert s) cv))

(define (s:solve-ud-rs rv s)
  (g:* rv (s:invert s)))

(define (s:transpose2 s)
  (cond ((2-down? s)
	 (Mnm->A_mn (m:transpose (A_mn->Mnm s))))
	((2-up? s)
	 (Mmn->A^mn (m:transpose (A^mn->Mmn s))))
	((up-of-downs? s)
	 (Mmn->A^m_n (A_m^n->Mnm s)))
	((down-of-ups? s)
	 (Mnm->A_m^n (A^m_n->Mmn s)))
	(else (error "s:transpose2" s))))

(define (s:transpose-up->down s)
  (vector->down (up->vector s)))

(define (s:transpose-down->up s)
  (vector->up (down->vector s)))

#|
(define (transpose-test left-multiplier thing right-multiplier)
  ;; Should produce numerical zero and a zero structure
  (list (- (* left-multiplier (* thing right-multiplier))
	   (* (* (s:transpose2 thing) left-multiplier) right-multiplier))
	(- (s:transpose left-multiplier thing right-multiplier)
	   (s:transpose2 thing))))

;;; down down
(transpose-test (up 'a 'b)
		(down (down 'c 'd) (down 'e 'f) (down 'g 'h))
		(up 'i 'j 'k))
#| (0 (down (down 0 0 0) (down 0 0 0))) |#

;;; up up
(transpose-test (down 'a 'b)
		(up (up 'c 'd) (up 'e 'f) (up 'g 'h))
		(down 'i 'j 'k))
(0 (up (up 0 0 0) (up 0 0 0)))

;;; up down
(transpose-test (up 'a 'b)
		(up (down 'c 'd) (down 'e 'f) (down 'g 'h))
		(down 'i 'j 'k))
#|
(0 (down (up 0 0 0) (up 0 0 0)))
|#

;;; down up
(transpose-test (down 'a 'b)
		(down (up 'c 'd) (up 'e 'f) (up 'g 'h))
		(up 'i 'j 'k))
#|
(0 (up (down 0 0 0) (down 0 0 0)))
|#
|#

(assign-operation 'invert             s:invert                 2-tensor?)

(assign-operation 'transpose          s:transpose2             2-tensor?)
(assign-operation 'transpose          s:transpose-up->down     up?)
(assign-operation 'transpose          s:transpose-down->up     down?)

(assign-operation '/   scalar/tensor  scalar?                  2-tensor?)

(assign-operation '/   s:solve-rs     row?                     2-down?)
(assign-operation '/   s:solve-cs     column?                  2-up?)
(assign-operation '/   s:solve-ud-rs  row?                     up-of-downs?)
(assign-operation '/   s:solve-du-cs  column?                  down-of-ups?)

(define (s:determinant s)
  (m:determinant (structure->matrix s)))

(define (s:trace s)
  (m:trace (structure->matrix s)))

(assign-operation 'determinant s:determinant 2-tensor?)
(assign-operation 'trace       s:trace       2-tensor?)

;;; This just changes the up and downness of a structure

(define (flip-indices s)
  (if (structure? s)
      (s:generate (s:length s)
		  (s:opposite s)
		  (lambda (i)
		    (flip-indices (s:ref s i))))
      s))
#|
(flip-indices (up 'a (up 'b 'c) (down 'd (up 'e 'f) 'g)))
#|
(down a (down b c) (up d (down e f) g))
|#
|#


(define (flip-outer-index s)
  (assert (structure? s))
  (s:generate (s:length s) (s:opposite s)
	      (lambda (i) (s:ref s i))))
#|
(flip-outer-index (up 'a (up 'b 'c) (down 'd (up 'e 'f) 'g)))
#|
(down a (up b c) (down d (up e f) g))
|#
|#

(define (typical-object s)
  (if (structure? s)
      (s:generate (s:length s) (s:same s)
		  (lambda (i)
		    (typical-object (s:ref s i))))
      (generate-uninterned-symbol 'x)))
#|
(typical-object (up 't (up 'u 'v) (down 'r 's) (up 'v1 'v2)))
#|
(up x328 (up x329 x330) (down x331 x332) (up x333 x334))
|#
|#

(define (structure->access-chains struct)
  (let lp ((struct struct) (chain '()))
    (if (structure? struct)
	(s:generate (s:length struct) (s:same struct)
		    (lambda (i)
		      (lp (s:ref struct i) (cons i chain))))
	(reverse chain))))
#|
(structure->access-chains (up 't (up 'u 'v) (down 'r 's) (up 'v1 'v2)))
#|
(up (0) (up (1 0) (1 1)) (down (2 0) (2 1)) (up (3 0) (3 1)))
|#
|#


(define (structure->prototype name struct)
  (s:map/r (lambda (chain)
	     (string->symbol
	      (string-append (symbol->string name)
	        (apply string-append
		       (map (lambda (el)
			      (string-append ":"
					     (number->string el)))
			    chain)))))
	   (structure->access-chains struct)))
#|
(structure->prototype 'foo (up 't (up 'u 'v) (down 'r 's) (up 'v1 'v2)))
#|
(up foo:0 (up foo:1:0 foo:1:1) (down foo:2:0 foo:2:1) (up foo:3:0 foo:3:1))
|#
|#


#|
;;; In src/mechanics/canonical.scm

(define (linear-function->multiplier F argument)
  ((derivative F) argument))

(define (compatible-zero s)
  (linear-function->multiplier (lambda (x) 0) s))
|#

(define (compatible-zero s)
  (if (structure? s)
      (flip-indices (s:zero-like s))
      (g:zero-like s)))
#|
(compatible-zero (up 't (up 'u 'v) (down 'r 's) (up 'v1 'v2)))
#|
(down 0 (down 0 0) (up 0 0) (down 0 0))
|#
|#

(define (compatible-shape s)
  (if (structure? s)
      (typical-object (flip-indices s))
      (typical-object s)))
#|
(compatible-shape (up 't (up 'u 'v) (down 'r 's) (up 'v1 'v2)))
#|
(down x335 (down x336 x337) (up x338 x339) (down x340 x341))
|#
|#

(define (s:transpose-outer struct)
  (s:generate (s:length (s:ref struct 0))
	      (s:same (s:ref struct 0))
	      (lambda (i)
		(s:generate (s:length struct)
			    (s:same struct)
			    (lambda (j)
			      (s:ref (s:ref struct j)
				     i))))))
#|
;;; used only in symmetrize-Christoffel in 
;;; src/calculus/covariant-derivative.scm

(define foo
  (down (down (up 'x 'y)
	      (up 'z 'w))
	(down (up 'a 'b)
	      (up 'c 'd))))

(s:transpose-outer foo)
#|
(down (down (up x y)
	    (up a b))
      (down (up z w)
	    (up c d)))
|#
|#

;;; contract assumes multi-index cubical structures

(define (s:contract struct index1 index2)
  (let ((scripts
	 (let lp ((s struct))
	   (if (structure? s)
	       (cons (cons (s:same s) (s:length s))
		     (lp (s:ref s 0)))
	       '()))))
    (assert (every (lambda (x)
		     (or (eq? (car x) 'up) (eq? (car x) 'down)))
		   scripts))
    (assert (not (eq? (car (list-ref scripts index1))
		      (car (list-ref scripts index2)))))

    (let lp ((scripts scripts) (indices '()) (index-number 0))
      (if (not (null? scripts))
	  (cond ((or (fix:= index-number index1)
		     (fix:= index-number index2))
		 (lp (cdr scripts)
		     (append indices (list 'index))
		     (fix:+ index-number 1)))
		(else
		 (s:generate (cdar scripts) (caar scripts)
			     (lambda (i)
			       (lp (cdr scripts)
				   (append indices (list i))
				   (fix:+ index-number 1))))))
	  (g:sigma (lambda (i)
		     (ref-internal struct
				   (list-with-substituted-coord 
				     (list-with-substituted-coord
				      indices index1 i)
				     index2 i)))
		   0
		   (fix:- (s:length struct) 1))))))

;;; beginning of ultra-flatten

(define (ultra-flatten s)
  (if (structure? s)
      (apply append
	     (map ultra-flatten
		  (vector->list (s:->vector s))))
      (list s)))

#|
(ultra-flatten (up 1 2 'a (down 3 4) (up (down 'c 'd) 'e)))
;Value 16: (1 2 a 3 4 c d e)

;;; similar to s:fringe (reverse)
|#

(define (s:dimension s)
  (if (structure? s)
      (reduce fix:+
	      0
	      (map s:dimension
		   (vector->list (s:->vector s))))
      1))

#|
(s:dimension (up 1 2 'a (down 3 4) (up (down 'c 'd) 'e)))
;Value: 8
|#

(define (ultra-unflatten shape list)
  (if (structure? shape)
      (let lp ((s '()) (i 0) (list list))
	(if (fix:< i (s:length shape))
	    (lp (cons (ultra-unflatten (s:ref shape i) list) s)
		(fix:+ i 1)
		(list-tail list (s:dimension (s:ref shape i))))
	    (s:structure (s:same shape) (list->vector (reverse s)))
	    ))
      (car list)))
#|
(ultra-unflatten
      (up 'x 'x 'x (down 'x 'x) (up (down 'x 'x) 'x))
      (list 1 2 'a 3 4 'c 'd 'e))
;Value 30: #(1 2 a (*down* #(3 4)) #((*down* #(c d)) e))

(pe (- (ultra-unflatten
	(up 'x 'x 'x (down 'x 'x) (up (down 'x 'x) 'x))
	(list 1 2 'a 3 4 'c 'd 'e))
       (up 1 2 'a (down 3 4) (up (down 'c 'd) 'e))))
(up 0 0 0 (down 0 0) (up (down 0 0) 0))
|#


#|
(define vs
  (velocity-tuple
   (velocity-tuple 'vx1 'vy1)
   (velocity-tuple 'vx2 'vy2)))

(define (L1 vs)
  (let ((v1 (ref vs 0))
	(v2 (ref vs 1)))
    (+ (* 1/2 'm1 (square v1))
       (* 1/2 'm2 (square v2)))))

(pe (((expt D 2) L1) vs))
(down (down (down (down m1 0) (down 0 0)) (down (down 0 m1) (down 0 0)))
      (down (down (down 0 0) (down m2 0)) (down (down 0 0) (down 0 m2))))
|#

(define *careful-conversion* #t)

(define (s->m ls ms rs)
  (if *careful-conversion*
      (assert (numerical-quantity? (g:* ls (g:* ms rs)))
              "Innapropriate s->m" ls ms rs))
  (let ((nrows (s:dimension ls))
	(ncols (s:dimension rs)))
    (m:generate nrows ncols
		(lambda (i j)
		  (g:* (ultra-unflatten
			ls
			(vector->list
			 (v:make-basis-unit nrows i)))
		       (g:* ms
			    (ultra-unflatten
			     rs
			     (vector->list
			      (v:make-basis-unit ncols j)))))))))

#|
(pe (s->m vs (((expt D 2) L1) vs) vs))
(matrix-by-rows (list m1 0 0 0)
                (list 0 m1 0 0)
                (list 0 0 m2 0)
                (list 0 0 0 m2))
|#


(define (m->s ls m rs)
  (let ((ncols (m:num-cols m))
	(col-shape (compatible-shape ls)))
    (let ((ms
           (ultra-unflatten (compatible-shape rs)
                            (let lp ((j 0))
                              (if (fix:= j ncols)
                                  '()
                                  (let ((colj (m:nth-col m j)))
                                    (cons (ultra-unflatten col-shape
                                                           (vector->list colj))
                                          (lp (fix:+ j 1)))))))))
      (if *careful-conversion*
          (assert (numerical-quantity? (g:* ls (g:* ms rs)))
                  "Innapropriate m->s" ls ms rs))
      ms)))

#|
(pe (m->s vs (s->m vs (((expt D 2) L1) vs) vs) vs))
(down (down (down (down m1 0) (down 0 0)) (down (down 0 m1) (down 0 0)))
      (down (down (down 0 0) (down m2 0)) (down (down 0 0) (down 0 m2))))
|#

;;; In the following procedures there are extra arguments, ls and rs.
;;; ls designates the shape of the object to be multiplied on the left
;;; of the result and rs designates the shape of the object to be
;;; multiplied on the right of the result.


;;; (* (* O |a>) |b>) = (* (* (s:flip O) |b>) |a>)

(define (s:flip rs ms ls)
  (m->s ls
	(m:transpose (s->m rs ms ls))
	rs))


;;; (* <a| O |b>) = (* <b| (s:transpose O) |a>)
#|
(define (s:transpose ls ms rs)
  (m->s (compatible-shape rs)
	(m:transpose (s->m ls ms rs))
	(compatible-shape ls)))
|#

(define (s:transpose ls ms rs)
   (m->s rs
         (m:transpose (s->m ls ms rs))
         ls))

(define (s:inverse ls ms rs)		;but see s:invert...
  (m->s (compatible-shape rs)
	(m:invert
	 (s->m ls ms rs))
	(compatible-shape ls)))


;;; Any one argument function of a structure can be seen
;;; as a matrix.  This is only useful if the function 
;;; has a linear multiplier (e.g. derivative)

(define ((as-matrix F) s)
  (let ((v (F s)))
    (s->m (compatible-shape (g:* v s)) v s)))

#|
(define ((D-as-matrix F) s)
  (s->m (compatible-shape (F s)) ((D F) s) s))


(define C-general
  (literal-function 'C
    (-> (UP Real
	    (UP Real Real)
	    (DOWN Real Real))
	(UP Real
	    (UP Real Real)
	    (DOWN Real Real)))))


(define s (up 't (up 'x 'y) (down 'px 'py)))

((as-matrix (D C-general)) s)
#|
(matrix-by-rows
 (list (((partial 0) C^0) (up t (up x y) (down px py)))
       (((partial 1 0) C^0) (up t (up x y) (down px py)))
       (((partial 1 1) C^0) (up t (up x y) (down px py)))
       (((partial 2 0) C^0) (up t (up x y) (down px py)))
       (((partial 2 1) C^0) (up t (up x y) (down px py))))
 (list (((partial 0) C^1^0) (up t (up x y) (down px py)))
       (((partial 1 0) C^1^0) (up t (up x y) (down px py)))
       (((partial 1 1) C^1^0) (up t (up x y) (down px py)))
       (((partial 2 0) C^1^0) (up t (up x y) (down px py)))
       (((partial 2 1) C^1^0) (up t (up x y) (down px py))))
 (list (((partial 0) C^1^1) (up t (up x y) (down px py)))
       (((partial 1 0) C^1^1) (up t (up x y) (down px py)))
       (((partial 1 1) C^1^1) (up t (up x y) (down px py)))
       (((partial 2 0) C^1^1) (up t (up x y) (down px py)))
       (((partial 2 1) C^1^1) (up t (up x y) (down px py))))
 (list (((partial 0) C^2_0) (up t (up x y) (down px py)))
       (((partial 1 0) C^2_0) (up t (up x y) (down px py)))
       (((partial 1 1) C^2_0) (up t (up x y) (down px py)))
       (((partial 2 0) C^2_0) (up t (up x y) (down px py)))
       (((partial 2 1) C^2_0) (up t (up x y) (down px py))))
 (list (((partial 0) C^2_1) (up t (up x y) (down px py)))
       (((partial 1 0) C^2_1) (up t (up x y) (down px py)))
       (((partial 1 1) C^2_1) (up t (up x y) (down px py)))
       (((partial 2 0) C^2_1) (up t (up x y) (down px py)))
       (((partial 2 1) C^2_1) (up t (up x y) (down px py)))))
|#

((- (D-as-matrix C-general)
    (as-matrix (D C-general)))
 s)
#|
(matrix-by-rows (list 0 0 0 0 0)
		(list 0 0 0 0 0)
		(list 0 0 0 0 0)
		(list 0 0 0 0 0)
		(list 0 0 0 0 0))
|#
|#