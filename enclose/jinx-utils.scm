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

;;; -*- Scheme -*-

(declare (usual-integrations))

;;; $Header: utils.scm,v 1.13 90/08/23 02:07:44 GMT jinx Exp $

;;;; Control utilities

;;; Allowing callees to default correctly.

(define make-default-object
  ;; Cacheing kludge
  (let ((result '()))
    (named-lambda (make-default-object)
      (if (null? result)
	  (set! result
		(cons
		 (unmap-reference-trap (make-unassigned-reference-trap))
		 '())))
      (car result))))

;;; Keyword parameter passing.

;; This should include an option that disables "Unknown option" errors
;; and passes left over parameters along.

(define (get-options all-options option-names recvr)
  (define (pass paired-options)
    (let loop ((names (reverse option-names))
	       (parameters '()))
      (cond ((null? names)
	     (apply recvr parameters))
	    ((assq (car names) paired-options)
	     =>
	     (lambda (pair)
	       (loop (cdr names)
		     (cons (cdr pair) parameters))))
	    (else
	     (loop (cdr names)
		   (cons (make-default-object) parameters))))))	  

  (let loop ((options all-options)
	     (acc '()))
    (cond ((null? options)
	   (pass acc))
	  ((not (memq (car options) option-names))
	   (error "get-options: Unknown option" (car options)))
	  ((null? (cdr options))
	   (error "get-options: No value" (car options)))
	  (else
	   (loop (cddr options)
		 (cons (cons (car options)
			     (cadr options))
		       acc))))))

;;;; Higher order utilities

(define make-reducer
  (let ()
    (define (construct-reducer binop direction min-args max-args null-value)
      (let ((reducer-maker
	     (if (eq? direction 'RIGHT)
		 (named-lambda (make-right-reducer handle-last)
		   (define (reducer next rest)
		     (if (null? rest)
			 (handle-last next)
			 (binop next
				(reducer (car rest) (cdr rest)))))
		   reducer)
		 (named-lambda (make-left-reducer handle-last)
		   (lambda (next rest)
		     (let loop ((accum (handle-last next))
				(rest rest))
		       (if (null? rest)
			   accum
			   (loop (binop accum (car rest))
				 (cdr rest))))))))
	    (handle1
	     (cond ((or (zero? max-args)
			(< max-args min-args))
		    (lambda (x) x))
		   ((eq? direction 'RIGHT)
		    (lambda (x)
		      (binop x null-value)))
		   (else
		    (lambda (x)
		      (binop null-value x))))))
	(case min-args
	  ((2)
	   (let ((process (reducer-maker handle1)))
	     (lambda (all)
	       (if (or (null? all) (null? (cdr all)))
		   (error "reducer: Too few arguments" all)
		   (process (car all) (cdr all))))))
	  ((1)
	   (if (= max-args 2)
	       (let ((process (reducer-maker handle1)))
		 (lambda (all)
		   (if (null? all)
		       (error "reducer: Too few arguments" all)
		       (process (car all) (cdr all)))))
	       (let ((process (reducer-maker (lambda (x) x))))
		 (lambda (all)
		   (cond ((null? all)
			  (error "reducer: Too few arguments" all))
			 ((null? (cdr all))
			  (handle1 (car all)))
			 (else
			  (process (car all) (cdr all))))))))
	  ((0)
	   (if (= max-args 2)
	       (let ((process (reducer-maker handle1)))
		 (lambda (all)
		   (if (null? all)
		       null-value
		       (process (car all) (cdr all)))))
	       (let ((process (reducer-maker (lambda (x) x))))
		 (lambda (all)
		   (cond ((null? all)
			  null-value)
			 ((null? (cdr all))
			  (handle1 (car all)))
			 (else
			  (process (car all) (cdr all))))))))
	  (else
	   (error "make-tail-collector: Inconsistency" min-args)))))

    (named-lambda (make-reducer binop . options)
      (get-options
       options
       '(DIRECTION MIN-ARGS MAX-ARGS-USING-NULL-VALUE NULL-VALUE)
       (lambda (direction min-args max-args null-value)
	 (let* ((min-args
		 (if (default-object? min-args)
		     2
		     min-args))
		(max-args
		 (cond ((not (default-object? max-args))
			max-args)
		       ((default-object? null-value)
			0)
		       (else
			2))))
	   ;; Paranoia check.
	   (cond ((or (not (integer? min-args))
		      (not (<= 0 min-args 2)))
		  (error "make-reducer: Bad min-args" min-args))
		 ((or (not (integer? max-args))
		      (not (<= 0 max-args 2)))
		  (error "make-reducer: Bad max-args-using-null-value"
			 max-args))
		 ((default-object? null-value)
		  (if (or (= min-args 0)
			  (>= max-args min-args))
		      (error "make-reducer: required null-value not supplied"
			     `((MIN-ARGS ,min-args)
			       (MAX-ARGS-USING-NULL-VALUE ,max-args)))))
		 ((< max-args min-args)
		  (error "make-reducer: null-value meaningless"
			 `((NULL-VALUE ,null-value)
			   (MIN-ARGS ,min-args)
			   (MAX-ARGS-USING-NULL-VALUE ,max-args)))))

	   (construct-reducer binop
			      (if (default-object? direction)
				  'LEFT
				  direction)
			      min-args
			      max-args
			      (if (default-object? null-value)
				  (make-default-object)
				  null-value))))))))

;;;; List utilities

(define (split-list list predicate recvr)
  (let split ((list list)
	      (recvr recvr))
    (if (not (pair? list))
	(recvr '() '())
	(split (cdr list)
	       (lambda (win lose)
		 (if (predicate (car list))
		     (recvr (cons (car list) win)
			    lose)
		     (recvr win
			    (cons (car list) lose))))))))

(define (find-infimum list predicate)
  (if (null? list)
      (error "find-infimum: empty list" list))
  (let loop ((current (car list))
	     (left (cdr list)))
    (cond ((null? left)
	   current)
	  ((predicate (car left) current)
	   (loop (car left) (cdr left)))
	  (else
	   (loop current (cdr left))))))

(define (subst new old where)
  (cond ((eq? where old)
	 new)
	((not (pair? where))
	 where)
	(else
	 (cons (subst new old (car where))
	       (subst new old (cdr where))))))

(define (delq-once element list)
  (cond ((null? list)
	 '())
	((eq? (car list) element)
	 (cdr list))
	(else
	 (cons (car list)
	       (delq-once element (cdr list))))))

;;;; Mapping and reducing

;; Important: All of these are iterative, so they won't run out of stack!

(define (map&reduce procedure combiner null-value list1 #!optional list2 . lists)
  ;; (reduce combiner null-value (map procedure list1 list2 . lists))
  (cond ((default-object? list2)
	 (let loop ((result null-value)
		    (l list1))
	   (if (null? l)
	       result
	       (loop (combiner (procedure (car l))
			       result)
		     (cdr l)))))
	((null? lists)
	 (let loop ((result null-value)
		    (l1 list1)
		    (l2 list2))
	   (if (or (null? l1) (null? l2))
	       result
	       (loop (combiner (procedure (car l1) (car l2))
			       result)
		     (cdr l1)
		     (cdr l2)))))
	(else
	 (let loop ((result null-value)
		    (l (cons* list1 list2 lists)))
	   (if (there-exists? l null?)
	       result
	       (loop (combiner (apply procedure (map car l))
			       result)
		     (map cdr l)))))))

(define (%append x y)
  (if (null? x)
      y
      (%reverse! (%reverse x '()) y)))
  
(define (%reverse! l #!optional tail)
  (let loop ((current l)
	     (new-cdr (if (default-object? tail)
			  '()
			  tail)))
    (if (pair? current)
	(let ((next (cdr current)))
	  (set-cdr! current new-cdr)
	  (loop next current))
	(begin
	  (if (not (null? current))
	      (error "%REVERSE!: Argument not a list" l))
	  new-cdr))))

(define (%reverse ol #!optional tail)
  (let loop ((l ol)
	     (accum (if (default-object? tail)
			'()
			tail)))
    (cond ((pair? l)
	   (loop (cdr l)
		 (cons (car l) accum)))
	  ((null? l)
	   accum)
	  (else
	   (error "%REVERSE: Argument not a list" ol)))))  

(define (%map f ol1 #| #!optional ol2 . rest |#)
  ;; Important: The circular list hack for multi-argument
  ;; map does not work here.
  (cond ((default-object? l2)
	 (%map-1 f ol1))
	((null? rest)
	 (%map-2 f ol1 ol2))
	(else
	 (let outer ((result '())
		     (ls (reverse (%map-1 reverse (cons* ol1 ol2 rest)))))
	   (cond ((pair? (car ls))
		  (let inner ((args (list (caar ls)))
			      (next (list (cdar ls)))
			      (rest (cdr ls)))
		    (cond ((null? rest)
			   (outer (cons (apply f args) result)
				  (reverse! next)))
			  ((not (pair? (car rest)))
			   (error "%map: Arguments have different lengths"
				  (cons* ol1 ol2 rest)))
			  (else
			   (inner (cons (caar rest) args)
				  (cons (cdar rest) next)
				  (cdr rest))))))
		 ((there-exists? ls (lambda (x) (not (null? x))))
		  (error "%map:Arguments have different lengths"))
		 (else
		  result))))))

(define-integrable (%map-1 f ol)
  (let loop ((result '()) (l1 (reverse ol)))
    (cond ((pair? l1)
	   (loop (cons (f (car l1)) result)
		 (cdr l1)))
	  ((null? l1)
	   result)
	  (else
	   (error "%map: Argument not a list" ol)))))      

(define-integrable (%map-2 f ol1 ol2)
  (let loop ((result '())
	     (l1 (reverse ol1))
	     (l2 (reverse ol2)))
    (cond ((and (pair? l1) (pair? l2))
	   (loop (cons (f (car l1) (car l2)) result)
		 (cdr l1)
		 (cdr l2)))
	  ((and (null? l1) (null? l2))
	   result)
	  (else
	   (error "%map: Arguments have different lengths"
		  ol1 ol2)))))

;;;; Set utilities

(define-integrable (eq-set/make-empty)
  '())

(define-integrable (eq-set/empty? set)
  (null? set))

(define-integrable (eq-set/member? element set)
  (memq element set))

(define-integrable (eq-set/adjoin element set)
  (if (eq-set/member? element set)
      set
      (cons element set)))

(define (eq-set/remove element set)
  (if (not (eq-set/member? element set))
      set
      (delq element set)))

;; Important: This will return set2 itself (rather than a copy) if the
;; union is set2.  Thus eq? can be used on the return value to
;; determine whether the set has grown.

(define (eq-set/union set1 set2)
  (define (loop set new-elements)
    (if (null? new-elements)
	set
	(loop (eq-set/adjoin (car new-elements) set)
	      (cdr new-elements))))

  ;; If set2 is smaller than set1, the union is guaranteed not to be set2.
  (if (< (length set2) (length set1))
      (loop set1 set2)
      (loop set2 set1)))

(define (eq-set/intersection set1 set2)
  (define (examine set1 set2)
    (let process ((set #| (reverse set1) |# set1)
		  (result (eq-set/make-empty)))
      (if (null? set)
	  result
	  (process (cdr set)
		   (if (eq-set/member? (car set) set2)
		       (cons (car set) result)
		       result)))))

  (if (< (length set2) (length set1))
      (examine set2 set1)
      (examine set1 set2)))

(define (eq-set/difference set1 set2)
  (if (null? set2)
      set1
      (let process ((set set1) (result (eq-set/make-empty)))
	(cond ((null? set)
	       result)
	      ((eq-set/member? (car set) set2)
	       (process (cdr set) result))
	      (else
	       (process (cdr set)
			(cons (car set) result)))))))

(define (eq-set/subset? set1 set2)
  (or (eq-set/empty? set1)
      (and (eq-set/member? (car set1) set2)
	   (eq-set/subset? (cdr set1) set2))))

(define (eq-set/equal? set1 set2)
  (or (eq? set1 set2)
      (and (eq-set/subset? set1 set2)
	   (eq-set/subset? set2 set1))))

;;;; Multi set utilities

(define-integrable (multi-set/empty)
  '())

(define-integrable (multi-set/adjoin element set)
  (cons element set))

(define-integrable (multi-set/empty? set)
  (null? set))

(define-integrable (multi-set/first set)
  (car set))

(define-integrable (multi-set/rest set)
  (cdr set))

(define-integrable (multi-set/remove element set)
  (delq-once element set))

(define-integrable (multi-set/element? element set)
  (memq element set))

(define-integrable (multi-set/union set1 set2)
  (%reverse set1 set2))

(define (multi-set/intersection set1 set2)
  (define (process set1 set2 result)
    (cond ((multi-set/empty? set1)
	   result)
	  ((not (multi-set/element? (multi-set/first set1) set2))
	   (process (multi-set/rest set1) set2 result))
	  (else
	   (process (multi-set/rest set1)
		    (multi-set/remove (multi-set/first set1)
				      set2)
		    (multi-set/adjoin (multi-set/first set1)
				      result)))))

  (if (< (length set2) (length set1))
      (process set2 set1 (multi-set/empty))
      (process set1 set2 (multi-set/empty))))

(define (multi-set/difference set1 set2)
  (define (process set1 set2 result)
    (cond ((multi-set/empty? set1)
	   result)
	  ((multi-set/element? (multi-set/first set1) set2)
	   (process (multi-set/rest set1)
		    (multi-set/remove (multi-set/first set1) set2)
		    result))
	  (else
	   (process (multi-set/rest set1)
		    set2
		    (multi-set/adjoin (multi-set/first set1)
				      result)))))
  (process set1 set2 (multi-set/empty)))

;;;; Random utilities: association tables (eq? based)

(let-syntax ((primitive
	      (lambda (name)
		`',(make-primitive-procedure name))))
  (define-integrable string-hash (primitive string-hash))
  (define-integrable symbol-print-name (primitive system-pair-car)))

(define (find-next-prime number)
  (let loop ((primes 
	      '(1009 2003 4001 8009 16001 32003 64007
		     128021 256019 512009 1024021)))
    (cond ((null? primes)
	   (error "find-next-prime: number too large" number))
	  ((< number (car primes))
	   (car primes))
	  (else (loop (cdr primes))))))

(define-integrable (size->table-size size)
 (find-next-prime (quotient (* 5 size) 4)))

(define-integrable (%table/make size)
  (make-vector size '()))

(define (table/make size)
  (%table/make (size->table-size size)))

(define (table/index table object)
  (modulo
   (cond ((exact-integer? object)
	  object)
	 ((symbol? object)
	  (string-hash (symbol-print-name object)))
	 ((string? object)
	  (string-hash object))
	 (else
	  (object-hash object)))
   (vector-length table)))

(define (table/association table key)
  (let ((bucket (vector-ref table (table/index table key))))
    (and (not (null? bucket))
	 (let ((place (assq key bucket)))
	   (and (pair? place)
		(cdr place))))))

(define (table/associate! table key value)
  (let* ((index (table/index table key))
	 (bucket (vector-ref table index))
	 (place (assq key bucket)))
    (if (not place)
	(begin
	  (vector-set! table index (cons (cons key value) bucket))
	  true)
	(begin
	  (set-cdr! place value)
	  false))))

;;;; Growing association tables

(define-structure (table+ (conc-name table+/)
			  (constructor %table+/make))
  (size false read-only false)
  (limit false read-only false)
  (entries false read-only false)
  (table false read-only false))

(define-integrable (table+-size->limit size)
  (quotient (* 4 size) 5))

(define (table+/make #!optional size)
  (let* ((size
	  ;; Bug in sf 4.8
	  (let ((input-size
		 (if (default-object? size)
		     100
		     size)))
	    (size->table-size input-size)))
	 (table (%table/make size)))
    (%table+/make size
		  (table+-size->limit size)
		  0
		  table)))

(define-integrable (table+/association table+ key)
  (table/association (table+/table table+) key))

(define (table+/associate! table+ key value)
  (if (table/associate! (table+/table table+) key value)
      (let ((entries (1+ (table+/entries table+))))
	(set-table+/entries! table+ entries)
	(if (> entries (table+/limit table+))
	    (table+/grow! table+)))))

(define (table+/grow! table+)
  (let* ((next-size (find-next-prime (table+/size table+)))
	 (new-table (%table/make next-size)))
    (for-each-vector-element
     (table+/table table+)
     (lambda (bucket)
       (for-each (lambda (pair)
		   (table/associate! new-table (car pair) (cdr pair)))
		 bucket)))
    (set-table+/limit! table+ (table+-size->limit next-size))
    (set-table+/table! table+ new-table)
    (set-table+/size! table+ next-size)))

;;;; I/O utilities

(define (warning . arguments)
  (apply warn arguments)
  (warning/stop-hook))

(define (warning/stop-hook/bkpt)
  (bkpt "Breakpoint at warning"))

(define (warning/stop-hook/default)
  unspecific)

(define warning/stop-hook
  warning/stop-hook/default)

(define message-tag
  (list '*MESSAGE-NOISE*))

(define (message/noise noise)
  (cons message-tag noise))

(define (message/noise? noise)
  (and (pair? noise)
       (eq? (car noise) message-tag)))

(define (message/pluralize string number #!optional suffix)
  (let ((result
	 (if (= number 1)
	     string
	     (string-append string "s"))))
    (message/noise
     (if (default-object? suffix)
	 result
	 (string-append result suffix)))))

(define (message string . values)
  (newline)
  (write-string string)
  (for-each (lambda (value)
	      (write-char #\Space)
	      (if (message/noise? value)
		  (write-string (cdr value))
		  (write value)))
	    values))

;;;; PP does a horrible job on lists.

(define (write-list all-l #!optional tab-position tagged?)
  (let ((port (current-output-port)))

    (let ((write-list/write-string
	   (output-port/operation/write-string port))
	  (write-list/write-char
	   (output-port/operation/write-char port)))

      (let ((port port)
	    (initial-tab-position (if (default-object? tab-position)
				      0
				      tab-position))
	    (size (output-port/x-size port)))
	(let* ((tab-position
		(if (default-object? tagged?)
		    initial-tab-position
		    (+ (1+ initial-tab-position)
		       (string-length (symbol->string (car all-l))))))
	       (prefix (make-string (1+ tab-position) #\Space)))

	  (define-integrable (write-string string)
	    (write-list/write-string port string))

	  (define-integrable (write-char char)
	    (write-list/write-char port char))

	  (define-integrable (newline)
	    (write-list/write-char port #\newline))

	  (newline)
	  (write-string (make-string initial-tab-position #\Space))
	  (write-char #\()
	  (let loop ((l all-l)
		     (start? true)
		     (left (- size (1+ initial-tab-position))))
	    (cond ((not (null? l))
		   (let* ((next (write-to-string (car l)))
			  (new-size
			   (if (not (null? (cdr l)))
			       (1+ (string-length next))
			       (+ (string-length next) 2))))
		     (cond ((<= new-size left)
			    (if start?
				(begin
				  (write-string next)
				  (loop (cdr l)
					false
					(- left (-1+ new-size))))
				(begin
				  (write-char #\Space)
				  (write-string next)
				  (loop (cdr l)
					false
					(- left new-size)))))
			   (start?
			    (write-string next)
			    (newline)
			    (write-string prefix)
			    (loop (cdr l)
				  true
				  left))
			   (else
			    (newline)
			    (write-string prefix)
			    (write-string next)
			    (loop (cdr l)
				  false
				  (- size (+ tab-position new-size)))))))
		  ((not (positive? left))
		   (newline)
		   (write-string prefix))))
	  (write-char #\))
	  (output-port/flush-output port)
	  unspecific)))))

;;;; Structure utilities

(define *structure-unparse-level* 1)

(define (structure/unparse state node description
			   #!optional show-false?)
  (let ((show-false?
	 (if (default-object? show-false?)
	     true
	     show-false?)))
    (unparse-string state "#[")
    (unparse-string state (car description))
    (unparse-char state #\Space)
    (unparse-object state (object-hash node))
    (if (not (zero? *structure-unparse-level*))
	(fluid-let ((*structure-unparse-level*
		     (-1+ *structure-unparse-level*)))
	  (for-each
	   (lambda (field)
	     (let ((value ((cadr field) node)))
	       (and (or show-false? value)
		    (begin
		      (unparse-char state #\Space)
		      (unparse-string state (car field))
		      (unparse-string state ": ")
		      (unparse-object state value)))))
	   (cdr description))))
    (unparse-string state "]")))

(define (@ object . selectors)
  (let loop ((object object)
	     (selectors selectors))
    (if (null? selectors)
	object
	(loop ((car selectors) object)
	      (cdr selectors)))))

(define (write-structures . structures)
  (for-each (lambda (obj)
	      (write-line
	       (cond ((not (number? obj))
		      obj)
		     ((valid-hash-number? obj)
		      (object-unhash obj))
		     (else
		      `(invalid-hash-number ,obj)))))
	    structures))

;;;; Random Scheme utilities

(define (errset thunk success failure)
  ((call-with-current-continuation
    (lambda (abort)
      (let ((result
	     (bind-condition-handler condition-type:error
		 (lambda (condition)
		   condition		; ignored
		   (abort failure))
	       thunk)))
	;; This could invoke abort,
	;; but there is no reason.
	(lambda ()
	  (success result)))))))
#|
;; Adapted from generate-uninterned-symbol in the runtime system
;; to allow string prefixes.

(define new-uninterned-symbol
  (let ((name-counter 0)
	(name-prefix "G"))
    (lambda (#!optional argument)
      (if (not (default-object? argument))
	  (cond ((symbol? argument)
		 (set! name-prefix (symbol->string argument)))
		((string? argument)
		 (set! name-prefix argument))
		((exact-nonnegative-integer? argument)
		 (set! name-counter argument))
		(else
		 (error "new-uninterned-symbol: Bad argument" argument))))
      (string->uninterned-symbol
       (string-append name-prefix
		      (number->string
		       (let ((result name-counter))
			 (set! name-counter (1+ name-counter))
			 result)))))))
|#

;; The version in the runtime system now allows strings as well.

(define-integrable new-uninterned-symbol
  generate-uninterned-symbol)