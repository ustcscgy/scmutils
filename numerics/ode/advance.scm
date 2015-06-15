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

;;; The advance-generator is used with 1-step adaptive integrators to
;;;  to advance a state by a given increment in the independent
;;;  variable, in the face of variable-stepsize advancers.

(declare (usual-integrations))

#| For examples

((advance-generator
  ((quality-control rk4 4)		;integration method
   (lambda (v) v)			;x' = x
   .0001))				;error tolerated
 #(1.0)					;initial state (at t = t0)
 1.0					;proceed to t = t0 + 1
 0.1					;first step no larger than .1
 0.5					;no step larger than .5
 (lambda (ns dt h cont)
   (pp ns)
   (cont))
 (lambda (ns dt sdt)
   ;; assert ns = #(2.718...)
   ;; assert dt = 1.000...+-
   (list ns dt sdt)))

((advance-generator
  (bulirsch-stoer-lisptran		;integrator
   (lambda (vin vout)			;x'= x
     (vector-set! vout 0
		  (vector-ref vin 0)))
   1					;1 dimension
   .0001))				;error tolerated
 #(1.0)
 1.0
 0.1
 0.5
 (lambda (ns dt h cont)
   (pp (list dt ns))
   (cont))
 (lambda (ns dt sdt)
   ;; assert ns = #(2.718...)
   ;; assert dt = 1.000...+-
   (list ns dt sdt)))
|#

#|
(declare (usual-integrations + - * / < > = zero? positive? negative?
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
	  #|(atan flonum-atan 1)
	  (atan flonum-atan2 2)|#))
|#

(define (advance-generator advancer)
  (define (advance start-state step-required h-suggested max-h continue done)
    ;; done = (lambda (end-state step-achieved h-suggested) ... )
    ;; continue = (lambda (state step-achieved h-taken next)
                    ;; next = (lambda () ...))
    (let lp ((state start-state)
	     (step-achieved 0.0)
	     (h (min-step-size step-required h-suggested max-h)))
      (if advance-wallp?
	  (pp `(advance: ,step-achieved ,state)))
      (continue state step-achieved h
	(lambda ()
          (advancer state h
            (lambda (new-state step-obtained h-suggested)
              (let ((ndt (+ step-achieved step-obtained)))
                (if (close-enuf? step-required ndt
                                 *independent-variable-tolerance*)
                    (done new-state ndt h-suggested)
                    (lp new-state
                        ndt
                        (min-step-size (- step-required ndt)
                                       h-suggested
                                       max-h))))))))))
  advance)

(define advance-wallp? false)

(define-integrable (first-with-sign-of-second x y)
  (if (> y 0.0)
      (abs x)
      (- (abs x))))

(define (min-step-size step-required h-suggested max-h)
  (let ((h-allowed (min (abs h-suggested) (abs max-h))))
    (if (< (abs step-required) h-allowed)
	step-required
	(first-with-sign-of-second h-allowed step-required))))


(define *independent-variable-tolerance* 1.0e-100)

#| For making a stream of states

((stream-of-states
  (advance-generator
   ((quality-control rk4 4)		;integration method
    (lambda (v) v)			;x' = x
    .0001)))				;error tolerated
 #(1.0)					;initial state (at t = t0)
 1.0					;proceed to t = t0 + 1
 0.1					;first step no larger than .1
 0.5					;no step larger than .5
 )
|#

(define (stream-of-states advance)
  (lambda (state step-required h-suggested max-h)
    (advance state step-required h-suggested max-h
	     (lambda (state step-achieved h cont)
	       (cons-stream state (cont)))
	     (lambda (state step-achieved h-suggested)
	       (cons-stream state the-empty-stream)))))


;;; Utilities for ODE integrators.

(define (vector-fixed-point-with-failure update start measure succeed fail)
  ;; update  = (lambda (x cont)
                 ;; cont = (lambda (nx fx) ...)
                 ;; ...)
  ;; succeed = (lambda (nx fx count) ...)
  ;; fail    = (lambda () ...)
  (let improve ((last-value start) (count 1))
    (if (> (maxnorm last-value) *vector-fixed-point-ridiculously-large*)
	(fail last-value last-value)
	(update last-value
		(lambda (value fvalue)
		  (let ((d (measure value last-value)))
		    (if *fixed-point-wallpaper*
			(write-line
			 `(vector-fixed-point ,count d ,d ,start ,value)))
		    (if (< d 2.0)
			(succeed value fvalue count)
			(if (fix:> count *vector-fixed-point-iteration-loss*)
			    (fail value fvalue)
			    (improve value (1+ count))))))))))



(define *vector-fixed-point-iteration-loss* 20)

(define *vector-fixed-point-ridiculously-large* 1.0e50)

(define *fixed-point-wallpaper* false)

;;; Error measures

(define (vector-error vh vf)
  (make-initialized-vector (vector-length vh)
    (lambda (i)
      (let ((yh (vector-ref vh i))
	    (yf (vector-ref vf i)))
	(if (and (zero? yh) (zero? yf))
	    0.0
	    (/ (- yh yf)
	       (+ (magnitude yh) (magnitude yf) 2.0)
	       0.5))))))

(define (maxnorm-relabs tolerance)
  (define (vector-error vh vf)
    (let ((n (vector-length vh)))
      (let lp ((i 0) (ans 0.0))
	(if (fix:= i n)
	    (/ (* 2.0 ans) tolerance)
	    (let ((yh (vector-ref vh i))
		  (yf (vector-ref vf i)))
	      (if (and (zero? yh) (zero? yf))
		  (lp (fix:1+ i) ans)
		  (lp (fix:1+ i)
		      (max (/ (magnitude (- yh yf))
			      (+ (magnitude yh) (magnitude yf) 2.0))
			   ans))))))))
  vector-error)

(define (maxnorm-relabs-each tolerances)
  (define (vector-error vh vf)
    (let ((n (vector-length vh)))
      (let lp ((i 0) (ans 0.0))
	(if (fix:= i n)
	    (* 2.0 ans)
	    (let ((yh (vector-ref vh i))
		  (yf (vector-ref vf i)))
	      (if (and (zero? yh) (zero? yf))
		  (lp (fix:1+ i) ans)
		  (lp (fix:1+ i)
		      (max (/ (magnitude (- yh yf))
			      (+ (magnitude yh) (magnitude yf) 2.0)
			      (vector-ref tolerances i))
			   ans))))))))
  vector-error)

(define (maxnorm-relabs-scaled tolerance scales)
  (define (vector-error vh vf)
    (let ((n (vector-length vh)))
      (let lp ((i 0) (ans 0.0))
	(if (fix:= i n)
	    (/ (* 2.0 ans) tolerance)
	    (let ((yh (vector-ref vh i))
		  (yf (vector-ref vf i)))
	      (if (and (zero? yh) (zero? yf))
		  (lp (fix:1+ i) ans)
		  (lp (fix:1+ i)
		      (max (/ (magnitude (- yh yf))
			      (+ (magnitude yh) (magnitude yf)
				 (* 2.0 (vector-ref scales i))))
			   ans))))))))
  vector-error)

(define (maxnorm-relabs-each-scaled tolerances scales)
  (define (vector-error vh vf)
    (let ((n (vector-length vh)))
      (let lp ((i 0) (ans 0.0))
	(if (fix:= i n)
	    (* 2.0 ans)
	    (let ((yh (vector-ref vh i))
		  (yf (vector-ref vf i)))
	      (if (and (zero? yh) (zero? yf))
		  (lp (fix:1+ i) ans)
		  (lp (fix:1+ i)
		      (max (/ (magnitude (- yh yf))
			      (+ (magnitude yh) (magnitude yf)
				 (* 2.0 (vector-ref scales i)))
			      (vector-ref tolerances i))
			   ans))))))))
  vector-error)

(define (parse-error-measure tolerance-specification #!optional multiplier)
  (if (default-object? multiplier) (set! multiplier 1.0))
  (cond ((number? tolerance-specification) ;uniform relative error -- scale = 1
	 (maxnorm-relabs (* multiplier tolerance-specification)))
	((vector? tolerance-specification) ;tolerance for each coord -- scale = 1
	 (maxnorm-relabs-each
	  (scalar*vector multiplier tolerance-specification)))
	((and (pair? tolerance-specification)
	      (number? (car tolerance-specification))
	      (vector? (cdr tolerance-specification)))
	 (maxnorm-relabs-scaled (* multiplier (car tolerance-specification))
				(cdr tolerance-specification)))
	((and (pair? tolerance-specification)
	      (vector? (car tolerance-specification))
	      (vector? (cdr tolerance-specification)))
	 (maxnorm-relabs-scaled-each
	  (scalar*vector multiplier (car tolerance-specification))
	  (cdr tolerance-specification)))
	((procedure? tolerance-specification) ;arbitrary user-supplied procedure
	 tolerance-specification)
	(else
	 (error "Unknown tolerance specification -- PARSE-ERROR-MEASURE"))))



;;; For integrators that need a partial Jacobian, we need to be able
;;; to clip and pad vectors.

;;; A dimension is either single, in which case no clipping and
;;; padding is necessary, or it is a triplet of 
;;;   (state-size start-index end-index)
;;; such that the Jacobian only applies to the coordinates from
;;; start-index (inclusive) to end-index (exclusive).

(define (vector-clipper dimension)
  (cond ((number? dimension)
	 (lambda (x) x))
	((list? dimension)
	 (let ((start (cadr dimension)) (end (caddr dimension)))
	   (lambda (v) (subvector v start end))))
	(else (error "Bad dimension -- VECTOR-CLIPPER" dimension))))

(define (vector-padder dimension)
  (cond ((number? dimension)
	 (lambda (x) x))
	((list? dimension)
	 (let* ((state-size (car dimension))
		(start (cadr dimension))
		(end (caddr dimension))
		(j-size (fix:- end start)))
	   (lambda (v)
	     (make-initialized-vector state-size
	       (lambda (i)
		 (if (and (or (fix:< start i)
			      (fix:= start i))
			  (fix:< i end))
		     (vector-ref v (fix:- i start))
		     0.0))))))
	(else (error "Bad dimension -- VECTOR-PADDER" dimension))))

(define (J-dimension dimension)
  (cond ((number? dimension) dimension)
	((list? dimension) (fix:- (caddr dimension) (cadr dimension)))
	(else (error "Bad dimension -- J-DIMENSION" dimension))))




(define integrator-table (make-table "integrator-table" assq))

(define (add-integrator! name maker-procedure needs)
  (adjoin-to-list! name integrator-table 'integrators)
  (put! integrator-table maker-procedure name 'maker)
  (put! integrator-table needs name 'needs))
	
