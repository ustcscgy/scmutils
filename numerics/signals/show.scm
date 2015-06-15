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

;;;; This makes the "oscilliscope"
(declare (usual-integrations))

;;; Needs space between traces, and way to undefault the scaling.

(define *the-scope*)

(define *ntraces* 5)
(define *trace-width-pixels* 512)
(define *trace-height-pixels* 100)
(define *trace-separation-pixels* 5)
(define *trace-shrink* .9)

(define *scope-left-edge* 5)
(define *scope-top-edge* 105)



(define (make-scope #!optional ntraces)
  (if (not (default-object? ntraces)) (set! *ntraces* ntraces))
  (set! *the-scope*
	(frame 0 +1 (- *ntraces*) 0
	       *trace-width-pixels*
	       (+ (* *trace-height-pixels* *ntraces*)
		  (* *trace-separation-pixels* (- *ntraces* 1)))
	       *scope-left-edge* *scope-top-edge*))
  (list 'new-scope *ntraces*))

(define (flush-scope)
  (graphics-close *the-scope*))


(define (plot-trace trace signal-function #!optional dots?)
  (if (default-object? dots?) (set! dots? #f))
  (let* ((span (sigfun:span signal-function))
	 (minx (exact->inexact (sigfun:min span)))
	 (maxx (exact->inexact (sigfun:max span)))
	 (xspan (- maxx minx))
	 (dx (/ xspan *nsamples*))
	 (xdata (generate-list *nsamples* (lambda (i) (+ minx (* i dx)))))
	 (ydata
	  (map (lambda (x)
		 (real-part ((sigfun:procedure signal-function) x)))
	       xdata))
	 (miny (apply min ydata))
	 (maxy (apply max ydata))
	 (yspan (/ (- maxy miny) *trace-shrink*))
	 (yoffset (/ (- 1 *trace-shrink*) 2)))
    (graphics-set-clip-rectangle *the-scope* 0 (- trace) 1 (- (- trace 1)))
    (graphics-clear *the-scope*)
    (graphics-bind-line-style *the-scope* 1 ; dashs
     (lambda ()
       (graphics-draw-line *the-scope* 0 (- trace) 1 (- trace))
       (graphics-draw-line *the-scope* 0 (- (- trace 1)) 1 (- (- trace 1)))))
    (cond ((= yspan 0)
	   (graphics-bind-line-style *the-scope* 7 ; dash dot dot
	    (lambda ()
	      (graphics-draw-line *the-scope* 0 (- .5 trace) 1 (- .5 trace)))))
	  ((= xspan 0)
	   (graphics-bind-line-style *the-scope* 6 ; center dash
	    (lambda ()
	      (graphics-draw-line *the-scope* 0 (- .5 trace) 1 (- .5 trace)))))
	  (else
	   (if (<= minx 0 maxx)
	       (graphics-bind-line-style *the-scope* 2 ; dots
		 (lambda ()
		   (let ((x0 (- (/ minx xspan))))
		     (graphics-draw-line *the-scope*
					 x0 (- trace) x0 (- (- trace 1)))))))
	   (let ((y0 (- yoffset (/ miny yspan) trace)))
	     (if (<= miny 0 maxy)
		 (graphics-bind-line-style *the-scope* 2 ; dots
		   (lambda ()
		     (graphics-draw-line *the-scope* 0 y0 1 y0))))
	     (for-each
	      (lambda (x y)
		(let ((xp (/ (- x minx) xspan))
		      (yp (- (+ (/ (- y miny) yspan) yoffset) trace)))
		  (if dots?
		      (graphics-draw-point *the-scope* xp yp)
		      (graphics-draw-line *the-scope* xp y0 xp yp))))
	      xdata
	      ydata))))
    (graphics-draw-text *the-scope* 0 (- .6 trace) (number->string trace))
    (list trace (list minx maxx miny maxy))))				  
