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

;;;;          New matcher language instructive sample

;;;  The language is implemented with new special forms, sigh.
;;;  The special forms are:

;;; MATCHER-PROCEDURE, :, ::, MATCHES-ONE-OF, and MATCH-ASSIGN.  Note
;;;  that ? and ??  never appear anywhere but in a matcher pattern, so
;;;  they need not be special forms.

;;; There are pattern procedures, as in Hewitt's PhD thesis.
;;; Pattern variables are either symbols or calls to pattern
;;; procedures.

;;; Some variables match elements -- these are indicated by "?".
;;; Others match segments (sublists) -- these are indicated by "??".

(define sincos-flush-ones
  (rule-system
   ((+ (?? a1) (? (sin-cos-sq)) (?? a2) (? (cos-sin-sq)) (?? a3))
    none
    (evaluate
     (apply +
	    (append (list (: extra-term))
		    (:: a1)
		    (:: a2)
		    (:: a3)))))
   ))

;;; The pattern variables may be restricted by predicates.
;;; Restrictions on the matched item may be placed in the pattern.
;;; There may also be restrictions that involve multiple pattern
;;; variables.  These come after the pattern.  If there is no such
;;; restriction we put in the placeholder "#t".

;;; Values of matcher variables may be accessed by ":".  The list of
;;; elements matched by a segment variable may be obtained by "::".

;;; A match may assign a match variable to a value, if it does not
;;; already have one.

(define sin-cos-sq
  (matcher-procedure
    (matches-one-of
       ((expt (sin (? x)) (? n at-least-two?))
	#t
	(match-assign opposite
		      (cos (: x)))
	(match-assign extra-term
		      (expt (sin (: x)) (- (: n) 2))))

       ((expt (cos (? x)) (? n at-least-two?))
	#t
	(match-assign opposite
		      (sin (: x)))
	(match-assign extra-term
		      (expt (cos (: x)) (- (: n) 2))))

       ((* (?? f1) (expt (sin (? x)) (? n at-least-two?)) (?? f2))
	#t
	(match-assign opposite
		      (cos (: x)))
	(match-assign extra-term
		      (apply *
			     (append (:: f1)
				     (list (expt (sin (: x))
						 (- (: n) 2)))
				     (:: f2)))))

       ((* (?? f1) (expt (cos (? x)) (? n at-least-two?)) (?? f2))
	#t
	(match-assign opposite
		      (sin (: x)))
	(match-assign extra-term
		      (apply *
			     (append (:: f1)
				     (list (expt (cos (: x))
						 (- (: n) 2)))
				     (:: f2)))))
       )))

(define (at-least-two? n)
  (and (number? n) (>= n 2)))


(define cos-sin-sq
  (matcher-procedure
   (matches-one-of
    ((expt (? opposite) (? m at-least-two?))
     (exact-zero?
	(rcf:simplify
	 (- (: extra-term)
	    (expt (: opposite) (- (: m) 2))))))
    ((* (?? fs1) (expt (? opposite) (? m at-least-two?)) (?? fs2))
     (exact-zero?
      (rcf:simplify
       (- (: extra-term)
	  (apply *
		 (append (:: fs1)
			 (list (expt (: opposite)
				     (- (: m) 2)))
			 (:: fs2)))))))
    )))
