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

;;;; Simplifier loader

;;; Canonical simplifiers are based on Rational Canonical Form,
;;;  which is, in turn, based on Polynomial Canonical Form.

;;;The following two files must be loaded in the given order.
(load "pcf"      scmutils-base-environment)
(load "rcf"      scmutils-base-environment)


;;; We need flattened polynomials to support rule-based simplifiers.

(load "fpf" scmutils-base-environment)


;;; Canonical simplifiers are glued together with SIMPLIFY.

(load "simplify" scmutils-base-environment)

(load "split-poly" scmutils-base-environment)

;;; Rule-based simplifiers

(load "symbenv" scmutils-base-environment)

(define rule-environment symbolic-environment)


;;; Rule interpreter

;;; (define (rule-memoize f) f)
(define (rule-memoize f) (linear-memoize-1arg f))
;;; (define (rule-memoize f) (linear-memoize f))
;;; (define (rule-memoize f) (hash-memoize f))

(load "bigsimp" scmutils-base-environment)


;;; Syntax support

(load (if (environment-bound? system-global-environment
			      'make-syntactic-closure)
	  "matchsyn"
	  "matchsyn-old")
      scmutils-base-environment)



;;; Rule systems

(load "sincos" rule-environment)

(load "rules" rule-environment)

(for-each (lambda (name)
	    #|
	    ;; This code exports by copying the binding:
	    (local-assignment scmutils-base-environment
			      name
			      (environment-lookup rule-environment
						  name))
	    |#
	    ;; This code shares the binding:
	    (environment-link-name scmutils-base-environment
				   rule-environment
				   name))
	  '(;; Useful simplifiers
	    ->poisson-form
	    easy-simplify
	    full-simplify
	    trigexpand
	    trigcontract

	    ;; Boolean simplifier controls.
	    log-exp-simplify
	    sqrt-expt-simplify
	    inverse-simplify
	    ignore-zero-simplify
	    commute-partials-simplify))

(load "sparse-load" scmutils-base-environment)
