#| -*-Scheme-*-

$Id: copyright.scm,v 1.5 2005/09/25 01:28:17 cph Exp $

Copyright 2005 Massachusetts Institute of Technology

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

(for-each cf-conditionally
	  '("pcf"
	    "rcf"
	    "simplify"
	    "split-poly"
	    "fpf"

	    "symbenv"
	    "bigsimp"
	    ))

(let ((fn
       (if (environment-bound? system-global-environment
			       'make-syntactic-closure)
	   "matchsyn"
	   "matchsyn-old")))
  (if (not (file-processed? fn "scm" "com"))
      (cf fn))
  (let ((environment (nearest-repl/environment)))
    (load fn environment)
    (for-each (lambda (name)
		(link-variables system-global-environment name
				environment name))
	      '(rule-system
		matcher-procedure
		matches
		matches-one-of
		match-assign
		:
		::))))

(for-each cf-conditionally
	  '("sincos"
	    "rules"
	    ))

(for-each cf-conditionally
	  '( "sparse"
	     "sparse-interpolate"
	     "sparse-gcd"
	     
	     "pcf-fpf"))

