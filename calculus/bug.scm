#| -*-Scheme-*-

Copyright (C) 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994,
    1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005,
    2006, 2007, 2008, 2009, 2010 Massachusetts Institute of Technology

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


(define bar
  '(+ (* x0 z0 (((partial 1) ((partial 2) g)) #(x0 y0 z0)))
      (* -1 x0 z0 (((partial 2) ((partial 1) g)) #(x0 y0 z0)))))
;Value: bar

(pp
 ((access new-simplify rule-environment)
  bar))
(+ (* x0 z0 (((partial 1) ((partial 2) g)) #(x0 y0 z0)))
   (* -1 x0 z0 (((partial 1) ((partial 2) g)) #(x0 y0 z0))))
;Unspecified return value

(pp (((access simplify-and-canonicalize rule-environment)
      (access canonicalize-partials rule-environment)
      (access simplify-and-flatten rule-environment))
     bar))
0
;Unspecified return value

(pp ((lambda (exp)
       ((compose ((access only-if rule-environment)
		  (memq 'partial (variables-in exp))
		  ((access simplify-and-canonicalize rule-environment)
		   (access canonicalize-partials rule-environment)
		   (access simplify-and-flatten rule-environment)))
		 (access simplify-and-flatten rule-environment))
	exp))
     bar))
0
;Unspecified return value

(ge rule-environment)
;Value: #[environment 20]

(pp ((lambda (exp)
       ((compose (only-if (memq 'partial (variables-in exp))
			  (simplify-and-canonicalize canonicalize-partials
						     simplify-and-flatten))
		 simplify-and-flatten)
	exp))
     (access bar generic-environment)))
0
;Unspecified return value

(pp (new-simplify (access bar generic-environment)))
(+ (* x0 z0 (((partial 1) ((partial 2) g)) #(x0 y0 z0)))
   (* -1 x0 z0 (((partial 1) ((partial 2) g)) #(x0 y0 z0))))
;Unspecified return value

;Value: new-simplify

(pp (new-simplify (access bar generic-environment)))
0
;Unspecified return value
