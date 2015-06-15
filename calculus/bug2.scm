#| -*-Scheme-*-

Copyright (C) 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994,
    1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005,
    2006, 2007, 2008, 2009, 2010, 2011 Massachusetts Institute of
    Technology

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


(cd "/usr/local/scmutils/src/calculus")
#| #[pathname 14 "/usr/local/scmutils/src/calculus/"] |#

(define (T V) V)
#| T |#

(declare-argument-types! T (list vector-field?))
#| done |#

(stop-scmutils-print)
;Value: #[compiled-procedure 20 (scmutils/repl-write "custom-repl" #x1) #x1a #x3a9e102]

((((covariant-derivative (literal-Cartan 'G R2-rect))
   (literal-vector-field 'X R2-rect))
  T)
 (literal-vector-field 'U R2-rect))
;Value: #[apply-hook 21]

(vector-field? #@21)
;Value: #t

((literal-1form-field 'omega R2-rect)
 ((((covariant-derivative (literal-Cartan 'G R2-rect))
    (literal-vector-field 'X R2-rect))
   T)
  (literal-vector-field 'U R2-rect)))
;Value: #[compiled-closure 31 (lambda "mathutil" #x1c) #x12ca #x3d3b442 #x3ab88400]

(((literal-1form-field 'omega R2-rect)
  ((((covariant-derivative (literal-Cartan 'G R2-rect))
     (literal-vector-field 'X R2-rect))
    T)
   (literal-vector-field 'U R2-rect)))
 (typical-point R2-rect))
;Bad point: rectangular
