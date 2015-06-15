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

;;;; Scmutils kernel loader

;;; Useful universal utilities

(load "numeric"  scmutils-base-environment)
(load "utils"    scmutils-base-environment)
(load "iterat"   scmutils-base-environment)
(load "express"  scmutils-base-environment)

;;; The following define the generic operators

(load "ghelper"  scmutils-base-environment)
(load "generic"  scmutils-base-environment)
(load "mathutil" scmutils-base-environment)

;;; Magic apply extension to allow application
;;;  of things, such as vectors and numbers, that
;;;  are not legal Scheme procedures.

(load "extapply" scmutils-base-environment)
;;; Disable for system debugging.
;;; (set! *enable-generic-apply* #f)
;;; Enable for mechanics.
(set! *enable-generic-apply* #t)

;;; GHELPER is needed to load specific types
;;;  Lookup is in reverse order, so put numbers last

;;; Support for loading types.
(load "types"    scmutils-base-environment)
;;(load "/usr/local/scmutils/src/kernel/types.scm.~21~"    scmutils-base-environment)

;;;(define (diff-memoize-1arg f) f)
;;;(define (diff-memoize-2arg f) f)
;;;(define (diff-memoize f) f)
(define (diff-memoize-1arg f) (linear-memoize-1arg f))
(define (diff-memoize-2arg f) (linear-memoize f))
(define (diff-memoize f) (linear-memoize f))
;;;(define (diff-memoize f) (hash-memoize f))

(load "diff"     scmutils-base-environment)
(load "deriv"    scmutils-base-environment)
(load "operator" scmutils-base-environment)
(load "function" scmutils-base-environment)

(load "matrices" scmutils-base-environment)

(load "modarith" scmutils-base-environment)
(load "numbers"  scmutils-base-environment)
(load "vectors"  scmutils-base-environment)
(load "quaternion" scmutils-base-environment)

(load "strutl"   scmutils-base-environment)
(load "pseries"  scmutils-base-environment)


;;; The following must be loaded late, to enable  special-case SQUARE.
(load "structs"  scmutils-base-environment)

;;; Literal-construction.
(load "numsymb"  scmutils-base-environment)

;;; must come after numsymb
(load "heuristic" scmutils-base-environment)

;;; Sets up generic environment
(load "genenv"   scmutils-base-environment)
