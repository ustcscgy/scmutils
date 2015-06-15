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

;;;; Scmutils top-level loader

(ge user-initial-environment)

(add-identification! "ScmUtils" "Mechanics " " Summer 2002")

(define scmutils-base-environment
  user-initial-environment)

(load "general/comutils" scmutils-base-environment)

(start-canonicalizing-symbols!)

;;; LOCAL-ASSIGNMENT should eventually be replaced with
;;; ENVIRONMENT-DEFINE when that has stabilized.

(local-assignment scmutils-base-environment
		  '*environment*
		  'scmutils-base-environment)

(local-assignment scmutils-base-environment
		  'derivative-symbol
		  (string->symbol "D"))

;;; EXTEND-IC-ENVIRONMENT should eventually be replaced with
;;; EXTEND-INTERPRETER-ENVIRONMENT.

(let ((numerical-environment
       (extend-ic-environment scmutils-base-environment)))
  (local-assignment scmutils-base-environment
		    'numerical-environment
		    numerical-environment)
  (local-assignment numerical-environment
		    '*environment*
		    'numerical-environment))

(define numerical-environment
  (access numerical-environment scmutils-base-environment))

(define (in-scmutils-directory relative-path thunk)
  (with-working-directory-pathname
      (merge-pathnames (pathname-as-directory relative-path)
		       (directory-pathname (current-load-pathname)))
    thunk))

(load-option 'hash-table)
(load-option 'synchronous-subprocess)

;;; This doesn't work because load must also get dll's into the microcode.
;;;(load-option 'swat)

(in-scmutils-directory "./general"
		       (lambda ()
			 (load "load" scmutils-base-environment)))
(in-scmutils-directory "./kernel"
		       (lambda ()
			 (load "load" scmutils-base-environment)))

(define generic-environment
  (access generic-environment scmutils-base-environment))

(in-scmutils-directory "./simplify"
		       (lambda ()
			 (load "load" scmutils-base-environment)))

(define symbolic-environment
  (access symbolic-environment scmutils-base-environment))
(define rule-environment
  (access rule-environment scmutils-base-environment))

(in-scmutils-directory "./display"
		       (lambda ()
			 (load "load" scmutils-base-environment)))

(in-scmutils-directory "./enclose"
		       (lambda ()
			 (load "load" scmutils-base-environment)))
(in-scmutils-directory "./numerics"
		       (lambda ()
			 (load "load" scmutils-base-environment)))
(in-scmutils-directory "./poly"
		       (lambda ()
			 (load "load" scmutils-base-environment)))


(start-preserving-case!)

(in-scmutils-directory "./kernel"
		       (lambda ()
			 (load "litfun" scmutils-base-environment)))
      
(in-scmutils-directory "./mechanics"
		       (lambda ()
			 (load "load" scmutils-base-environment)))

(ge generic-environment)


