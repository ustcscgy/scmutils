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

;;; Extension of Scheme for application of non-procedures in operator position.

(define *enable-generic-apply* true)

(define inapplicable-object/operator
  (condition-accessor condition-type:inapplicable-object 'DATUM))

(define (apply-extension-init)
  (bind-default-condition-handler
   (list condition-type:inapplicable-object)
   (lambda (condition)
     (if *enable-generic-apply*
	 ((stack-frame->continuation
	   (stack-frame/next
	    (stack-frame/next
	     (stack-frame/next
	      (continuation->stack-frame
	       (condition/continuation condition))))))
	  (lambda args
	    (g:apply (inapplicable-object/operator condition) args)))))))

(define (once-only! thunk name)
  (if (lexical-unbound? system-global-environment name)
      (begin (thunk)
	     ;; Create NAME in SGE
	     (eval `(define ,name #t) system-global-environment)
	     'done)
      'already-done))

(once-only! apply-extension-init 'apply-extension-init)
