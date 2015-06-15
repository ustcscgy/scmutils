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

(define *opaque-procedure-table* '())

(define (make-opaque name #!optional function-type)
  (if (default-object? function-type)
      (set! function-type (default-function-type 1)))
  (cond ((assq name *opaque-procedure-table*)
	 =>
	 (lambda (entry)
	   (let ((value (environment-lookup generic-environment name)))
	     (if (and (not (literal-function? value))
		      (or (not (eq? value (cadr entry)))
			  (not (equal? function-type (caddr entry)))))
		 (set-car! (cdr entry) value)))))
	((environment-bound? generic-environment name)
	 (set! *opaque-procedure-table*
	       (cons (list name
			   (environment-lookup generic-environment name)
			   function-type)
		     *opaque-procedure-table*)))
	(else (error "Cannot find procedure to opaqify" name)))
  (eval `(define ,name (literal-function ',name ',function-type))
	generic-environment)
  name)

(define (make-transparent name)
  (cond ((assq name *opaque-procedure-table*)
	 =>
	 (lambda (entry)
	   (environment-assign! generic-environment name (cadr entry))))
	(else 'done)))
	 
(define (compile-opaque name)
  (eval `(define ,name) numerical-environment)
  (cond ((assq name *opaque-procedure-table*)
	 =>
	 (lambda (entry)
	   (let ((procedure (cadr entry))
		 (function-type (caddr entry)))
	     ;; Must use function type
	     (let ((arity (procedure-arity procedure)))
	       (if (not (eq? (car arity) (cdr arity)))
		   (error "I don't know how to compile this kind of procedure"
			  name))
	       (let ((cproc
		      (lambda->numerical-procedure
		       (lambdafy (car arity)
			 (lambda names
			   (simplify (apply procedure names)))))))
		 (environment-assign! numerical-environment name cproc))))))
	(else (error "No opaque definition for procedure" name))))
