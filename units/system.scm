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

;;;; Unit systems

(define (define-unit-system system-name #!rest base-units)
  (if (environment-bound? scmutils-base-environment system-name)
      (write-line `(clobbering ,system-name)))
  (let ((n (length base-units)))    
    (environment-define scmutils-base-environment '*unitless*
			(make-vector n 0))
    (environment-define scmutils-base-environment 'unitless
			(make-unit system-name
				   (access *unitless*
					   scmutils-base-environment)
				   1))
    (environment-define scmutils-base-environment '*angular* ; may be reset by SI-units
			(access *unitless* scmutils-base-environment))
    (environment-define scmutils-base-environment 'angular
			(make-unit system-name
				   (access *angular*
					   scmutils-base-environment)
				   1))
    (let ((base-specs
	   (map (lambda (base-spec i)
		  (let* ((unit-name (car base-spec))
			 (exponents
			  (make-initialized-vector n
						   (lambda (j)
						     (if (fix:= i j) 1 0))))
			 (unit (make-unit system-name exponents 1)))
		    (if (environment-bound? scmutils-base-environment
					    unit-name)
			(write-line `(clobbering ,unit-name)))
		    (environment-define scmutils-base-environment
					unit-name
					unit)
		    (append base-spec (list unit))))
		base-units
		(iota n))))
      (environment-define scmutils-base-environment
			  system-name
			  (list '*unit-system*
				system-name
				base-specs          ;base units
				'()	            ;derived units
				'()	            ;additional units
				))))
  system-name)


(define (unit-system? system)
  (and (pair? system)
       (eq? (car system) '*unit-system*)))

(define (unit-system-name system)
  (cadr system))

(define (base-units system)
  (caddr system))

(define (derived-units system)
  (cadddr system))

(define (alternate-units system)
  (car (cddddr system)))

(define (define-derived-unit system unit-name tex description content
	                     #!optional scale-factor)
  (assert (unit-system? system))
  (if (environment-bound? scmutils-base-environment unit-name)
      (write-line `(clobbering ,unit-name)))
  (if (default-object? scale-factor)
      (set! scale-factor 1))
  (set! content
	(make-unit (unit-system-name system)
		   (unit-exponents content)
		   (* scale-factor (unit-scale content))))
  (let ((unit-spec (list unit-name tex description content)))
    (define-derived-unit! system unit-spec)
    (environment-define scmutils-base-environment unit-name content)
    unit-name))

(define (define-derived-unit! system unit-spec)
  (set-car! (cdddr system)
	    (append (cadddr system)
		    (list unit-spec))))




(define (define-additional-unit system unit-name tex description content
	                        #!optional scale-factor)
  (assert (unit-system? system))
  (if (environment-bound? scmutils-base-environment unit-name)
      (write-line `(clobbering ,unit-name)))
  (if (default-object? scale-factor)
      (set! scale-factor 1))
  (set! content
	(make-unit (unit-system-name system)
		   (unit-exponents content)
		   (* scale-factor (unit-scale content))))
  (let ((unit-spec (list unit-name tex description content)))
    (define-additional-unit! system unit-spec)
    (environment-define scmutils-base-environment unit-name content)
    unit-name))

(define (define-additional-unit! system unit-spec)
  (set-car! (cddddr system)
	    (append (car (cddddr system))
		    (list unit-spec))))



(define (define-multiplier name tex-string value)
  (if (environment-bound? scmutils-base-environment name)
      (write-line `(clobbering ,name)))
  (environment-define scmutils-base-environment name value))

(define *constants* '())

(define *numerical-constants* #t)

(define (define-constant name tex-string description value units
	                #!optional uncertainty)
  (if (environment-bound? scmutils-base-environment name)
      (write-line `(clobbering ,name)))
  (let ((constant (literal-number name)))
    (cond ((with-units? value)
	   (assert (equal? (u:units value) units)))
	  ((unitless? units))
	  (else (set! value (with-units value units))))
    (add-property! constant 'numerical-value value)
    (add-property! constant 'units units)
    (add-property! constant 'tex-string tex-string)
    (add-property! constant 'description description)
    (if (not (default-object? uncertainty))
	(add-property! constant 'uncertainty uncertainty))
    (if *numerical-constants*
	(environment-define scmutils-base-environment name value)
	(environment-define scmutils-base-environment name constant))
    (set! *constants* (cons constant *constants*))
    name))


;;; & is used to attach units to a number, or to check that a number
;;; has the given units.

(define (& value u1 #!optional u2)
  (let ((units (if (default-object? u2) u1 u2))
	(scale (if (default-object? u2) 1 u1)))
    (assert (and (not (units? value)) (number? scale) (units? units)))
    (if (with-units? value)
	(if (equal? (unit-exponents units)
		    (unit-exponents (u:units value)))
	    value
	    (error "Units do not match: &" value units))
	(with-units (g:* scale (unit-scale units) value)
	  (make-unit (unit-system units)
		     (unit-exponents units)
		     1)))))

(define *unit-constructor* '&)

(define (with-units->expression system num)
  (assert (unit-system? system))
  (cond ((with-units? num)
	 (let ((value (g:* (unit-scale (u:units num)) (u:value num)))
	       (vect (unit-exponents (u:units num))))
	   (list *unit-constructor*
		 value
		 (exponent-vector->unit-expression system vect))))
	((units? num)
	 (list *unit-constructor*
	       (unit-scale num)
	       (exponent-vector->unit-expression system
						 (unit-exponents num))))
	(else num)))

(define (exponent-vector->unit-expression system vect)
  (or (find-unit vect (base-units system))
      (find-unit vect (derived-units system))
      (unit-expresson (vector->list vect)
		      (map car (base-units system)))))

(define (find-unit vect ulist)
  (let ((v
	 (find-matching-item ulist
	   (lambda (entry)
	     (equal? (unit-exponents (list-ref entry 3))
		     vect)))))
    (if v (car v) #f)))

(define (unit-expresson exponents base-unit-names)
  (cons '*
	(apply append
	       (map (lambda (exponent base-name)
		      (cond ((g:zero? exponent) '())
			    ((g:one? exponent) (list base-name))
			    (else
			     (list (list 'expt base-name exponent)))))
		    exponents
		    base-unit-names))))


#|
(with-units->expression SI foot)
;Value: (& .3048 meter)

(with-units->expression SI (& 2 foot))
;Value: (& .6096 meter)

(with-units->expression SI (/ (* :k (& 300 kelvin)) :e))
;Value: (& .02585215707677003 volt)

(with-units->expression SI :c)
;Value: (& 299792458. (* meter (expt second -1)))

(with-units->expression SI :h)
;Value: (& 6.6260755e-34 (* (expt meter 2) kilogram (expt second -1)))
|#