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

(define (generic-environment-maker)
  (let ((e (extend-ic-environment scmutils-base-environment)))
    (let ((d (lambda (name value)
	       (local-assignment e name value))))
	(d '*environment* 'generic-environment)

	;; Unary operators from generic.scm

	(d 'type g:type)
	(d 'type-predicate g:type-predicate)
	(d 'arity g:arity)

	(d 'inexact? g:inexact?)

	(d 'zero-like g:zero-like)
	(d 'one-like g:one-like)
	(d 'identity-like g:identity-like)
	
	(d 'zero? g:zero?)
	(d 'one? g:one?)
	(d 'identity? g:identity?)

	(d 'negate g:negate)
	(d 'invert g:invert)

	(d 'square g:square)
	(d 'cube   g:cube)

	(d 'sqrt g:sqrt)

	(d 'exp g:exp)
	(d 'log g:log)

	(d 'exp2  g:exp2)
	(d 'exp10 g:exp10)
	(d 'log2  g:log2)
	(d 'log10 g:log10)

	(d 'sin g:sin)
	(d 'cos g:cos)
	(d 'tan g:tan)
	(d 'sec g:sec)
	(d 'csc g:csc)

	(d 'asin g:asin)
	(d 'acos g:acos)

	(d 'sinh g:sinh)
	(d 'cosh g:cosh)
	(d 'tanh g:tanh)
	(d 'sech g:sech)
	(d 'csch g:csch)

	(d 'asinh g:asinh)
	(d 'acosh g:acosh)
	(d 'atanh g:atanh)

	(d 'abs g:abs)

	(d 'determinant g:determinant)
	(d 'trace g:trace)

	(d 'derivative g:derivative)

	;; Binary (and nary) operators from generic.scm

	(d '= g:=)

	(d '+ g:+)
	(d '- g:-)
	(d '* g:*)
	(d '/ g:/)

	(d 'dot-product g:dot-product)

	(d 'expt g:expt)
	(d 'gcd g:gcd)


        ;; Complex operators from generic.scm

	(d 'make-rectangular g:make-rectangular)
	(d 'make-polar g:make-polar)

	(d 'real-part g:real-part)
	(d 'imag-part g:imag-part)
	(d 'magnitude g:magnitude)
	(d 'angle g:angle)

	(d 'conjugate g:conjugate)


	;; Wierd operators from generic.scm

	(d 'atan g:atan)

	(d 'partial-derivative g:partial-derivative)
	(d 'partial g:partial)

	(d 'apply g:apply)


	;; Compound operators from mathutil.scm

	(d 'arg-scale g:arg-scale)
	(d 'arg-shift g:arg-shift)

	(d 'sigma g:sigma)

        (d 'ref   g:ref)
	(d 'size  g:size)

	(d 'compose g:compose)
	)
    e))

(define generic-environment
  (generic-environment-maker))

(define generic-numerical-operators
  '(	
	zero-like
	one-like
	identity-like

	negate
	invert

	square
	cube

	sqrt

	exp
	log

	exp2
	exp10
	log2
	log10

	sin
	cos
	tan
	sec
	csc

	asin
	acos

	sinh
	cosh
	tanh
	sech
	csch

	abs

	+
	-
	*
	/

	expt
	gcd

	make-rectangular
	make-polar

	real-part
	imag-part
	magnitude
	angle

	conjugate

	atan))





