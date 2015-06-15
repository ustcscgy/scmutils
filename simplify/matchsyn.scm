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

;;;; Syntax definitions for rule and matcher language

(declare (usual-integrations))

;;;   This file (or a compiled version of it) must be loaded into
;;;   the compiler before any match language code is compiled.

;;;  The match special forms defined are in upper case.
;;;  The variables assumed to be provided are:
;;;   *expression* *dictionary* *fail* *succeed*

(define-syntax rule-system
  (sc-macro-transformer
   (lambda (form environment)
     environment
     `(rule-simplifier ',(cdr form)))))


(define-syntax matcher-procedure
  (sc-macro-transformer
   (lambda (form environment)
     (if (syntax-match? '(EXPRESSION) (cdr form))
	 `(lambda (*expression* *dictionary* *fail* *succeed*)
	    ,(make-syntactic-closure environment
				     '(*expression* *dictionary* *fail* *succeed*)
				     (cadr form)))
	 (ill-formed-syntax form)))))

(define-syntax matches
  (rsc-macro-transformer
   (lambda (form environment)
     (if (syntax-match? '(DATUM EXPRESSION * EXPRESSION) (cdr form))
	 (match-maker (cadr form)
		      (close-syntax (caddr form) environment)
		      (map (lambda (action)
			     (close-syntax action environment))
			   (cdddr form))
		      environment)
	 (ill-formed-syntax form)))))

(define (match-maker pattern predicate actions environment)
  `(,(close-syntax 'match environment)
    ',pattern *expression* *dictionary*
    *fail*
    (,(close-syntax 'lambda environment)
     (*dictionary* *fail*)
     ,(if (eq? predicate #t)
	  (concat-actions actions environment)
	  `(,(close-syntax 'if environment)
	    (,(close-syntax 'not environment) ,predicate)
	    (*fail*)
	    ,(concat-actions actions environment))))))

(define (concat-actions actions environment)
  (cond ((null? actions)
	 `(*succeed* *dictionary* *fail*))
	((null? (cdr actions))
	 `(,(car actions)
	   *fail*
	   *succeed*))
	(else
	 `(,(car actions)
	   *fail*
	   (,(close-syntax 'lambda environment)
	    (*dictionary* *fail*)
	    ,(concat-actions (cdr actions) environment))))))

(define-syntax matches-one-of
  (rsc-macro-transformer
   (lambda (form environment)
     (match-disjunction-maker (cdr form) environment))))

(define (match-disjunction-maker match-clauses environment)
  (cond ((null? match-clauses)
	 `(*fail*))
	((null? (cdr match-clauses))
	 (let ((clause (car match-clauses)))
	   (let ((pattern (car clause))
		 (predicate (cadr clause))
		 (actions (cddr clause)))
	     (match-maker pattern predicate actions environment))))
	(else
	   `(,(close-syntax 'let environment)
	     ((*fail*
	       (,(close-syntax 'lambda environment)
		()
		,(match-disjunction-maker (cdr match-clauses)
					  environment))))
	     ,(let ((clause (car match-clauses)))
		(let ((pattern (car clause))
		      (predicate (cadr clause))
		      (actions (cddr clause)))
		  (match-maker pattern predicate actions environment)))))))

(define-syntax match-assign
  (sc-macro-transformer
   (lambda (form environment)
     (if (syntax-match? '(DATUM EXPRESSION) (cdr form))
	 `(match-assign-element ',(cadr form)
				,(close-syntax '*dictionary* environment)
				,(close-syntax (caddr form) environment))
	 (ill-formed-syntax form)))))

(define-syntax :
  (sc-macro-transformer
   (lambda (form environment)
     (if (syntax-match? '(DATUM) (cdr form))
	 `(match-get-value ',(cadr form)
			   ,(close-syntax '*dictionary* environment))
	 (ill-formed-syntax form)))))

(define-syntax ::
  (sc-macro-transformer
   (lambda (form environment)
     (if (syntax-match? '(DATUM) (cdr form))
	 `(match-get-value ',(cadr form)
			   ,(close-syntax '*dictionary* environment))
	 (ill-formed-syntax form)))))
