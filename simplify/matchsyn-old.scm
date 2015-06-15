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

;;;; Syntax definitions for rule and matcher language

(declare (usual-integrations))

;;;   This file (or a compiled version of it) must be loaded into
;;;   the compiler before any match language code is compiled.

;;;  The match special forms defined are in upper case.
;;;  The variables assumed to be provided are:
;;;   *expression* *dictionary* *fail* *succeed*

(syntax-table-define system-global-syntax-table 'RULE-SYSTEM
  (lambda rules
    `(rule-simplifier ',rules)))


(define (concat-actions actions)
  (cond ((null? actions)
	 `(*succeed* *dictionary* *fail*))
	((null? (cdr actions))
	 `(,(car actions)
	   *fail*
	   *succeed*))
	(else
	 `(,(car actions)
	   *fail*
	   (lambda (*dictionary* *fail*)
	     ,(concat-actions (cdr actions)))))))


(syntax-table-define system-global-syntax-table 'MATCHER-PROCEDURE
  (lambda (match-action)
    `(lambda (*expression* *dictionary* *fail* *succeed*)
       ,match-action)))

(define (match-maker pattern predicate actions)
  `(match ',pattern *expression* *dictionary*
	  *fail*
	  (lambda (*dictionary* *fail*)
	    ,(if (eq? predicate #t)
		 (concat-actions actions)
		 `(if (not ,predicate)
		      (*fail*)
		      ,(concat-actions actions))))))


(syntax-table-define system-global-syntax-table 'MATCHES
  (lambda (pattern #!optional predicate #!rest actions)
    (match-maker pattern predicate actions)))

(define (match-disjunction-maker match-clauses)
  (cond ((null? match-clauses)
	 `(*fail*))
	((null? (cdr match-clauses))
	 (let ((clause (car match-clauses)))
	   (let ((pattern (car clause))
		 (predicate (cadr clause))
		 (actions (cddr clause)))
	     (match-maker pattern predicate actions))))
	(else
	   `(let ((*fail*
		   (lambda ()
		     ,(match-disjunction-maker (cdr match-clauses)))))
	      ,(let ((clause (car match-clauses)))
		 (let ((pattern (car clause))
		       (predicate (cadr clause))
		       (actions (cddr clause)))
		   (match-maker pattern predicate actions)))))))


(syntax-table-define system-global-syntax-table 'MATCHES-ONE-OF
  (lambda match-clauses
    (match-disjunction-maker match-clauses)))



(syntax-table-define system-global-syntax-table 'MATCH-ASSIGN
  (lambda (match-variable value-expression)
    `(match-assign-element ',match-variable
			   *dictionary*
			   ,value-expression)))

(syntax-table-define system-global-syntax-table ':
  (lambda (match-variable)
    `(match-get-value ',match-variable *dictionary*)))

(syntax-table-define system-global-syntax-table '::
  (lambda (match-variable)
    `(match-get-value ',match-variable *dictionary*)))

