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

(declare (usual-integrations))

;;; Default settings

(define *ode-integration-method* 'BULIRSCH-STOER)

;;; (set! *ode-integration-method* 'QCRK4)
;;; (set! *ode-integration-method* 'BULIRSCH-STOER)
;;; (set! *ode-integration-method* 'QCCTRAP2)

(define *first-step-scale* 1.0)

(define *corrector-convergence-margin* 1.0e-1)

(define *progress-monitor* #f)

(define *last-state*)

(define (ode-advancer sysder state dt local-error-tolerance)
  (case *ode-integration-method*
    ((QCRK4 qcrk4)
     (qcrk4-advancer sysder state dt local-error-tolerance))
    ((QCCTRAP2 qcctrap2)
     (qc-ctrap-advancer sysder state dt local-error-tolerance))
    ((BULIRSCH-STOER bulirsch-stoer Bulirsch-Stoer)
     (bs-advancer sysder state dt local-error-tolerance))
    (else
     (error "Unknown ode integrator" *ode-integration-method*))))


(define (advance-monitor ns step-achieved step-suggested cont)
  (if *progress-monitor* (pp `(,ns ,step-achieved ,step-suggested)))
  (set! *last-state* ns)
  (cont))

(define (final-step-monitor ns step-achieved step-suggested)
  (if *progress-monitor* (pp `(,ns ,step-achieved ,step-suggested)))
  (set! *last-state* ns)
  ns)

(define (bs-advancer sysder state dt local-error-tolerance)
  ((advance-generator
    (bulirsch-stoer-lisptran		;integrator
     (system-derivative->lisptran-derivative sysder)
     (vector-length state)
     local-error-tolerance))
   state
   dt
   (* *first-step-scale* dt)
   dt
   advance-monitor
   final-step-monitor))

(define (qcrk4-advancer sysder state dt local-error-tolerance)
  ((advance-generator
    ((quality-control rk4 4)
     sysder	
     local-error-tolerance))
   state
   dt
   (* *first-step-scale* dt)
   dt
   advance-monitor
   final-step-monitor))

(define (qc-ctrap-advancer sysder state dt local-error-tolerance)
  ((advance-generator
    ((quality-control c-trapezoid 2)
     sysder			
     local-error-tolerance 
     (* *corrector-convergence-margin*
	local-error-tolerance)))
   state
   dt 
   (* *first-step-scale* dt)
   dt
   advance-monitor
   final-step-monitor))
