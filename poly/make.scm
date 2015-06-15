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

;;; Polynomial stuff

(error "This file (poly/make.scm) shouldn't be loaded.")

(if (lexical-unreferenceable? user-initial-environment
			      'scmutils-base-environment)
    (with-working-directory-pathname
     "../kernel/"
     (lambda ()
       (load "make" user-initial-environment))))

(in-package user-initial-environment

  (define (polysys:make subsystem-name)

    (case subsystem-name
      ((root-finder poly->roots poly:->roots) 
       (if (lexical-unreferenceable? scmutils-base-environment 'root-finder-package)
	   (load "make-rootfinder" scmutils-base-environment)))
      (else
       (error "I don't know about" subsystem-name)))))

