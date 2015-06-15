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

(define (canonicalize-input-filename filename)
  (->namestring (canonicalize-input-pathname filename)))

(define (canonicalize-input-pathname filename)
  (let ((pathname (->pathname filename)))
    (or (pathname->input-truename pathname)
	(canonicalize-input-pathname
	 (error:file-operation pathname
			       "find"
			       "file"
			       "file does not exist"
			       canonicalize-input-pathname
			       (list filename))))))

(define (pathname->input-truename pathname)
  (let ((pathname (merge-pathnames pathname)))
    (and (eq? true (file-exists? pathname))
	 pathname)))

(define (canonicalize-output-filename filename)
  (->namestring (canonicalize-output-pathname filename)))

(define (canonicalize-output-pathname filename)
  (pathname->output-truename (->pathname filename)))

(define (pathname->output-truename pathname)
  (merge-pathnames pathname))

(define (canonicalize-overwrite-filename filename)
  (->namestring (canonicalize-overwrite-pathname filename)))

(define (canonicalize-overwrite-pathname filename)
  (pathname->overwrite-truename (->pathname filename)))

(define (pathname->overwrite-truename pathname)
  (merge-pathnames pathname))

(define (pathname-components pathname receiver)
  (receiver (pathname-host pathname)
	    (pathname-device pathname)
	    (pathname-directory pathname)
	    (pathname-name pathname)
	    (pathname-type pathname)
	    (pathname-version pathname)))

(define (pathname-relative? pathname pathname*)
  (let ((diff (enough-pathname pathname pathname*)))
    (and (not (equal? (pathname-directory pathname) (pathname-directory diff)))
	 (make-pathname (pathname-host diff)
			(pathname-device diff)
			(pathname-directory diff)
			(pathname-name pathname)
			(pathname-type pathname)
			(pathname-version pathname)))))

(define home-directory-pathname user-homedir-pathname)
(define init-file-truename init-file-pathname)
(define pathname->absolute-pathname merge-pathnames)
(define pathname->string ->namestring)
(define pathname-directory-path directory-pathname)
(define pathname-directory-string directory-namestring)
(define pathname-name-path file-pathname)
(define pathname-name-string file-namestring)
(define string->pathname ->pathname)