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

;;;; Current stable of ODE integrators

;;; Needs GENERAL/TABLE 
(load "advance" scmutils-base-environment)

;;; Needs NUMERICS/LINEAR/GAUSS-JORDAN
(load "qc" scmutils-base-environment)

(load "bulirsch-stoer" scmutils-base-environment)

;;; Needs NUMERICS/LINEAR/FULL-PIVOT
(load "gear" scmutils-base-environment)

(load "ode-advancer" scmutils-base-environment)

(load "interface" scmutils-base-environment)

