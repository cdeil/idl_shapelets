function comb,n,r

;$Id: comb.pro, v2$
;
; Copyright © 2005 Richard Massey and Alexandre Refregier.
;
; This file is a part of the Shapelets analysis code.
; www.astro.caltech.edu/~rjm/shapelets/
;
; The Shapelets code is free software; you can redistribute it and/or
; modify it under the terms of the GNU General Public Licence as published
; by the Free Software Foundation; either version 2 of the Licence, or
; (at your option) any later version.
;
; The Shapelets code is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
; GNU General Public Licence for more details.
;
; You should have received a copy of the GNU General Public Licence
; along with the Shapelets code; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
;
;+
; NAME:
;       COMB
;
; PURPOSE:
;       Calculates binomial coefficients C(n,r) = (n)
;                                                 (r)
;
; CATEGORY:
;       Mathematics.
;
; CALLING SEQUENCE:
;       Result = comb(n,r)
;
; INPUTS:
;       n:    A non-negative scalar or array of values.
;       r:    A non-negative scalar or array of values.
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       comb: Floating point scalar or array.
;
; EXAMPLE:
;       Calculate number of ways of choosing 2 elements from an (unordered)
;       set of 5.
;         result=comb(5,2)
;
;       There are 10 distinct ways of choosing these elements.
;
; NOTES:
;       Will return a floating point number even if integers would do.
;       This is beause of factorial.pro.
;
; PROCEDURES USED:
;       factorial - calculate x!
;
; MODIFICATION HISTORY:
;       Apr 2005 - Generalised to accept arrays of coefficients by RM
;       Jul 2003 - Header added by RM
;       Feb 2002 - Written by Richard Massey
;-

COMPILE_OPT HIDDEN, idl2

if (min(n)<min(r)<0) then message, 'Values for n and r must be non-negative.'
if min((n+1.)/(r+1.)) lt 1 then message, 'Values for r must be less than those for n.'

return,factorial(n)/factorial(n-r)/factorial(r)

end
