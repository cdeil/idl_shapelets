pro shapelets_exponentiate_operations, operation, $
                                       structure, $
				       amount,    $
				       order,     $
				       extend=extend, $
				       _extra = ex

;$Id: shapelets_exponentiate_operations.pro, v2$
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
;      SHAPELETS_EXPONENTIATE_OPERATIONS
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Applies a (pure) shear to an object, using ladder operators in
;      (Cartesian) shapelet space.
;
; EXAMPLE USE:
;      shapelets_exponentiate_operations, "shapelets_shear", decomp, gamma, 3
;
; INPUTS:
;      OPERATION - String containing the routine to be called.
;      STRUCTURE - A shapelet decomp structure or shapecat structure.
;      AMOUNT    - Amount of shear/translation/whatever the operation is.
;      ORDER     - Maximum order in the expansion.
;
; OPTIONAL INPUTS:
;      None. 
;
; KEYWORD PARAMETERS:
;      None. 
;
; OUTPUTS:
;      Performs a higher-order operation on the object(s) in structure.
;
; TO DO:
;      Error propagation.
;
; MODIFICATION HISTORY:
;      Apr 05 - Written by Richard Massey
;-

COMPILE_OPT idl2, HIDDEN

; Extend n_max to cope with the higher order crosstalk 
case strupcase(operation) of 
  "SHAPELETS_SHEAR": n_talk=2
  "SHAPELETS_DILATE": n_talk=2
  "SHAPELETS_TRANSLATE": n_talk=1
  else: n_talk=1
endcase
if not keyword_set(extend) then extend=0
extension=(fix(order)*n_talk)>extend
shapelets_extend_nmax, structure, extension

; Perform operation to first-order multiple times
case fix(order) of
  1: begin
       ; 1+gS = 1+gS (shouldn't really be used)
       call_procedure, operation, structure, amount, order=1, _extra= ex, /NOHISTORY
     end
  2: begin
       ; 1+gS+gS/2 = 1/2 + 1/2(1+gS)^2
       term1=structure.coeffs/2.
       for i=1,2 do call_procedure, operation, structure, amount, order=1, _extra= ex, /NOHISTORY
       term2=structure.coeffs/2
       structure.coeffs=term1+term2
     end
  3: begin
       ; 1+gS+gS/2!+gS/3! = 1/3 + 1/2(1+gS) + 1/6(1+gS)^3
       term1=structure.coeffs/3.
       call_procedure, operation, structure, amount, order=1, _extra= ex, /NOHISTORY
       term2=structure.coeffs/2.
       for i=2,3 do call_procedure, operation, structure, amount, order=1, _extra= ex, /NOHISTORY
       term3=structure.coeffs/6.
       structure.coeffs=term1+term2+term3
     end
  4: begin
       ; 1+gS+gS/2!+gS/3!+gS/4! = 3/8 + 1/3(1+gS) + 1/4(1+gS)^2 + 1/24(1+gS)^4
       term1=structure.coeffs*3/8.
       call_procedure, operation, structure, amount, order=1, _extra= ex, /NOHISTORY
       term2=structure.coeffs/3.
       call_procedure, operation, structure, amount, order=1, _extra= ex, /NOHISTORY
       term3=structure.coeffs/4.
       for i=3,4 do call_procedure, operation, structure, amount, order=1, _extra= ex, /NOHISTORY
       term4=structure.coeffs/24.
       structure.coeffs=term1+term2+term3+term4
     end
  else: message," Higher order transformations not yet coded beyond 4th order!"
endcase

; Reduce n_max back to input level 
shapelets_extend_nmax, structure, extend-extension

end

