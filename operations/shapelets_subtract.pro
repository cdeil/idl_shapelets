pro shapelets_subtract, decomp_in1, decomp_in2, NOHISTORY=nohistory

;$Id: shapelets_subtract.pro, v2$
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
;      SHAPELETS_SUBTRACT
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Subtracts one image from another (when both are stored as a shapelet
;      decomposition).
;
; INPUTS:
;      DECOMP_IN1 - A Cartesian shapelet decomposition.
;      DECOMP_IN2 - A Cartesian shapelet decomposition.
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      DECOMP_OUT - A Cartesian shapelet decomposition, DECOMP_IN1-DECOMP_IN2.
;
; NOTES:
;      Assumes all sorts of things like the same beta, pixellisation: only
;      useful for comparing one object before & after an operation.
;
; MODIFICATION HISTORY:
;      Jul 05 - Polar decomp structures incorporated by RM.
;      Apr 05 - Converted from a function to a procedure by RM.
;      Apr 02 - Writen by Richard Massey.
;-

; Check that the structure types (polar or Cartesian) match
if not tag_exist(decomp_in1,"polar") then decomp_in1=create_struct(decomp_in1,"polar",0B)
if not tag_exist(decomp_in2,"polar") then decomp_in2=create_struct(decomp_in2,"polar",0B)
if ((decomp_in1.polar+decomp_in2.polar) mod 2 ne 0) then $
  shapelets_polar_convert, decomp_in2, polar=decomp_in1.polar, cartesian=1-decomp_in1.polar,/SILENT

; Copy decomp structure
decomp_out=decomp_in1

; Subtract the other's shapelet coefficients (the model is linear)
decomp_out.coeffs=decomp_in1.coeffs-decomp_in2.coeffs

; Do the same to the (statistically independent) errors on the coefficients
decomp_out.coeffs_error=decomp_out.coeffs*sqrt((decomp_in1.coeffs_error/decomp_in1.coeffs)^2+(decomp_in2.coeffs_error/decomp_in2.coeffs)^2)

; Add operation to object's history record
if not keyword_set(nohistory) then shapelets_update_history,decomp_out,decomp_in2.name+" subtracted"

; Return answer
decomp_in1=decomp_out

end

