pro shapelets_polar_convert, input,              $
                             POLAR=polar,        $
                             CARTESIAN=cartesian,$
                             C2P=c2p,            $
                             P2C=p2c,            $
			     PROPAGATE=propagate,$
                             REDUCE=reduce,      $ 
                             EXPAND=expand,      $
                             SILENT=silent,      $
                             MATRIX=matrix

;$Id: shapelets_polar_convert.pro, v2$
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
;      SHAPELETS_POLAR_CONVERT
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Converts Cartesian shapelet coefficients to their equivalent polar
;      shapelet representation, or vice versa.
;
; INPUTS:
;      INPUT - Either a shapelets shapecat or decomp structure,...
;              ...or just the array of coefficients as a floating point/
;              complex vector.
;
; KEYWORD PARAMETERS:
;      POLAR     - Convert into polar shapelet representation.
;      C2P       - A synonym for POLAR.
;      CARTESIAN - Convert into Cartesian shapelet representation.
;      P2C       - A synonym for CARTESIAN.
;      REDUCE    - Removes degenerate coefficients from a polar shapelet
;                  representation to keep the number of coefficients
;                  the same as in the Cartesian case, using
;                  shapelets_polar_reduce.pro. OBSOLESCENT
;      EXPAND    - Recovers the degenerate coefficients removed by REDUCE.
;      SILENT    - Operates silently.
;      PROPAGATE - Calculate errors correctly but irrevesibly. See note below.
;
; OUTPUTS:
;      Either the converted shapelets shapecat or decomp structure,...
;      ...or the converted array of coefficients.
;
; OPTIONAL OUTPUTS:
;      MATRIX    - Optionally returns conversion matrix to speed up multiple 
;                  conversions. This matrix is independent of beta etc, so once
;                  it has been calculated within the current scope, simply use:
;                  polar_coeffs=matrix#cartesian_coeffs
;
; NOTES:
;      A vector of polar shapelet coefficients may be indexed using:
;        shapelets_make_nvec,n_max,n_r,n_l,n_coeffs
;        n=n_r+n_l
;        m=n_r-n_l
;      or
;        shapelets_make_nvec,n_max,n,m,n_coeffs,/POLAR
;
; KNOWN BUGS:
;      The conversion of errors on coefficients is incorrect by default.
;      -----------------------------------------------------------------
;      For example, polar coeffs with m ne 0 involve differences of Cartesian
;      coefficients, say p=c1-c2. The error on p should be the errors on c1 and
;      c2 added in quadrature, but here we calculate the difference of the 
;      errors, thus underestimating the true error. This is done to keep the 
;      conversion reversible. The correct operation can be performed by setting
;      the PROPAGATE flag. However, this involves multiplication by a matrix
;      which is singular matrix if n_max>3, and which cannot be inverted later.
;      It is probably easiest to just use Cartesian shapelets to form errors.
;
; MODIFICATION HISTORY:
;      Jul 05 - PROPAGATE keyword added by RM.
;      Apr 05 - Combined into one general routine by RM.
;      Aug 04 - Shapecat structure conversion code written by RM.
;      Feb 02 - Decomp structure conversion code written by Richard Massey.
;-

COMPILE_OPT idl2

if keyword_set(propagate) then message,"PROPAGATE option not yet coded up!"

if shapelets_structure_type(input,message=message) then begin

  ; Maintain backwards compatibility
  if not tag_exist(input,"polar") then input=create_struct(input,"polar",0B)

  ; Decide which way to convert
  c2p=byte(keyword_set(c2p) or keyword_set(polar))
  p2c=byte(keyword_set(p2c) or keyword_set(cartesian))
  if (c2p+p2c) eq 0B then begin
    p2c=input.polar
    c2p=1B-p2c
  endif
 
  ; Convert from Cartesian to polar shapelet form
  if c2p then begin
  
    if input.polar then begin
      message,"Input "+input.type+" structure is already in polar shapelet form!",/info,noprint=silent
      output=input
    endif else begin
    
      message,"Converting "+input.type+" structure to polar shapelet form!",/info,noprint=silent
      if not keyword_set(reduce) then begin
        ; Calculate conversion matrix and new coefficients
        case strupcase(input.type) of
          "DECOMP":   begin 
                        matrix=shapelets_polar_matrix(input.n_max,/c2p)
                        polar_coeffs=matrix#input.coeffs
			;if keyword_set(propagate) then begin
			;  matrix_error=complex(abs(float(matrix)),abs(imaginary(matrix)))
                        ;  polar_coeffs_error=matrix_error#input.coeffs_error
		        ;endif else begin
                          polar_coeffs_error=matrix#input.coeffs_error
		        ;endelse
                     end
          "SHAPECAT": begin
                        matrix=shapelets_polar_matrix(input.maxn_max,/c2p)
                        polar_coeffs=transpose(matrix#transpose(input.coeffs))
                        ;if keyword_set(propagate) then begin
			;  matrix_error=complex(abs(float(matrix)),abs(imaginary(matrix)))
                        ;  polar_coeffs_error=transpose(matrix_error#transpose(input.coeffs_error))
			;endif else begin
                          polar_coeffs_error=transpose(matrix#transpose(input.coeffs_error))
			;endelse
                      end
          else: message,"Cannot convert a "+input.type+" structure from Cartesian into polar shapelets!"
        endcase
        
        ; Make a new structure, containing the complex coefficients.
        input_tagnames=tag_names(input)
        output={name:input.name,type:input.type}
        for i=0,n_tags(input)-1 do begin
          case strupcase(input_tagnames[i]) of
            "NAME":
            "TYPE":
            "N1": begin
	            output=create_struct(output,"n",input.n1+input.n2)
                    output=create_struct(output,"m",input.n1-input.n2)
                    output=create_struct(output,"nl",input.n2)
		    output=create_struct(output,"nr",input.n1)
                  end
            "N2": 
            "COEFFS": output=create_struct(output,"coeffs",polar_coeffs)
            "COEFFS_ERROR": output=create_struct(output,"coeffs_error",polar_coeffs_error)
            else: output=create_struct(output,input_tagnames[i],input.(i))
          endcase
        endfor
        output.polar=1B 
      endif else begin
        
        ; Convert coefficients from Cartesian to polar shapelet form, but
        ; keep only the independent parameters (as n_coeffs real numbers).
        if strupcase(input.type) ne "SHAPECAT" then message,"The REDUCE option cannot be used with this structure type!"
        message,"The REDUCE option is obsolescent, and should not be used by new code!",/info
        for i=0,input.n-1 do begin
          temp_coeffs=shapelets_c2p(reform(input.coeffs[i,*]),/REDUCE,MATRIX=matrix)
          input.coeffs[i,*]=transpose(temp_coeffs)
          temp_coeffs_error=shapelets_c2p(reform(input.coeffs_error[i,*]),/REDUCE,MATRIX=matrix)
          input.coeffs_error[i,*]=transpose(temp_coeffs_error)
        endfor
        
      endelse
    endelse
  
  endif else begin
  
    if input.polar then begin
      
      ; Calculate conversion matrix and new coefficients
      message,"Converting "+input.type+" structure to Cartesian shapelet form!",/info,noprint=silent
      case strupcase(input.type) of
        "DECOMP":   begin 
                      matrix=shapelets_polar_matrix(input.n_max,/p2c)
                      cartesian_coeffs=float(matrix#input.coeffs)
                      ;if keyword_set(propagate) then begin
		      ;  matrix_error=complex(abs(float(matrix)),-abs(imaginary(matrix)))
                      ;  cartesian_coeffs_error=float(matrix_error#input.coeffs_error)
		      ;endif else begin
                        cartesian_coeffs_error=float(matrix#input.coeffs_error)
		      ;endelse
                    end
        "SHAPECAT": begin
                      matrix=shapelets_polar_matrix(input.maxn_max,/p2c)
                      cartesian_coeffs=float(transpose(matrix#transpose(input.coeffs)))
                      ;if keyword_set(propagate) then begin
		      ;  matrix_error=complex(abs(float(matrix)),-abs(imaginary(matrix)))
                      ;  cartesian_coeffs_error=float(transpose(matrix_error#transpose(input.coeffs_error)))
                      ;endif else begin
                        cartesian_coeffs_error=float(transpose(matrix#transpose(input.coeffs_error)))
		      ;endelse
                    end
        else: message,"Cannot convert a "+input.type+" structure from polar into Cartesian shapelets!"
      endcase
      
      ; Make a new structure, containing just the real coefficients.
      input_tagnames=tag_names(input)
      output={name:input.name,type:input.type}
      for i=0,n_tags(input)-1 do begin
        case strupcase(input_tagnames[i]) of
          "NAME":
          "TYPE":
          "N": if (input.type eq "shapecat") then output=create_struct(output,"n",input.n)
          "M":
          "NR": output=create_struct(output,"n2",input.nl)
          "NL": output=create_struct(output,"n1",input.nr)
          "N1": output=create_struct(output,"n1",input.nr)
          "N2": output=create_struct(output,"n2",input.nl)
          "COEFFS": output=create_struct(output,"coeffs",cartesian_coeffs)
          "COEFFS_ERROR": output=create_struct(output,"coeffs_error",cartesian_coeffs_error)
          else: output=create_struct(output,input_tagnames[i],input.(i))
        endcase
      endfor
      output.polar=0B 
      
    endif else begin
      message,"Input "+input.type+" structure is already in Cartesian shapelet form!",/info,noprint=silent
      output=input
    endelse
  
  endelse


endif else begin

  ; Work on just a coefficent list. This isn't very convenient, as we don't
  ; know things like n_max or the original set of basis functions a priori.
  ; A lot of assumptions are needed in this case!
  n_dimensions=size(input,/N_DIMENSIONS)
  type=size(input,/TYPE)
  if n_dimensions eq 1 and (type eq 4 or type eq 6) then begin

    ; Check vector is in the correct orientation for matrix operation to work
    ;arrsize=size(cartesian_coeffs) & while arrsize(0) gt 1 do reform(cartesian_coeffs)
    input_coeffs=reform(input)

    ; Re-insert degenerate coefficients
    if keyword_set(expand) then input_coeffs=shapelets_polar_expand(input_coeffs)

    ; Work out what n_max was
    n_coeffs = n_elements(input_coeffs)
    n_max=fix(0.5*(sqrt(8.*n_coeffs+1.)-1.))-1
    
    ; Quick check to see if everything is going to plan
    shapelets_make_nvec, n_max, n1, n2, n_n
    if n_n ne n_coeffs then message,"The coefficient list has an unrecognised number of entries!"
    
    ; Convert coefficients
    if keyword_set(p2c) or keyword_set(cartesian) then begin
      matrix=shapelets_polar_matrix(n_max,/p2c)
      output=float(matrix#input_coeffs)
    endif else begin
      matrix=shapelets_polar_matrix(n_max,/c2p)
      output=matrix#input_coeffs
    endelse
        
    ; Remove degenerate coefficients
    if keyword_set(reduce) then output=shapelets_polar_reduce(output)
  
  endif else message,"Cannot parse input variable type!"
endelse

; Replace the old coefficients by the new ones
input=output

end
