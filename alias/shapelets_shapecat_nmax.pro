pro shapelets_shapecat_nmax, shapecat, n_max, SILENT=silent

; NAME:
;      SHAPELETS_SHAPECAT_NMAX
;
; PURPOSE:
;      Alter n_max in a shapelet catalogue.
;
;
; CATEGORY:
;       Shapelet catalogue manipulation.
;
; INPUTS:
;       shapecat - IDL shapecat structure.
;       n_max    - Desired new n_max (can be the same as, less than or
;                  greater than that already used in shapecat).
;
; OUTPUTS:
;       shapecat - IDL shapecat structure.
;
; MODIFICATION HISTORY:
;       Apr 05 - Obsoleted by RM
;       Aug 04 - Rendered compatabile with polar shapelet catalogues by RM.
;       May 04 - Written by Richard Massey

COMPILE_OPT idl2, OBSOLETE

; Do we need to change anything after all?
if (n_max ne shapecat.maxn_max) then begin

  ; Determine what we want n_max to be
  shapelets_make_nvec, n_max, n1, n2, n_coeffs
  n=n1+n2
  
  ; Truncate coefficients to reduce n_max
  if n_max lt shapecat.maxn_max then begin
    message,'Truncating coefficients above '+$
            'n_max='+strtrim(string(n_max),2),/info,noprint=silent
    keep=where(n le n_max,n_coeffs)
    coeffs=shapecat.coeffs[*,keep]
    coeffs_error=shapecat.coeffs_error[*,keep]
  endif else begin
  
  ; Pad with zeros to increase n_max
    message,'Padding coefficients with zeros to increase from '+$
            'n_max='+strtrim(string(shapecat.maxn_max),2)+$
            ' to n_max='+strtrim(string(n_max),2),/info,noprint=silent
    if shapecat.polar then begin
      coeffs=[[shapecat.coeffs],[complexarr(shapecat.n,n_coeffs-shapecat.n_coeffs)]]
      coeffs_error=[[shapecat.coeffs_error],[complexarr(shapecat.n,n_coeffs-shapecat.n_coeffs)]]
    endif else begin
      coeffs=[[shapecat.coeffs],[fltarr(shapecat.n,n_coeffs-shapecat.n_coeffs)]]
      coeffs_error=[[shapecat.coeffs_error],[fltarr(shapecat.n,n_coeffs-shapecat.n_coeffs)]]
    endelse
  endelse
  
  ; Adjust shapecat struture to reflect new length of coefficient list
  shapecat_tagnames=tag_names(shapecat)
  temp_shapecat={name:shapecat.name,type:shapecat.type}
  for i=0,n_tags(shapecat)-1 do begin
    case strupcase(shapecat_tagnames[i]) of
      "NAME":
      "TYPE":
      "COEFFS": temp_shapecat=create_struct(temp_shapecat,"coeffs",coeffs)
      "COEFFS_ERROR": temp_shapecat=create_struct(temp_shapecat,"coeffs_error",coeffs_error)
      else: temp_shapecat=create_struct(temp_shapecat,shapecat_tagnames[i],shapecat.(i))
    endcase
  endfor
  shapecat=temp_shapecat
  shapecat.maxn_max=n_max
  shapecat.n_coeffs=n_coeffs
  
endif

end
