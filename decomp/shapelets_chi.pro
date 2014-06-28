function shapelets_chi, nm, x1, x2,          $
                        INTEGRATE=integrate, $
                        ARRAY=array,         $
                        BETA=beta,           $
                        DEADZONE=deadzone

;$Id: shapelets_phi.pro, v2$
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
;      SHAPELTS_CHI
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Compute (dimensionless) 2D polar shapelet basis functions chi_n_m(x,y),
;      based on Laguerre polynomials. Dimensional basis functions
;      can be calculated with shapelets_chi(nm,x1/beta[,x2/beta])/beta.
;
; INPUTS:
;      nm - Basis function order. Integer vector [n,m].
;
; OPTIONAL INPUTS:
;      x1   - Grid of x coordinates (e.g. created by shapelets_make_xarr.pro).
;      x2   - Grid of y coordinates.
;      BETA - Shapalet scale size beta.
;
; KEYWORD PARAMETERS:
;      /INTEGRATE - Default behaviour is to simply evaluate phi at the centre
;                   of each pixel. If /INTEGRATE is set, the routine will use
;                   the recursion relation in Shapelets III to integrate the
;                   the Cartesian basis functions within pixels, then convert
;                   these to polar shapelets via the matrix in Shapelets I.
;      /ARRAY     - Instead of just returning one basis function, if /AARRAY is
;                   set, the routine will return an array containing all of
;                   the (pixellated) basis functions up to n_max=n.
;
; OUTPUTS:
;      Polar shapelet basis function chi_n1[_n2](x[,y]), in a complex array.
;
; EXAMPLE USE:
;      plot,r,shapelets_chi(4,r,beta=1),psym=-3
;      tvscl,float(shapelets_chi([2,2],/integrate))
;
; MODIFICATION HISTORY:
;      Jul 05 - 1D basis functions enabled by RM.
;      Jul 05 - Minus sign in definition of polar basis functions fixed by RM.
;      Jul 05 - BETA input added by RM to make dimensionful basis functions.
;      Oct 03 - Written by R. Massey.
;-

COMPILE_OPT idl2

if nm[0] lt 0 then message,"n must be non-negative!"
if not keyword_set(beta) then beta=1.

; Decide whether to calcuate 1D or 2D basis functions.
case n_elements(nm) of

 1: begin

    ; Calculate 1D basis functions (the radial part of m=0 polar shapelets).
    if keyword_set(integrate) then message,"WARNING: polar integration not yet coded. Evaluating at centres of pixels!",/info
    if (nm[0] mod 2) then message,"n must be even!"
    
    ; Form default coordinate grid
    if not keyword_set(x1) then x1=(findgen(201)/20.)*beta ; Default x array.
    
    ; Determine size of grid
    n_pix=(size(x1,/DIMENSIONS))[0]

    if keyword_set(array) then begin
      ; Calculate all the basis functions up to n
      BasisF=fltarr(n_pix, (nm[0]/2)+1, /nozero)
      gaussian=exp(-.5*(x1/beta)^2)/(beta*sqrt(!pi))
      for i=0,nm[0]/2 do BasisF[*,i]=(-1)^i*gaussian*laguerre((x1/beta)^2,i)
    endif else begin
      ; Calculate just that one basis function.
      BasisF=(-1)^(nm[0]/2)/sqrt(!pi)/beta*$
   		  laguerre((x1/beta)^2,nm[0]/2)*exp(-.5*(x1/beta)^2)
    endelse
    
    
    end

 2: begin

    ; Calculate 2D basis functions.
    if ((nm[0]+nm[1]) mod 2) then message,"n and m must both be odd or both be even!"

    ; Form default Cartesian coordinate grid
    if not keyword_set(x1) or not keyword_set(x2) then begin
      shapelets_make_xarr,[128,128],x1,x2
      x1=x1*beta/8.
      x2=x2*beta/8.
    endif
    
    ; Determine size of grid
    n_pix_x=(size(x1,/DIMENSIONS))[0]
    n_pix_y=(size(x1,/DIMENSIONS))[1]
    
    if not keyword_set(integrate) then begin
      
      ; Prepare vectors of n, m, n_l and n_r
      if not keyword_set(array) then begin
  	n_coeffs=1
  	n=nm[0]
  	m=nm[1]
      endif else shapelets_make_nvec,nm[0],n,m,n_coeffs,/POLAR

      ; Evaluate only those pixels that are likely to be nonzero, by looking
      ;   at only those within 10 th_max of the origin
      BasisF = complexarr(n_pix_x,n_pix_y,n_elements(n))
      close  = where((x1^2)>(x2^2) lt 10*beta^2*(n+1),n_close)
      if n_close gt 0 then begin
  	close_l= min(close mod n_pix_x,max=close_r)
  	close_b= min(close/n_pix_x,max=close_t)
    
  	; Form polar coordinate grid.
  	rsq    = (x1[close_l:close_r,close_b:close_t]^2+x2[close_l:close_r,close_b:close_t]^2)/beta^2
  	theta  = atan(x2[close_l:close_r,close_b:close_t],x1[close_l:close_r,close_b:close_t])
  	rootpi = sqrt(!pi)
    
  	; Draw basis functions
  	const  = ((-1)^((n-abs(m))/2)) * sqrt(factorial((n-abs(m))/2)/factorial((n+abs(m))/2)) / beta / rootpi
  	gauss  = exp(-rsq/2.)
  	for i=0,n_coeffs-1 do BasisF[close_l:close_r,close_b:close_t,i] =  $
  			      (const[i] * rsq^(abs(m[i]/2.)) )  	 * $
  			      laguerre(rsq,(n[i]-abs(m[i]))/2,abs(m[i])) * $
  			      gauss					 * $
  			      exp(complex(0,-m[i]*(theta)))
      endif
    
    endif else begin

      ; Integrate phi_n(x,y) within pixels...
      BasisC=shapelets_phi([nm[0],0],x1,x2,/ARRAY,/INTEGRATE,BETA=beta,DEADZONE=deadzone)
      shapelets_make_nvec,nm[0],n,m,n_coeffs,/POLAR
      BasisC=reform(BasisC,n_pix_x*n_pix_y,n_coeffs)

      ; ...then convert to chi_nm(x,y)
      ; THIS IS SLOW!
      matrix=shapelets_polar_matrix(nm[0],/C2P)
      if not keyword_set(array) then begin
  	sel=where(n eq nm[0] and m eq nm[1])
  	BasisF=reform(matrix[sel,*]) ## BasisC 
  	BasisF=reform(BasisF,n_pix_x,n_pix_y)
      endif else begin
  	BasisF=transpose(matrix) ## BasisC
  	BasisF=reform(BasisF,n_pix_x,n_pix_y,n_coeffs)
      endelse

    endelse
    end
else: message,"Currently works only for 1D or 2D basis functions!"
endcase

return, BasisF

end


