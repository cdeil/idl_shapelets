function shapelets_phi, n, x1, x2,           $
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
;      SHAPELETS_PHI
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Compute (dimensionless) Cartesian shapelet basis functions phi(x)
;      in 1D or 2D, based on Hermite polynomials. Dimensionful basis functions
;      can be calculated with shapelets_phi(n,x1/beta[,x2/beta])/beta.
;
; INPUTS:
;      n - Basis function order. Integer n OR vector [n1,n2].
;      If you supply an n with 1 dimension, it will return a 1D basis function.
;      If you supply an n with 2 dimensions, it will return a 2D basis function.
;
; OPTIONAL INPUTS:
;      x1   - Grid of x coordinates (e.g. created by shapelets_make_xarr.pro).
;      x2   - Grid of y coordinates if a 2D basis function is required.
;      BETA - Shapelet scale size beta.
;
; KEYWORD PARAMETERS:
;      /INTEGRATE - Default behaviour is to simply evaluate phi at the centre
;                   of each pixel. If /INTEGRATE is set, the routine will use
;                   the recursion relation in Shapelets III to integrate the
;                   the basis functions within pixels (currently 2D only).
;      /ARRAY     - Instead of just returning one basis function, if /AARRAY is
;                   set, the routine will return an array containing all of
;                   the (pixellated) basis functions up to n_max=n[0]+n[1].
;
; OUTPUTS:
;      Cartesian shapelet basis function phi_n1[_n2](x[,y])
;
; EXAMPLE USE:
;      plot,x,shapelets_phi(4,x,beta=1),psym=-3
;      tvscl,float(shapelets_phi([2,2],/integrate))
;
; TO DO:
;      Truncate pixellation to within pixels near the centre only.
;      Model unresponsive dead zones around the edges of CCD pixels by 
;       rescaling pix_size.
;
; MODIFICATION HISTORY:
;      Jul 05 - BETA input added by RM to make dimensionful basis functions.
;      Feb 05 - Bug fixed when using /ARRAY but not /INTEGRATE by RM.
;      Aug 04 - Error catching for 1D arrays of 2D functions improved by RM.
;      Oct 03 - Speed up by efficient calculation of polynomials by RM.
;      Oct 03 - Bug fixed in integration limits by Alain Bonissent.
;      Sep 03 - Analytic integration option for 1D shapelets added by RM.
;      Nov 01 - Analytic integration option for 2D shapelets added by RM.
;      Sep 01 - 1D/2D version written by Richard Massey.
;      Jul 99 - Written by Alexandre Refregier.
;-


if min(n) lt 0 then message,'n must be non-negative!'
if not keyword_set(beta) then beta=1.

; Decide whether to calcuate 1D or 2D basis functions.
case n_elements(n) of

  1: begin

      if not keyword_set(x1) then x1=(findgen(401)/20.-10.)*beta ; Default x array.
      n_pix=(size(x1,/DIMENSIONS))[0]

      if keyword_set(integrate) then begin

        ; Find the extremeties of each cell, to be the limits of the integrals.
        if n_pix eq 1 then message,"Can't really integrate within one point..."
        pix_size= x1[1]-x1[0]           ; Width of individual grid cells.
        left    = (x1-pix_size/2.)/beta ; Half a pixel width in each direction
        right   = (x1+pix_size/2.)/beta ;  gets us to the edges of the pixel.
        roottwo = sqrt(2.)
        piquart = !pi^0.25
        
        ; Evaluate integrals.
        Integral  = fltarr(n_pix,n+1,/nozero)
        ; I0...
        Integral[*,0] = piquart/roottwo * $
                        (erf(right[*]/roottwo)-erf(left[*]/roottwo))

        ; I1...
        if n ge 1 then begin
          gaussianl=exp(-.5*(left^2))
          gaussianr=exp(-.5*(right^2))
          Integral[*,1] = ( -1. * roottwo / piquart ) * $
                          ( shapelets_phi(0,right) - shapelets_phi(0,left) )

    	  ; I2 and above...
    	  i=2
    	  while i le n do begin
    	    Integral[*,i] = sqrt((float(i)-1.)/float(i)) * Integral[i-2,*] $
    			   + (-1. * sqrt(2./float(i))) * ( shapelets_phi(i-1,right) - shapelets_phi(i-1,left) )
    	    i=i+1
    	  endwhile
        endif

        ; Put these integrals into basis functions.
        if keyword_set(array) then BasisF = Integral / grid_width $
          else BasisF = reform( Integral[*,n] ) / grid_width

      endif else begin

        ; Simply calculate basis functions phin_n(x) at the centres of bins.

        if keyword_set(array) then begin
          ; Calculate all the basis functions up to n_max.
          BasisF=fltarr(n_pix, n+1, /nozero)
          gaussian=exp(-.5*((x1/beta)^2))/piquart/sqrt(beta)
          for i=0,n do BasisF[*,i]=shapelets_hermite(i,(x1/beta))*gaussian/$
             sqrt(2.^i*factorial(i))
        endif else begin
          ; Calculate just that one basis function.
          constfact=1/sqrt((2.^n[0]*sqrt(!pi)*factorial(n[0]))*beta)
          BasisF=shapelets_hermite(n[0],(x1/beta))*exp(-.5*((x1/beta)^2))*constfact
        endelse

      endelse

    end

 2: begin

  ; Calculate 2D basis functions.

  ; Basic groundwork.
  n_max=n[0]+n[1]
  shapelets_make_nvec,n_max,n1,n2,n_a
  
  ; Default x and y arrays.
  if not keyword_set(x1) or not keyword_set(x2) then begin
      shapelets_make_xarr,[128,128],x1,x2
    x1=x1*beta/8
    x2=x2*beta/8
  endif

  ; Calculate basis functions
  if keyword_set(integrate)         and $
     size(x1,/N_DIMENSIONS) eq 2    and $  ; Error catching for 1xn arrays
     min(size(x1,/DIMENSIONS)) gt 1 then begin

    ; Assuming the grid is evenly spaced and aligned to the x- and y-axes
    ; lets us calculate the (separable) basis functions only at the edges
    ; then quickly replicate the values to fill the entire array.

    ; Find the extremeties of each cell, to be the limits of the integrals.
    n_pix=size(x1,/DIMENSIONS)
    n_pix_x=n_pix[0]
    n_pix_y=n_pix[1]
    pix_size_x = x1[1,0]-x1[0,0] ; Width of individual grid cells
    pix_size_y = x2[0,1]-x2[0,0] ; Height of individual grid cells
    if keyword_set(deadzone) then message,"Dead zones in pixels not yet implemented!"
    left   = (x1[*,0]-0.5*pix_size_x)/beta ; Half a pixel width in each direction 
    right  = (x1[*,0]+0.5*pix_size_x)/beta ;  gets us to the edges of the pixel.
    bottom = (x2[0,*]-0.5*pix_size_y)/beta ;
    top    = (x2[0,*]+0.5*pix_size_y)/beta ;

    ; Evaluate integrals.
    roottwo    = sqrt(2.)
    piquart    = !pi^0.25
    Integral_x = fltarr(n_pix_x,n_max+1,/nozero)
    Integral_y = fltarr(n_pix_y,n_max+1,/nozero)

    ; I0...
    coeff1 = piquart / roottwo
    Integral_x[*,0]=coeff1*(erf(right/roottwo)-erf(left/roottwo))
    Integral_y[*,0]=coeff1*(erf(top/roottwo)-erf(bottom/roottwo))

    ; I1...
    if n_max ge 1 then begin
      gaussianl=exp(-.5*(left^2))
      gaussianr=exp(-.5*(right^2))
      gaussianb=exp(-.5*(bottom^2))
      gaussiant=exp(-.5*(top^2))
      coeff2 = -1 * roottwo / piquart  ; coeff1=/=0 but I(-1) not defined...
      Integral_x[*,1] = coeff2 * ( shapelets_hermite(0,right)*gaussianr - shapelets_hermite(0,left)*gaussianl )
      Integral_y[*,1] = coeff2 * ( shapelets_hermite(0,top)*gaussiant - shapelets_hermite(0,bottom)*gaussianb )

      ; I2 and above...
      i=2
      while i le n_max do begin
        coeff1 = sqrt((float(i)-1.)/float(i))
        coeff2 = -1. * sqrt(2./float(i)) / sqrt(2.^(i-1)*sqrt(!pi)*factorial(i-1))
        Integral_x[*,i] = coeff1 * Integral_x[*,i-2] $
          + coeff2 * ( gaussianr*shapelets_hermite(i-1,right) $
          - gaussianl*shapelets_hermite(i-1,left)   )
        Integral_y[*,i] = coeff1 * Integral_y[*,i-2] $
          + coeff2 * ( shapelets_hermite(i-1,top)*gaussiant   $
          - shapelets_hermite(i-1,bottom)*gaussianb )
        i=i+1
      endwhile
    endif

    ; Cope with limits when beta wasn't known.
    ; (essentially multiplies by beta, but we divide again outside this routine)
    Integral_x=Integral_x/pix_size_x*sqrt(beta)
    Integral_y=Integral_y/pix_size_y*sqrt(beta)


    ; Put these integrals into basis functions.
    if keyword_set(array) then begin
      BasisF = fltarr(n_pix_x,n_pix_y,n_a,/nozero)
      i=0
      for j=0,n_max do begin
        for n1i=0,j do begin
          n2i=j-n1i
          BasisF[*,*,i] = Integral_x[*,n1i] # replicate(1.,n_pix_y) * $
                          replicate(1.,n_pix_x) # Integral_y[*,n2i]
          i=i+1
        endfor
      endfor
    endif else begin
      BasisF = Integral_x[*,n[0]] # replicate(1.,n_pix_y) * $
               replicate(1.,n_pix_x) # Integral_y[*,n[1]]
    endelse

  endif else begin

    ; Simply evaluate phi_n(x,y) at the grid points (usually the centres 
    ; of pixels). This will also be useful if we have an irregular grid, 
    ; so I won't make the shortcut of assuming that we can just calculate
    ; polynomials along the sides of the array, then replicate them across,
    ; as we do for the integration-within-pixels case.

    ; Error catching
    if keyword_set(integrate) then $
      message,"Can't integrate within a one-pixel array!",/info;,noprint=silent

    ; Get number of pixels
    n_pix=size(x1,/DIMENSIONS)
    n_pix_x=n_pix[0]
    if size(x1,/N_DIMENSIONS) gt 1 then n_pix_y=n_pix[1]

    if keyword_set(array) then begin
      ; Calculate all the basis functions up to nmax.
      BasisF=fltarr(n_pix_x, n_pix_y, n_a, /nozero)
      for i=0,n_a-1 do begin
        gaussian=exp(-.5*(x1^2+x2^2)/beta^2)/sqrt(!pi)/beta
        BasisF[*,*,i]=shapelets_hermite(n1[i],(x1/beta)) * $ ; This could be stored for speed
                      shapelets_hermite(n2[i],(x2/beta)) * $ ; This could be stored for speed
                      gaussian /                           $
                      sqrt( 2.^(n1[i]+n2[i]) *             $
                      factorial(n1[i]) *                   $
                      factorial(n2[i]) )
      endfor
    endif else begin
      ; Calculate just that one basis function.
      BasisF = shapelets_hermite(n[0],(x1/beta)) * $
               shapelets_hermite(n[1],(x2/beta)) * $
               exp(-.5*(x1^2+x2^2)/beta^2)       / $
               beta                              / $
               sqrt( 2.^(n[0]+n[1]) * !pi * factorial(n[0]) * factorial(n[1]) )
    endelse

  endelse

 end
 else: message,"Currently works only for 1D or 2D basis functions!"
endcase

return, BasisF

end
