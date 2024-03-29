function shapelets_decomp, IMAGE, BETA, N_MAX,    $
                           CENTRE=centre,	        $
                           NAME=name,		          $
                           PSF=psf,		            $
                           NOISE=noise, 	        $
                           RECOMP=recomp,         $
                           OVERSAMPLE=oversample, $
                           OVERLAP=overlap,	      $
                           INTEGRATE=integrate,   $
                           SKY=sky,		            $
                           NON1=non1,		          $
                           POLAR=polar, 	        $
                           DIAMOND=diamond,	      $
                           FULL_ERROR=full_error, $
                           SILENT=silent,	        $
                           LS=ls,		              $
                           X0=x0

;$Id: shapelets_decomp.pro, v2$
;
; Copyright � 2005 Richard Massey and Alexandre Refregier.
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
;      SHAPELETS_DECOMP
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Decompose an image into a weighted sum of shapelet basis functions.
;
; INPUTS:
;      IMAGE     - Shapelets image/pstamp structure, or 2D image array.
;      BETA      - Shapelet basis function scale size.
;      N_MAX     - Maximum order of the basis coefficients used in model.
;
; OPTIONAL INPUTS:
;      PSF       - A decomp structure representing the local PSF, to be
;                  decomvolved.
;      CENTRE    - Coordinates of the centre of the basis functions.
;                  DEFAULT: from structure, or centre of the image.
;      NAME      - Name of object (string).
;                  DEFAULT: from structure, or empty.
;      NOISE     - Inverse variance map of pixels, or constant noise value.
;                  Assumes zero covariance between adjacent pixels.
;                  DEFAULT: from structure, or constant unity.
;      OVERSAMPLE- Force an oversampling factor to evaluate the basis fns.
;                  the basis functions (over=1: no oversampling;
;                  DEFAULT: over=round(1/theta_min)
;      SKY       - 1: fit sky background with a constant value
;                  2: fit sky gradient with a plane
;            	   DEFAULT: no sky subtraction
;	      
; KEYWORD PARAMETERS:
;      POLAR     - Perform decomposition using polar shapelet basis fns.
;                  DEFAULT: Cartesian shapelet basis functions.
;      OVERLAP   - Obtain coefficients via Fourier-style overlap integrals.
;                  DEFAULT: least-squares fitting.
;      FULL_ERROR- Return full covariance matrix of coefficients.
;            	   DEFAULT: return error on coefficients in error variable.
;      INTEGRATE - Use basis functions integrated within pixels.
;            	   DEFAULT: integration (for central pixel value, set int=0)
;      NON1      - n=1 coefficients are forced to be zero, Kuijken-like.
;      DIAMOND   - Alternate truncation scheme for shapelet basis functions.
;                  DEFAULT: usual tringular scheme, with n<=n_max
;
; OUTPUTS:
;      DECOMP    - A shapelet decomp structure is returned, containing the
;                  evaluated shapelet coefficeints, meta-parameters and other
;                  information.
;
; OPTIONAL OUTPUTS
;      RECOMP    - Repixellated image of the shapelet model (still containing
;                  the sky background and convolved with the original PSF. To
;                  remove this, use shapelets_recomp,decomp,recomp).
;
; CALLING SEQUENCE:
;      decomp=shapelets_decomp(image, beta, n_max)
;
; KNOWN BUGS:
;      Overlap integrals give slightly different results for polar and
;      Cartesian shapelet basis functions. Results are identical when fitting,
;      and I'd've thought that they would also be for overlap. Should they?
;
; MODIFICATION HISTORY:
;      Sep 05 - Changed from a procedure to a function by RM.
;      Sep 05 - POLAR and DIAMOND options added by RM.
;      Aug 05 - Images accepted in a wider variety of input formats by RM.
;      Aug 05 - Creation of LS fitting matrix internalised and whole routine
;               reorganised in an attempt to speed things up by RM.
;      Jul 05 - External decomp creation routine incorporated by RM.
;      Jan 05 - RM optimised gamma and n_max_gamma during PSF deconvolution.
;      Apr 04 - RM added ability to mask pixels (by setting weight map to zero)
;               but still calculate the correct number of DOFs in the fit
;      Sep 03 - RM changed name of routine to shapelet_decomp.pro
;      Jul 03 - RM added "name" element to decomp structure
;      Apr 02 - RM added PSF deconvolution
;      Mar 02 - AR computed Chi^2 and fixed normalisation. Also added
;               option to return the recomposed (still convolved) image
;      Dec 01 - RM implemented least squares fitting of basis fns to data
;      Nov 01 - RM added background subtraction and sky fitting
;      Nov 01 - Richard Massey added numerical integration of basis fns
;      Feb 01 - AR incorporated oversampling of the basis functions
;      Jul 00 - AR changed to order by n=n1+n2
;      Jul 99 - mk_decomp.pro written by Alexandre Refregier
;-

COMPILE_OPT idl2
ON_ERROR,2


;
; Parse inputs and set default parameters
;
if n_params() ne 3 or size(n_max,/TYPE) ge 4 then message,"Usage: decomp=shapelets_decomp(image,beta,n_max)"
if keyword_set(ls) then message,"Least-squares fitting now the default. Use /OVERLAP for linear overlap method."
if keyword_set(x0) then message,"Keyword X0 is obsolete. Please use CENTRE.",/INFO
if keyword_set(centre) then x0=centre           ; Centre of the basis functions
if keyword_set(overlap) then ls=0B else ls=1B   ; Default is to fit rather than do overlap integrals
if n_elements(integrate) eq 0 then integrate=1B ; Interate basis functions within pixels
if n_max le 1 then non1=0B                      ; Can't avoid n=1 coeffs if there are no higher ones available!
if keyword_set(oversample) then begin           ; Decide pixel subsampling rate
  if oversample eq 1 then over=sqrt(n_max)/float(beta) else over=oversample
  over=round(over)>1
endif else over=1
if shapelets_structure_type(image,message=message,/SILENT) then begin
  if strupcase(image.type) eq "IMAGE" or strupcase(image.type) eq "PSTAMP" then begin
    data=image.image
    if not keyword_set(name) then name=image.name
    if not keyword_set(noise) then if tag_exist(image,"noise") then noise=image.noise
    if keyword_set(noise) then if n_elements(noise) gt 1 then if tag_exist(image,"mask") then begin
      masked_pixels=where(image.mask eq 1,n_masked_pixels)
      if n_masked_pixels gt 0 then noise[masked_pixels]=0
    endif
    if not keyword_set(x0) then if tag_exist(image,"xo") and tag_exist(image,"yo") then begin
      x0=[image.xo,image.yo]
    endif
  endif else message,"Cannot apply shapelet transform to a "+image.type+" structure!"
endif else if size(image,/n_dimensions) eq 2 then begin
  data=image
endif else if size(image,/n_dimensions) eq 1 then begin
  message,"One dimensional shapelet transform not yet implemented"
endif else message,"Input image format not recognised!"
if keyword_set(psf) then begin
  if not shapelets_structure_type(psf,message=message,/SILENT) then message,message
  shapelets_polar_convert,psf,/P2C,/SILENT
endif
if not keyword_set(name) then name=""

;
; Reformat data and pixel weight (inverse variance) maps into a vector
;
fsize       = size(data)
n_pixels_x  = fsize[1]
n_pixels_y  = fsize[2]
n_pixels    = n_pixels_x*n_pixels_y
n_pixels_x_o= n_pixels_x*over
n_pixels_y_o= n_pixels_y*over
n_pixels_o  = n_pixels_x_o*n_pixels_y_o
data_o      = reform(rebin(data,n_pixels_x_o,n_pixels_y_o,/sample), n_pixels_o, /overwrite)
if keyword_set(noise) then begin
  if n_elements(noise) gt 1 then begin
    noise_o = reform(rebin(noise>0,n_pixels_x_o,n_pixels_y_o,/sample), n_pixels_o, /overwrite)
  endif
endif else noise_o=1


;
; Finalise shapelet meta-parameters
;
if not keyword_set(x0) then x0=[n_pixels_x,n_pixels_y]/2.
x0_o=x0*over
beta_o=float(beta)*over


;
; Calculate PSF convolution matrix
;
if keyword_set(psf) then begin
  if keyword_set(overlap) then begin
    ; Force values of shapelet meta-parameters so that P_nm is square
    n_max_gamma=n_max
    gamma=beta_o
  endif
  P_nm=shapelets_convolution_matrix(psf,beta_o,n_max,gamma,n_max_gamma)
endif else begin
  n_max_gamma=n_max
  gamma=beta
endelse


;
; Calculate basis functions
;
shapelets_make_xarr, [n_pixels_x_o, n_pixels_y_o], x1, x2, x0=x0_o
Basis=shapelets_phi([n_max_gamma,0],x1,x2,beta=gamma,integrate=integrate,/array)
MatrixT=transpose(reform(Basis,n_pixels_o,(size(Basis,/DIMENSIONS))[2]))


;
; Convolve basis functions with the PSF
;
if keyword_set(psf) and not keyword_set(overlap) then begin
  MatrixT=transpose(P_nm)#MatrixT
endif


;
; Convert complex polar shapelet basis functions into real components
;
if keyword_set(polar) or keyword_set(diamond) then begin
  c2p_matrix=shapelets_polar_matrix(n_max,/C2P)  
  MatrixT=transpose(transpose(c2p_matrix) ## transpose(MatrixT))
  shapelets_make_nvec, n_max, n, m, n_coeffs, /polar
  m_negative=where(m lt 0,n_m_negative)
  if n_m_negative gt 0 then  MatrixT[m_negative,*]=MatrixT[m_negative,*]*complex(0,1)
  MatrixT=float(MatrixT)
endif


;
; Remove basis functions for the diamond truncation scheme
;
if keyword_set(diamond) then begin
  diamond_coeffs=where(n+abs(m) le n_max)
  MatrixT=MatrixT[diamond_coeffs,*]
endif


;
; Add extra basis functions to model the sky background
;
if keyword_set(sky) then begin
  MatrixT=[MatrixT,replicate(1.,1,n_pixels_o)]                        ; Constant
  if sky eq 2 then begin
    MatrixT=[MatrixT,reform(x1,1,n_pixels_o),reform(x2,1,n_pixels_o)] ; Plane
    sky_size=3
  endif else sky_size=1
endif else sky_size=0


;
; Remove some basis functions if we are requiring a_01 and a_10 to be zero
;
if keyword_set(non1) then begin
  MatrixT[2,*]=MatrixT[0,*]
  MatrixT=MatrixT[2:*,*]
endif
; Old version of code to get to (roughly) this point
;shapelets_make_ls_matrix, MatrixT, n_max_gamma, gamma, $
;  [n_pixels_x_o, n_pixels_y_o], x0_o, $
;  integrate=integrate, polar=polar, diamond=diamond, sky=sky, non1=non1


;
; Perform shapelet decomposition
;
Matrix=transpose(MatrixT)
if keyword_set(overlap) then begin
  ; Perform "overlap integral" type of decomposition
  coeffs=MatrixT#data_o
  ; We don't know the errors on coefficients using this method, so set them to zero
  if keyword_set(polar) then begin
    coeffs_error=fltarr(n_elements(coeffs))
  endif else begin
    coeffs_error=complexarr(n_elements(coeffs))
  endelse
endif else begin
  ; Perform least-squares linear algebra fit
  ; (Matrix is transpose(M) of Anton p460+ or M in Lupton p84)
  if keyword_set(noise) then begin
    if n_elements(noise_o) gt 1 then begin
      for i=0L,n_pixels_o-1 do MatrixT[*,i]=MatrixT[*,i]*noise_o[i] ; / since V^-1 / since INVERSE variance
    endif else begin
      MatrixT=temporary(MatrixT)*noise_o
    endelse
  endif
  ; NB: this might be more numerically stable using SVD (cf Berry, Hobson & Withington 2004)
  MatrixI=invert(MatrixT#Matrix) ; size of matrix to invert is only n_coeffs^2
  coeffs=MatrixI#MatrixT#data_o
endelse


;
; Make a recomposed image (but still PSF convolved and including sky fit)
;
if keyword_set(overlap) then begin
  recomp_o=Matrix[*,0:(size(Matrix,/DIMENSIONS))[1]-1-sky_size]#coeffs[0:n_elements(coeffs)-1-sky_size]
endif else begin
  recomp_o=Matrix#coeffs
endelse
recomp=rebin(reform(recomp_o,n_pixels_x_o,n_pixels_y_o),n_pixels_x,n_pixels_y)


;
; Compute Chi^2 residual (directly from the original image and the recomposed
;  model, but exactly the same calculation as that done with matrices in Lupton)
;
chisq=total((data_o-recomp_o)^2*noise_o)
masked_pixels=where(noise_o eq 0.,n_masked_pixels)
dof=n_pixels_o-n_masked_pixels-n_elements(coeffs)
; Warn if the number of coefficients is larger than the number of pixels
if dof le 0 then begin
  message,"WARNING: No of coefficients "+strtrim(n_pixels-dof,1)+$
    " > No of pixels "+strtrim(n_pixels,1),/info,noprint=silent
  ; Could prevent least-squares from overfitting and then failing
  ; by checking this earlier? Then should set n_max to a lower number,
  ; but store n_max and recover it afterwards, setting a(^) to zero.
  chisq=[0.,0.]
endif else begin
  chisq_reduced=chisq/float(dof)
  chisq=[chisq,chisq_reduced]
endelse


;
; Calculate the covariance matrix of coefficients
;
if keyword_set(overlap) then begin
  if keyword_set(cov) then begin
    coeffs_error=fltarr(n_elements(coeffs),n_elements(coeffs))
  endif else begin
    coeffs_error=fltarr(n_elements(coeffs))
  endelse
endif else begin
  if keyword_set(cov) then begin
    coeffs_error=MatrixI
  endif else begin
    coeffs_error=fltarr(n_elements(coeffs))
    for i=0,n_elements(coeffs)-1 do coeffs_error[i]=sqrt(MatrixI[i,i])
  endelse
endelse


;
; Reinsert n=1 coefficients if they had been omitted
;
if keyword_set(non1) then begin
  coeffs = [coeffs[0],0.,0.,coeffs[1:*]]
  coeffs_error = [coeffs_error[0],0.,0.,coeffs_error[1:*]]
endif


;
; Discard sky background fit if it had been done
;
skyfit=fltarr(3)
if keyword_set(sky) then begin
  first_sky_coeff=n_elements(coeffs)-sky_size
  skyfit[0:sky_size-1]=coeffs[first_sky_coeff:*]
  coeffs=coeffs[0:first_sky_coeff-1]
  coeffs_error=coeffs_error[0:first_sky_coeff-1]
  if keyword_set(overlap) then skyfit[0]=skyfit[0]/n_pixels_o
endif


;
; Reinsert coefficients excluded by diamond truncation scheme
;
if keyword_set(diamond) then begin
  temp_coeffs=fltarr(n_coeffs)
  temp_coeffs_error=fltarr(n_coeffs)
  temp_coeffs[where(n+abs(m) le n_max)]=coeffs
  temp_coeffs_error[where(n+abs(m) le n_max)]=coeffs_error
  coeffs=temp_coeffs
  coeffs_error=temp_coeffs_error
endif


;
; Expand polar shapelet coefficients back into complex form
;
if keyword_set(polar) or keyword_set(diamond) then begin
  temp_coeffs=coeffs/2
  temp_coeffs_error=coeffs_error/2
  coeffs=complex(coeffs)
  coeffs_error=complex(coeffs)
  for nn=1,n_max do begin
    for mm=2-(nn mod 2),nn,2 do begin
      m_positive=where(n eq nn and m eq mm)
      m_negative=where(n eq nn and m eq -mm)
      coeffs[m_positive]=complex(temp_coeffs[m_positive],temp_coeffs[m_negative])
      coeffs[m_negative]=complex(temp_coeffs[m_positive],-temp_coeffs[m_negative])
      coeffs_error[m_positive]=complex(temp_coeffs_error[m_positive],temp_coeffs_error[m_negative])
      coeffs_error[m_negative]=complex(temp_coeffs_error[m_positive],-temp_coeffs_error[m_negative])
    endfor
  endfor
endif


;
; Convert back to Cartesian shapelet coefficients
;
if keyword_set(diamond) and not keyword_set(polar) then begin
  p2c_matrix=shapelets_polar_matrix(n_max,/P2C)
  coeffs=p2c_matrix#coeffs
  coeffs_error=p2c_matrix#coeffs_error
endif


;
; Deconvolve from the PSF
;
if keyword_set(psf) and keyword_set(overlap) then begin
  if keyword_set(polar) then begin
    p2c_matrix=shapelets_polar_matrix(n_max,/P2C)
    coeffs=p2c_matrix#coeffs
    coeffs=invert(P_nm)#coeffs
    coeffs=c2p_matrix#coeffs
  endif else begin
    coeffs=invert(P_nm)#coeffs
  endelse
endif


;
; Store results in a structure
;
decomp=shapelets_create_decomp(n_max,polar=polar)
decomp.name=name
decomp.x=x0
decomp.beta=beta
decomp.coeffs=coeffs
decomp.coeffs_error=coeffs_error
decomp.n_pixels=[n_pixels_x,n_pixels_y]
decomp.over=over
decomp.integrate=integrate
decomp.chisq=chisq
decomp.sky_level=skyfit[0]
decomp.sky_slope=skyfit[1:2]


;
; Tell the world
;
return,decomp

end

