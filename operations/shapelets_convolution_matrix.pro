;$Id: shapelets_convolution_matrix.pro, v2 $
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
; ***********************************************************************
; ***********************************************************************

function shapelets_b3_tensor, a1, a2, a3, n_max_gamma, n_max_alpha, n_max_beta

; NAME:
;      SHAPELETS_B3_TENSOR
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Compute the 3-product tensor B^(3)_lmn(alpha,beta,gamma).
;      This matrix is involved in convolutions, noise estimates and
;      deprojections.
;
; INPUTS:
;      a1,a2,a3         - inverse shapelet scale sizes 
;                         (in the specific order 1/gamma,1/alpha,1/beta)
;
; OPTIONAL INPUTS:
;      n_max            - size of array. Unfortunately has to be square at
;                         the moment; it would be quicker to let the size vary 
;                         in n,m,l directions.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      B^3_lnm(gamma,alpha,beta) tensor from shapelets: II S3.
;
; MODIFICATION HISTORY:
;      Apr 05 - Further optimised by RM to calculate only essential entries.
;      Mar 05 - For loops optimised to use IDL array operations by Will High.
;      Mar 02 - n_max keyword added by Richard Massey.
;      Feb 01 - Written by David Bacon.

COMPILE_OPT idl2, HIDDEN

; Initial declarations
nu=1./sqrt(a1^(-2)+a2^(-2)+a3^(-2))
roottwo=sqrt(2.)
rootpi=sqrt(!pi)
a=roottwo*nu/a1
b=roottwo*nu/a2
c=roottwo*nu/a3

; Construct L_lmn(a,b,c) tensor
L_lmn=fltarr(n_max_gamma+3,n_max_alpha+3,n_max_beta+3)
L_lmn[2,2,2]=1
; Make first row. (Generate l's where m=n=0.)
for l=4,n_max_gamma+2 do L_lmn[l,2,2]=2*(l-3)*(a^2-1)*L_lmn[l-2,2,2]
; Make first l x m array where n=0. (Generate m's where n=0.)
l=lindgen(n_max_gamma+1)
for m=3,n_max_alpha+2 do L_lmn[l+2,m,2]=$
     2*( (m-3)*(b^2-1.)*L_lmn[l+2,m-2,2] + l*a*b*L_lmn[l+1,m-1,2] )
; Make the rest. (Generate n's).
lind=lindgen(n_max_gamma+1)
l=rebin(lindgen(n_max_gamma+1,1),n_max_gamma+1,n_max_alpha+1,/sample)
mind=lindgen(n_max_alpha+1)
m=rebin(lindgen(1,n_max_alpha+1),n_max_gamma+1,n_max_alpha+1,/sample)
for n=3,n_max_beta+2 do L_lmn[lind+2,mind+2,n]=$
    2*( (n-3)*(c^2-1) * L_lmn[lind+2,mind+2,n-2] + $
          l * c * a   * L_lmn[lind+1,mind+2,n-1] + $
          m * b * c   * L_lmn[lind+2,mind+1,n-1]    )
; Get the useful part of the tensor.
L_lmn=L_lmn[2:n_max_gamma+2,2:n_max_alpha+2,2:n_max_beta+2]

; Slow but more explicit method using nested for-loops
;bt=fltarr(n_max+1,n_max+1,n_max+1)
;bt[0,0,0]=1.
;for i=0,n_max do begin
;  for l=0,i do begin
;    for m=0,i do begin
;      for n=0,i do begin
;                
;        if (l ge m) and (l ge n) then begin
;          c1=0
;          c2=0
;          c3=0
;          if ((l-1.) ge 0.) and ((l-2.) ge 0.) then c1=2.*(l-1.)*(a^2-1.)*bt[l-2,m,n]
;          if ((l-1.) ge 0.) and ((m-1.) ge 0.) then c2=2.*m*a*b*bt[l-1,m-1,n]
;          if ((l-1.) ge 0.) and ((n-1.) ge 0.) then c3=2.*n*a*c*bt[l-1,m,n-1]
;        endif
;        
;        if (m ge l) and (m ge n) then begin
;          c1=0
;          c2=0
;          c3=0
;          if ((m-1.) ge 0.) and ((m-2.) ge 0.) then c1=2.*(m-1.)*(b^2-1.)*bt[l,m-2,n]
;          if ((m-1.) ge 0.) and ((l-1.) ge 0.) then c2=2.*l*a*b*bt[l-1,m-1,n]
;          if ((m-1.) ge 0.) and ((n-1.) ge 0.) then c3=2.*n*b*c*bt[l,m-1,n-1]
;        endif
;        
;        if (n ge l) and (n ge m) then begin
;          c1=0
;          c2=0
;          c3=0
;          if ((n-1.) ge 0.) and ((n-2.) ge 0.) then c1=2.*(n-1.)*(c^2-1.)*bt[l,m,n-2]
;          if ((n-1.) ge 0.) and ((l-1.) ge 0.) then c2=2.*l*a*c*bt[l-1,m,n-1]
;          if ((n-1.) ge 0.) and ((m-1.) ge 0.) then c3=2.*m*b*c*bt[l,m-1,n-1]
;        endif
;        
;        bt[l,m,n]=c1+c2+c3
;        
;        if (l eq 0) and (m eq 0) and (n eq 0) then bt[l,m,n]=1.
;                
;      endfor
;    endfor
;  endfor
;endfor

; Evaluate prefactors of the B^(3) tensor
l=rebin(lindgen(n_max_gamma+1,1,1),n_max_gamma+1,n_max_alpha+1,n_max_beta+1,/sample)
m=rebin(lindgen(1,n_max_alpha+1,1),n_max_gamma+1,n_max_alpha+1,n_max_beta+1,/sample)
n=rebin(lindgen(1,1,n_max_beta+1),n_max_gamma+1,n_max_alpha+1,n_max_beta+1,/sample)
prefactor=nu/sqrt(2.^(l+m+n-1)*rootpi*factorial(l)*factorial(m)*factorial(n)*a1*a2*a3)

; Construct B^(3) tensor
bthree=prefactor*temporary(L_lmn)

; Return results to higher scope
return,bthree

end






; ***********************************************************************
; ***********************************************************************

function shapelets_convolution_matrix,psf,alpha,n_max_alpha,gamma,n_max_gamma

;+
; NAME:
;      SHAPELETS_CONVOLUTION_MATRIX
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Create an array of which unconvolved basis functions (with scale
;      alpha) combine to create a basis function convolved with a (decomp
;      structure) PSF (and scale gamma).
;
;          CBF_m(x;a) = Sigma(n=0->inf) { P_nm * BF_n(x;g) }
;       where    P_nm = Sigma(l=0->n_max_beta) { C_nml(a,b,g) * g_l(b) }
;
; INPUTS:
;      PSF         - Cartesian decomp structure representing the local PSF
;      ALPHA       - Scale size of unsmeared object
;      N_MAX_ALPHA - Original BFs go up to this n_max
;
; OPTIONAL INPUTS:
;      GAMMA       - Scale size of convolved basis functions
;      N_MAX_GAMMA - Convolved BFs will be made up of unconvoled BFs up to this
;                    large number (ideally infinite, since noiseless)
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      P_nm        - "PSF matrix" of PaperII eqn 12
;                    Each column represents which BFs(g) a BF(a) transforms 
;                    into under convolution with a PSF(b)
;
; MODIFICATION HISTORY:
;      Jul 05 - Use of 0.5 (1D) in n_max_gamma & gamma changed to 1 (2D) by RM.
;      Apr 05 - Choice of n_max_gamma protected against negative numbers by RM.
;      Mar 05 - Sign error corrected by WH.
;      Mar 05 - Default choice of gamma improved by RM.
;      Mar 05 - For loops optimised to use IDL array operations by Will High.
;      Mar 02 - Written by Richard Massey.
;-


; Initialise properties of the original, unsmeared image f(x;a) or basis 
; functions. Alpha is read in from from the command line (this is the scale
; size for a decomposition that we typically want to find). n_max_alpha is 
; similarly read in from the command line.
shapelets_make_nvec, n_max_alpha, n1_alpha, n2_a, n_coeffs_alpha


; Properties of the PSF g(x;b). Counting variable is l.
beta=psf.beta
n_max_beta=psf.n_max
shapelets_make_nvec, n_max_beta, n1_beta, n2_beta, n_coeffs_beta
;g=psf.coeffs


; Properties of the convolved, observed image h(x;g). Counting variable is n.
if not keyword_set(gamma) then $
  gamma=sqrt(alpha^2+beta^2)
;  gamma=(( (alpha^2*(n_max_alpha+1.)+beta^2*(n_max_beta+1.)) * $
;           (alpha^2*(n_max_beta+1.)+beta^2*(n_max_alpha+1.)) / $
;	   ((n_max_alpha+1.)*(n_max_beta+1.))		    )^0.25)>0
if n_elements(n_max_gamma) eq 0 then $
  n_max_gamma=n_max_alpha>n_max_beta
;  n_max_gamma=round(sqrt( (alpha^2*(n_max_alpha+1.)+beta^2*(n_max_beta+1.)) / $
;		      (alpha^2*(n_max_beta+1.)+beta^2*(n_max_alpha+1.)) * $
;		      ((n_max_alpha+1.)*(n_max_beta+1.))  )-1)>0
shapelets_make_nvec, n_max_gamma, n1_gamma, n2_gamma, n_coeffs_gamma


; Work out B3 tensor using subroutine
bthree=shapelets_b3_tensor(1./gamma,1./alpha,1./beta,n_max_gamma,n_max_alpha,n_max_beta)


; Construct 1-dimensional convolution matrix c (PaperII eqn 7)
; Fast version (by Will High) 
n=rebin(lindgen(n_max_gamma+1,1,1),n_max_gamma+1,n_max_alpha+1,n_max_beta+1,/sample)
m=rebin(lindgen(1,n_max_alpha+1,1),n_max_gamma+1,n_max_alpha+1,n_max_beta+1,/sample)
l=rebin(lindgen(1,1,n_max_beta+1),n_max_gamma+1,n_max_alpha+1,n_max_beta+1,/sample)
evens=where((l+m+n) mod 2 eq 0, n_evens) ; Only evens are nonzero.
c=fltarr(n_max_gamma+1,n_max_alpha+1,n_max_beta+1)
if n_evens gt 0 then c[evens]=sqrt(2.*!pi)*((-1)^((3*n+m+l)[evens]/2.))*bthree[evens] 
;; Slower but more explicit alternative #1 (by David Bacon) 
;c=fltarr(n_max_gamma+1,n_max_alpha+1,n_max_beta+1)
;for n=0,n_max_gamma do begin
;  for m=0,n_max_alpha do begin
;    for l=0,n_max_beta do begin
;      if ((n+m+l) mod 2. eq 0.) then $
;        c[n,m,l]=(2.*!pi)^.5*(-1.)^(1.5*n+.5*m+.5*l)*bthree[n,m,l]
;    endfor
;  endfor
;endfor
;; Slower but more explicit alternative #2 (by Richard Massey)
;; Calculates c(i,j,k)=(2.*!dpi)^.5*(-1.)^i*complex(0,1)^(i+j+k)*bt(i,j,k)
;; and lets c be a complex variable
;c_complex=complexarr(n_max_gamma+1,n_max_alpha+1,n_max_beta+1)
;for n=0,n_max_gamma do begin
;  for m=0,n_max_alpha do begin
;    for l=0,n_max_beta do begin
;      c_complex(n,m,l)=(2.*!dpi)^.5*(-1.)^n*complex(0,1)^(n+m+l)*bthree(n,m,l)
;    endfor
;  endfor
;endfor
;c=c_complex


; Construct 2-dimensional convolution matrix cc (PaperI eqn 52)
n=rebin(lindgen(n_coeffs_gamma,1,1),n_coeffs_gamma,n_coeffs_alpha,n_coeffs_beta,/sample)
m=rebin(lindgen(1,n_coeffs_alpha,1),n_coeffs_gamma,n_coeffs_alpha,n_coeffs_beta,/sample)
l=rebin(lindgen(1,1,n_coeffs_beta),n_coeffs_gamma,n_coeffs_alpha,n_coeffs_beta,/sample)
cc=c[n1_gamma[n],n1_alpha[m],n1_beta[l]]*c[n2_gamma[n],n2_a[m],n2_beta[l]]
;; Slower, older method
;cc=fltarr(n_coeffs_gamma,n_coeffs_alpha,n_coeffs_beta)
;for n=0,n_coeffs_gamma-1 do begin
;  for m=0,n_coeffs_alpha-1 do begin
;    for l=0,n_coeffs_beta-1 do begin
;      cc[n,m,l]=c[n1_gamma[n],n1_alpha[m],n1_beta[l]]*c[n2_gamma[n],n2_a[m],n2_beta[l]]
;    endfor
;  endfor
;endfor


; Calculate the "PSF matrix" (PaperII eqn 12)
P=fltarr(n_coeffs_gamma,n_coeffs_alpha,/NOZERO)
for n=0,n_coeffs_gamma-1 do begin
  for m=0,n_coeffs_alpha-1 do begin
    P[n,m]=total(cc[n,m,*]*psf.coeffs)
  endfor
endfor


; Tell the world our results
return,P

end
