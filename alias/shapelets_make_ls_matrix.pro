;$Id: shapelets_make_ls_matrix.pro, v2$
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

pro shapelets_add_zeros,matrix,matrix2

; NAME:
;       SHAPELETS_ADD_ZEROS
;
; CATEGORY:
;       A component of shapelets_make_ls_matrix.pro.
;
; PURPOSE:
;       Add back in two rows of zeros to the covariance matrix if
;       a_01 and a_10 had previsiously  been contrained to be zero 
;       (probably using shapelets_remove_zeros: see below).
;
; INPUTS:
;       matrix  - The (as yet incomplete) SQUARE covariance matrix
;                 (MatrixI in shapelets_decomp.pro).
;       matrix2 - The matrix of basis functions and pixel values.
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       matrix  - The same matrix, but with two extra rows and two
;                 extra columns of zeros.
;
; MODIFICATION HISTORY:
;       Sep 03 - Combined with shapelets_make_ls_matrix.pro.
;       Apr 02 - Second matrix added, to allow reconstruction of image.
;       Nov 01 - Written by Richard Massey.

COMPILE_OPT idl2, HIDDEN

n=size(matrix)
n=n[1]+2

; add (easy) horizontal row
zero=fltarr(2,n-2)
matrix=[matrix[0,*],zero,matrix[1:*,*]]

;add (difficult) vertical column - must be a better notation
;like zero=fltarr(n,2) etc...? SEE BELOW!

matrix_new=fltarr(n,n)
for i=0,n-1 do begin
 matrix_new[i,0]=matrix[i,0]
 matrix_new[i,1]=0.
 matrix_new[i,2]=0.
 for j=3,n-1 do matrix_new[i,j]=matrix[i,j-2]
endfor
matrix=matrix_new

; add (difficult) vertical column to second matrix
; nice new method is to transpose the matrix, do the same as before,then
;  transpose back!
n=total(size(matrix2)*[0,1])
zero=fltarr(2,n)
matrix2=transpose([transpose(matrix2[*,0]),zero,transpose(matrix2[*,1:*])])

end




; ***********************************************************************
; ***********************************************************************

pro shapelets_extend_zeros,matrix,sky

; NAME:
;       SHAPELETS_EXTEND_ZEROS
;
; CATEGORY:
;       A component of shapelets_make_ls_matrix.pro.
;
; PURPOSE:
;       Extend the PSF matrix to allow sky fitting in shapelets_decomp.pro
;
; INPUTS:
;       Matrix - Least squares fitting matrix.
;
; KEYWORD PARAMETERS:
;       /SKY - 1: simulataneously fit a flat constant sky value.
;              2: simulataneously fit a plane with variable slope.
;
; OUTPUTS:
;       Matrix - Modified least squares fitting matrix. 
;                If SKY=2, the output will be:
;     
;                  [              |   ]
;                  [     P_nm     | 0 ]
;                  [              |   ]
;                  [--------------|---]
;                  [              |100]
;                  [       0      |010]
;                  [              |001]
;
; MODIFICATION HISTORY:
;       Nov 01 - Written by Richard Massey.
      
COMPILE_OPT idl2, HIDDEN

if sky eq 1 then n=1 else n=3

msize=size(matrix)
h=msize[1]
v=msize[2]

newmatrix=fltarr(h+n,v+n)

for i=0,h-1 do begin
  for j=0,v-1 do begin
    newmatrix[i,j]=matrix[i,j]
  endfor
endfor

for i=0,n-1 do newmatrix[h+i,v+i]=1.

matrix=newmatrix

return

end



; ***********************************************************************
; ***********************************************************************

pro shapelets_remove_zeros,matrix,trans=trans

; NAME:
;       SHAPELETS_REMOVE_ZEROS
;
; CATEGORY:
;       A component of shapelets_make_ls_matrix.pro.
;
; PURPOSE:
;       Remove two columns from the PSF matrix which correspond to n=1 states.
;       These will then not be fit in shapelets_decomp.pro.
;
; INPUTS:
;       Matrix - Least squares fitting matrix.
;
; KEYWORD PARAMETERS:
;       /TRANS - Takes transpose first, to remove rows instead of columns.
;                Internal flag, used if the input matrix is really the 
;                transpose of that expected.
; OUTPUTS:
;       Matrix - Modified least squares fitting matrix:
;
;          [                   ]
;          [  P_nm or MatrixT  ]
;          [                   ]
;           | XX        |
;           | XX        |
;           V XX        V
;          [ ]  [              ]
;          [ ]..[     P_nm     ]
;          [ ]  [              ]
;            \        |
;             \       |
;              _|     V
;            [                 ]
;            [      P_nm       ]
;            [                 ]
;
; MODIFICATION HISTORY:
;       Nov 01 - Written by Richard Massey.

COMPILE_OPT idl2, HIDDEN

if keyword_set(trans) then matrix=transpose(matrix)

  matrix[2,*]=matrix[0,*]
  matrix=matrix[2:*,*]

if keyword_set(trans) then matrix=transpose(matrix)

return

end





; ***********************************************************************
; ***********************************************************************

pro shapelets_make_ls_matrix, MATRIX, N_MAX, BETA, N_PIXELS, X0, $
                              NON1=non1,                         $
                              SKY=sky,                           $
			      NOSKYGRAD=noskygrad,               $
			      INTEGRATE=integrate,               $
			      PSF=psf,                           $
			      TRANS=trans,                       $
			      MATRIX2=matrix2,                   $
			      POLAR=polar,                       $
			      DIAMOND=diamond

;+
; NAME:
;       SHAPELETS_MAKE_LS_MATRIX
;
; CATEGORY:
;       Linear algebra.
;
; PURPOSE:
;       Create an array with which to compute the least-squares linear algebra
;       fit of shapelets basis functions to data. Can also use this array to
;       perform the overlap integrals calculation.
;
; CALLING PROCEDURE:
;       shapelets_make_ls_matrix, MatrixI, n_max, beta, n_pixels, x0
;
; INPUTS:
;       n_max      - maximum order of the hermite ceofficients to compute
;       beta       - scalar basis function scale size
;       n_pixels   - [x,y] integer number of pixels in image
;       x0         - [x,y] vector centre coordinates of the basis functions
;
; OPTIONAL INPUTS:
;       psf   - Cartesian shapelet decomp structure of the local PSF
;               for deconvolution
;
; KEYWORD PARAMETERS:
;       /INT  - use basis functions integrated within square pixels
;       /SKY  - fit sky background with plane
;       /NOSKY- fit sky with just a constant
;
; OUTPUTS:
;       Matrix- least squares fitting matrix
;
; NOTES:
;       Perhaps it would have been easier to just pass the grid coordinate
;       positions of pixel centres (usually called x1 and x2), rather than
;       the number of pixels and their centre. The grids now have to be
;       reconstructed on the first line. Oh well.
;
; MODIFICATION HISTORY:
;       Aug 05 - Obsoleted by RM.
;       Apr 05 - Hacked to allow users to call internal subroutines by RM.
;       Sep 03 - Other routines involving the addition or subtraction of
;                rows contaiing zeros (to constrain some fit parameters)
;                combined into one file by RM.
;       Nov 01 - Written by Richard Massey.
;-

COMPILE_OPT idl2, OBSOLETE

; allow user to call subroutines from outside
if keyword_set(add_zeros) then begin
  shapelets_add_zeros,matrix,matrix2
  return
endif
if keyword_set(extend_zeros) then begin
  shapelets_extend_zeros,matrix,sky
  return
endif
if keyword_set(remove_zeros) then begin
  shapelets_remove_zeros,matrix,trans=trans
  return
endif

; construct x-y coordinate array
shapelets_make_xarr, [n_pixels[0], n_pixels[1]], x1, x2, x0=x0

; calculate basis function arrays
if keyword_set(polar) then begin
  Basis=shapelets_chi([n_max,n_max],x1/beta,x2/beta,integrate=integrate,/array)/beta
  shapelets_make_nvec, n_max, n, m, n_coeffs, polar=polar
  Basis[*,*,where(m lt 0)]=Basis[*,*,where(m lt 0)]*complex(0,1)
  Basis=float(Basis)
  if keyword_set(diamond) then Basis=Basis[*,*,where(n+abs(m) le n_max,n_coeffs)]
endif else begin
  Basis=shapelets_phi([n_max,0],x1/beta,x2/beta,integrate=integrate,/array)/beta
  n_coeffs=(size(basis,/dimensions))[2]
endelse
  help,n_max,n_coeffs

; work out ns and if necessary add extra fitting parameter for sky background
n_coeffs_tot=n_coeffs
if keyword_set(sky) then begin
  n_coeffs_tot=n_coeffs_tot+1
  if ( sky eq 2 ) then n_coeffs_tot=n_coeffs_tot+2
endif

; create LS matrix
Matrix=fltarr(n_coeffs_tot, n_pixels[0]*n_pixels[1], /nozero)

; put basis functions into Matrix
; SLOW!
j=0L                         ; loop over all pixels in image
for f2i=0,n_pixels[1]-1 do begin
  for f1i=0,n_pixels[0]-1 do begin
    for i=0L,n_coeffs-1 do begin   ; loop over all basis functions
      Matrix[i,j]=Basis[f1i,f2i,i]
    endfor
    j=j+1
  endfor
endfor

; add extra basis functions to map sky background
if keyword_set(sky) then begin
  Matrix[n_coeffs,*]=1.                   ; constant
  if ( sky eq 2 ) then begin
    j=0L
    ; SLOW!
    for f2i=0,n_pixels[1]-1 do begin ; loop over all pixels in image
      for f1i=0,n_pixels[0]-1 do begin
        Matrix[n_coeffs+1,j]=x1[f1i,f2i]  ; x (for a plane)
        Matrix[n_coeffs+2,j]=x2[f1i,f2i]  ; y
        j=j+1
      endfor
    endfor
  endif
endif

; remove some basis functions if we are requiring a_01 and a_10 to be zero
if keyword_set(non1) and n_max ge 2 then shapelets_remove_zeros,Matrix

end

