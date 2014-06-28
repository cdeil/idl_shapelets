function shapelets_interpolate_psf_lsmatrix, degree, x, y,                              $
                                             XRANGE=xrange,                             $
                                             YRANGE=yrange,                             $
                                             MAX_P=max_p,                               $
                                             PLOTIT=plotit,                             $
                                             N_BASIS=n_basis,                           $
                                             BASIS_FUNCTION_VALUE=basis_function_value, $
                                             BASIS_FUNCTION_NAME=basis_function_name,   $
                                             BASIS_FUNCTION_ORDER=basis_function_order


; Currently makes basis functions only for an even+odd Fourier series.
; Could easily be generalised to polynomials/Bessel functions/other bases.
;
; Mar 05 - Written by Richard Massey


; Set boundary conditions
n=n_elements(x)<n_elements(y)
if not keyword_set(xrange) then xrange=[min(x[0:n-1]),max(x[0:n-1])]
if not keyword_set(yrange) then yrange=[min(x[0:n-1]),max(x[0:n-1])]


; Determine number of basis functions that will be needed in order to fit
; to the desired order Fourier series
if keyword_set(max_p) then begin
  n_basis = 2 * ( degree^2 + degree ) +1
endif else begin
  n_basis = 4 * ( degree^2 + degree ) +1
endelse
basis_function_value = dblarr(n_basis, n, /nozero)
basis_function_name  = strarr(n_basis,2)
basis_function_order = intarr(n_basis,2)


; Fill each column of basis functions
k = 0L
for i=0,degree do for j=0,degree do begin ;Fill each column of basis
  if keyword_set(max_p) and (i+j gt degree) then continue
  if k eq 0 then begin
    basis_function_value[k,*] = replicate(1, 1, n)
    basis_function_name[k,*]  = "cos"
    basis_function_order[k,*] = reform([i,j],1,2)
    k=1
  endif else if i eq 0 then begin
    basis_function_value[k,*]     = reform( sin(j*(y[0:n-1]-yrange[0])/(yrange[1]-yrange[0])*!pi), 1, n)
    basis_function_name[k,*]      = ["cos","sin"]
    basis_function_value[k+1,*]   = reform( cos(j*(y[0:n-1]-yrange[0])/(yrange[1]-yrange[0])*!pi), 1, n)
    basis_function_name[k+1,*]    = ["cos","cos"]
    basis_function_order[k:k+1,*] = [replicate(i,2),replicate(j,2)]
    k = k + 2
  endif else if j eq 0 then begin
    basis_function_value[k,*]     = reform( sin(i*(x[0:n-1]-xrange[0])/(xrange[1]-xrange[0])*!pi), 1, n)
    basis_function_name[k,*]      = ["sin","cos"]
    basis_function_value[k+1,*]   = reform( cos(i*(x[0:n-1]-xrange[0])/(xrange[1]-xrange[0])*!pi), 1, n)
    basis_function_name[k+1,*]    = ["cos","cos"]
    basis_function_order[k:k+1,*] = [replicate(i,2),replicate(j,2)]
    k = k + 2
  endif else begin
    basis_function_value[k,*]     = reform( sin(i*(x[0:n-1]-xrange[0])/(xrange[1]-xrange[0])*!pi) * $
                                            sin(j*(y[0:n-1]-yrange[0])/(yrange[1]-yrange[0])*!pi), 1, n)
    basis_function_name[k,*]      = ["sin","sin"]
    basis_function_value[k+1,*]   = reform( cos(i*(x[0:n-1]-xrange[0])/(xrange[1]-xrange[0])*!pi) * $
                                            sin(j*(y[0:n-1]-yrange[0])/(yrange[1]-yrange[0])*!pi), 1, n)
    basis_function_name[k+1,*]    = ["cos","sin"]
    basis_function_value[k+2,*]   = reform( sin(i*(x[0:n-1]-xrange[0])/(xrange[1]-xrange[0])*!pi) * $
                                            cos(j*(y[0:n-1]-yrange[0])/(yrange[1]-yrange[0])*!pi), 1, n)
    basis_function_name[k+2,*]    = ["sin","cos"]
    basis_function_value[k+3,*]   = reform( cos(i*(x[0:n-1]-xrange[0])/(xrange[1]-xrange[0])*!pi) * $
                                            cos(j*(y[0:n-1]-yrange[0])/(yrange[1]-yrange[0])*!pi), 1, n)
    basis_function_name[k+3,*]    = ["cos","cos"]
    basis_function_order[k:k+3,*] = [replicate(i,4),replicate(j,4)]
    k = k + 4
  endelse
endfor


; Make Least-Squares matrix (basis_function_value is M of Anton p460+ or Transpose(M) in Lupton p84)
;if ( keyword_set(noise) ) then $ 
;for i=0L,n_pix_o-1 do basis_function_value[*,i]=basis_function_value[*,i]*noise[i] ; / since V^-1 / since INVERSE variance 
help,basis_function_value
ls_matrix = invert(basis_function_value # transpose(basis_function_value)) # basis_function_value


; c.f. shapelets_decomp.pro
; (Matrix is transpose(M) of Anton p460+ or M in Lupton p84)
;if ( keyword_set(noise) ) then $ 
;for i=0L,n_pix_o-1 do MatrixT[*,i]=MatrixT[*,i]*noise_o[i] ; / since V^-1 / since INVERSE variance 
;MatrixI=invert(MatrixT#Matrix) ; size of matrix to invert is only n_coeffs^2 
;coeffs=MatrixI#MatrixT#data_o 



; Plot the basis functions if requested
if keyword_set(plotit) then begin
  for i=0,n_basis-1 do begin
    !P.region=0
    !P.multi=[0,1,2]
    plot,x,basis_function_value[i,*],psym=1,$
         title="!6 f!d"+strtrim(string(i),2)+"!n(!8x!6,!8y!6)= "+$
               basis_function_name[i,0]+"("+strtrim(string(basis_function_order[i,0]),2)+"!8x!6)"+$
               basis_function_name[i,1]+"("+strtrim(string(basis_function_order[i,1]),2)+"!8y!6)",$
         xtitle="!6x [pixels]",ytitle="!6Basis function #"+strtrim(string(i),2)
    plot,y,basis_function_value[i,*],psym=1,$
         xtitle="!6y [pixels]",ytitle="!6Basis function #"+strtrim(string(i),2)
    read,junk
  endfor
endif


; Tell the outside world
return,ls_matrix

end
