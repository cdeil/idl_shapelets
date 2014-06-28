;+
; NAME:
;       MRD_STRUCT
; PURPOSE:
;       Return a structure as defined in the names and values data.
; CALLING SEQUENCE:
;       struct = MRD_STRUCT(NAMES, VALUES, NROW, STRUCTYP='name' )
; INPUT PARAMETERS:
;       NAMES   = A string array of names of structure fields.
;       VALUES  = A string array giving the values of the structure
;                 fields.  See examples below.
;       NROW    = The number of elements in the structure array.
;       
; RETURNS:
;       A structure as described in the parameters or 0 if an error
;       is detected.
;
; OPTIONAL KEYWORD PARAMETERS:
;       STRUCTYP = The structure type.  Since IDL does not allow the
;                  redefinition of a named structure it is an error
;                  to call MRD_STRUCT with different parameters but
;                  the same STRUCTYP in the same session.  If this
;                  keyword is not set an anonymous structure is created.
; COMMON BLOCKS:
;       MRD_COMMON
; SIDE EFFECTS:                                                            
;       May create a temporary file if the structure definition is too long 
;       for the EXECUTE function and using old style structures
;
; RESTRICTIONS:
;       By default, the program defines the structure in a long string
;       which is executed with CREATE_STRUCT within a single EXECUTE statement.
;
;       If program is being run in the IDL Virtual machine (EXECUTE statement
;       not allowed), then a separate CREATE_STRUCT statement is called
;       for each tag.   This mode does not have the full capabilities of the
;       normal mode, but should be sufficient for use with MRDFITS().
; PROCEDURE:
;       A structure definition is created using the parameter values.
;       MRD_NSTRUCT is called  and generates the structure in pieces using the
;       execute and create_struct keywords.
;
; EXAMPLES:
;       (1) str = mrd_struct(['fld1', 'fld2'], ['0','dblarr(10,10)'],3)
;           print, str(0).fld2(3,3)
;
;       (2) str = mrd_struct(['a','b','c','d'],['1', '1.', '1.d0', "'1'"],1)
;               ; returns a structure with integer, float, double and string
;               ; fields.
; PROCEDURE CALLS:
;       GETTOK() - needed for virtual machine mode only
; MINIMUM IDL VERSION:
;       V5.3 (uses STRSPLIT)
; MODIFICATION HISTORY:
;       Created by T. McGlynn October, 1994.
;       Modified by T. McGlynn September, 1995.
;          Added capability to create substructures so that structure
;          may contain up to 4096 distinct elements.  [This can be
;          increased by futher iteration of the process used if needed.]
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Removed V4.0 reference to common block  October 1997
;       Allowed unlimited number of structure elements if the version
;       is greater than 5.0.  Put back in code to handle prior versions.
;       The [] will need to be translated back to () for this to
;       work.  T. McGlynn December 15 1998.
;       Add MRD_NSTRUCT since IDL has mysterious problems compiling
;       very large structures.
;       Removed TEMPDIR and OLD_STRUCT keywords  W. Landsman October 2003   
;       Alternate pathway without EXECUTE for V6.0 virtual machine, D. Lindler
;       Removed limit on EXECUTE statement.  W. Landsman  October 203
;-

; Check that the number of names is the same as the number of values.

function mrd_struct, names, values, nrow, $
    structyp=structyp,  tempdir=tempdir, silent=silent, old_struct=old_struct

; Keywords TEMPDIR, SILENT and OLD_STRUCT no longer do anything but are kept
; for backward compatibility.

noexecute = 0b
if !VERSION.RELEASE GE '6.0' then if lmgr(/vm) then noexecute = 1b

if noexecute then begin

    ntags = n_elements(names)
    for i=0,ntags-1 do begin
;
; create variable with the specified data type
;
	case values[i] of 
;
; scalar values
;
	    '0B': v = 0B
	    '0' : v = 0
	    '0L': v = 0L
	    '0LL' : v = 0LL
	    '0.': v = 0.0
	    '0.d0': v = 0.0d0
	    'complex(0.,0.)': v = complex(0.,0.)
	    'dcomplex(0.d0,0.d0)': v = dcomplex(0.d0,0.d0)
;
; strings and arrays
;
	    else: begin
	        value = values[i]
	        type = gettok(value,'(')
		if type eq 'string' then $
			junk = gettok(value,',')      ;remove "replicate(32b"
		dimen_string = gettok(value,')')	
		dimen = long(strsplit(dimen_string,',',/extract))
		case type of
		    'bytarr': v = make_array(dimen=dimen,/byte)
		    'intarr': v = make_array(dimen=dimen,/int)
		    'fltarr': v = make_array(dimen=dimen,/float)
		    'lonarr': v = make_array(dimen=dimen,/long)
		    'lon64arr': v = make_array(dimen=dimen,/long64)
		    'dblarr': v = make_array(dimen=dimen,/double)
		    'complexarr': v = make_array(dimen=dimen,/complex)
		    'dcomplexarr': v = make_array(dimen=dimen,/dcomplex)
                    'ptr_new': v = ptr_new()
		    'string': begin
		    		ndimen = n_elements(dimen)-1
				if ndimen gt 0 then begin
					v = make_array(dimen=dimen[1:*],/string)
					v[*] = string(replicate(32B,dimen[0]))
		    		end else v = string(replicate(32B,dimen[0]))
			      end
		endcase
	    end
	endcase     	
	if i eq 0 then struct = create_struct(name=structyp,names[i],v) $
		  else struct = create_struct(temporary(struct),names[i],v)
    end; for i    

endif else begin


; Check that the number of names is the same as the number of values.

nel = n_elements(names)
if nel ne n_elements(values) then return, 0

; Start formatting the string.
strng = "struct={ " + names[0] + ':' + values[0]

    comma = ","

if nel GE 1 then for i=1,nel-1 do  begin
  
    ; Now for each element put in a name/value pair.
    strng = strng + comma+names[i] + ':' + values[i]
	 
endfor

if strlen(strng) gt 3 then begin
    strng = strng + "}"
    res = execute(strng)
    if  res eq 0 then return, 0
endif

if keyword_set(structyp) then begin
     struct = create_struct(temporary(struct), name=structyp)
endif

  
if nrow le 1 then return, struct

endelse
return, replicate(struct, nrow)

end

