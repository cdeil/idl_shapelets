function shapelets_split_sexcat, sexcat_in, selected
COMPILE_OPT IDL2, HIDDEN, OBSOLETE
sexcat_out=sexcat_in
shapelets_split,sexcat_out,selected
if n_tags(sexcat_out) ne n_tags(sexcat_in) then print,$
    "WARNING: not all tags transferred to separated catalogue!"
return,sexcat_out
end
