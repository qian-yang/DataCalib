pro get_aper_data, aper
aper = 6
data = get_data(aper)
objs = mask_new2(data)
nall = long(n_elements(objs))
print, objs[0]
struct_add_field, objs, 'OIND', lindgen(nall)
;
sort_index = sort(objs.frameid)
objs = objs[sort_index]
struct_add_field, objs, 'framerank', lindgen(nall)
;
sort_index = sort(objs.oind)
objs = objs[sort_index]

mwrfits, objs, 'data_aper6_175_amp.fits'
end
