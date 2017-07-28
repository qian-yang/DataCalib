pro get_nightindex_g
file = 'sdssrm-bok_tref.fits'
data = mrdfits(file, 1)
struct_add_field, data, 'dayindex', -1
struct_add_field, data, 'nightindex', -1
utdate = data.utDate
index_sort = sort(utdate)
utdate = utdate[index_sort]
index_uniq = uniq(utdate)
date = utdate[index_uniq]
nday = n_elements(date)
data.imtype = STRCOMPRESS(data.imtype, /REMOVE_ALL)
data.filter = STRCOMPRESS(data.filter, /REMOVE_ALL)
ij = 0
for i = 0, nday - 1 do begin
  index = where(data.utDate eq date[i])
  data[index].dayindex = i
  ; ind = where((data.dayindex eq i) and (data.good eq 'T') and (data.imtype eq 'object') and (data.filter eq 'g'), none)
  ind = where((data.dayindex eq i) and (data.imtype eq 'object') and (data.filter eq 'g'), none)
  help, ind
  if (none gt 0) then begin
    data[ind].nightindex = ij
    ij += 1
  endif
endfor
mwrfits, data, 'sdssrm-bok_tref_nightg_2.fits'
end
