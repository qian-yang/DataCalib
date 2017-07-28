pro get_tref
file = 'sdssrm-bok.fits'
data = mrdfits(file, 1)
struct_add_field, data, 'tref', 0.0
struct_add_field, data, 'dt', 0.0
utdate = data.utDate
index_sort = sort(utdate)
utdate = utdate[index_sort]
index_uniq = uniq(utdate)
date = utdate[index_uniq]
nday = n_elements(date)
ij = 0
for i = 0, nday - 1 do begin
  index = where(data.utDate eq date[i])
  year = strmid(date[i], 0, 4)
  month = strmid(date[i], 4, 2)
  day = strmid(date[i], 6, 2)
  print, date[i], year, month, day
  tref = date_conv(year + '-' + month + '-' + day, 'J') - 2400000.0
  data[index].tref = tref
endfor
data.dt = data.mjd - data.tref
mwrfits, data, 'sdssrm-bok_tref.fits'
end
