pro get_gframes
file = 'sdssrm-bok_tref_nightg.fits'
data = mrdfits(file, 1)
index = where((data.nightindex ge 0) and (data.good eq 'T') and (data.frameindex le 8069), nf)
print, nf
data = data[index]
index_sort = sort(data.frameindex)
data = data[index_sort]
struct_add_field, data, 'i1', 0l
struct_add_field, data, 'i2', 0l
struct_add_field, data, 'number', 0l
file = 'data_aper6_175_amp.fits'
objs = mrdfits(file, 1)
index_sort = sort(objs.framerank)
objs = objs[index_sort]
for i = 0, nf - 1 do begin
  ind = where(objs.frameid eq data[i].frameindex, none)
  data[i].number = none
  data[i].i1 = min(ind)
  data[i].i2 = max(ind)
  print, i, data[i].frameindex, data[i].i1, data[i].i2, none, data[i].i2-data[i].i1+1
endfor
filew = 'frames_g.fits'
spawn, 'rm -rf '+filew
mwrfits, data, filew
end
