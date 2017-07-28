pro get_ak_indx_all_newmask_nightly, objs, framedata, objs_out = objs_out
n_night = max(framedata.nightindex) + 1
objs.aindx = -1
objs.kindx = -1
sort_index = sort(objs.framerank)
objs = objs[sort_index]
ki = 0
for j = 0, n_night - 1 do begin
  print, "night ", j
  index = where((framedata.nightindex eq j), nfs)
  ;only check when whole night is not masked
  if (framedata[index[0]].wholenight eq 0) then begin
    i1 = framedata[index[0]].i1
    i2 = framedata[index[nfs - 1]].i2
    id = lindgen(i2 - i1 + 1) + i1
    objs[id].kindx = ki
    objs[id].aindx = ki * 4 + objs[id].CCDNUM - 1
    ki += 1
  endif else begin
    print, "the night was masked"
  endelse
endfor
sort_index = sort(objs.oind)
objs = objs[sort_index]
objs_out = objs
end
