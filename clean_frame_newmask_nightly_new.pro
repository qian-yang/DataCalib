pro clean_frame_newmask_nightly_new, objs, framedata, framedata_out = framedata_out, objs_out = objs_out, nc = nc, masknight = masknight
offset = objs.mcal - objs.sdss_mag
n_night = max(framedata.nightindex) + 1
;rank by frameid
sort_index = sort(objs.framerank)
objs = objs[sort_index]
offset = offset[sort_index]
nc = 0
for j = 0, n_night - 1 do begin
    print, "night ", j
    index = where((framedata.nightindex eq j), nfs)
    ;only check when whole night is not masked
    if (framedata[index[0]].wholenight eq 0) then begin
    for i = 0, nfs - 1 do begin
    ;check when the frame is not selfmask non-photometric masked
    if (framedata[index[i]].selfmask eq 0) then begin
        print, "frame", framedata[index[i]].frameindex
        i1 = framedata[index[i]].i1
        i2 = framedata[index[i]].i2
        id = lindgen(i2 - i1 + 1) + i1
        obj = objs[id]
        oft = offset[id]
        tp = where((obj.cal eq 1) or ((obj.flagq mod 13) eq 0) or ((obj.flagq mod 331) eq 0), ninde)
        if (ninde gt 0) then begin
            framedata[index[i]].cal_mask = ninde
            if (ninde gt 20) then begin
                ind = where(((obj.cal eq 1) and ((oft gt 0.1))) or ((obj.flagq mod 13) eq 0) or (((obj.flagq mod 331) eq 0) and ((oft gt 0.1))), n1)
                framedata[index[i]].non_number = n1
                framedata[index[i]].non_ratio = float(n1)/float(ninde)
                if ((n1 gt 0.5*ninde)) then begin
                    print, "!! mask frame -- ", 'number & ratio', n1, ninde, float(n1)/float(ninde), ' ; frameid ', framedata[index[i]].frameindex
                    nc += 1
                    framedata[index[i]].photometric = 0
                    framedata[index[i]].selfmask = 1
                    ; mask all the frame
                    tparr = objs[id].flagq mod 317
                    tpind = where(tparr gt 0, tpn)
                    if (tpn gt 0) then begin
                    objs[id[tpind]].flagq *= 317
                    objs[id[tpind]].cal = 0
                    objs[id[tpind]].stat = 0
                    endif
                    ; objs[id].remain = 1
                endif else begin
                    print, "Normal for calibration", ninde
                endelse
            endif else begin
                print, "!! mask frame < 20 -- ", ninde, ' ; frameid ', framedata[index[i]].frameindex
                nc += 1
                framedata[index[i]].photometric = 0
                framedata[index[i]].selfmask = 1
                tparr = objs[id].flagq mod 23
                tpind = where(tparr gt 0, tpn)
                if (tpn gt 0) then begin
                objs[id[tpind]].flagq *= 23
                objs[id[tpind]].cal = 0
                objs[id[tpind]].stat = 0
                endif
                ; objs[id].remain = 1
            endelse
        endif else begin
            print, "!No useful objects in this frame, shouldn't happen when selfmask = 0"
        endelse
    endif
    endfor
    endif
endfor
masknight = 0
for j = 0, n_night - 1 do begin
    print, "night ", j
    index = where((framedata.nightindex eq j), nfs)
    ;only check when whole night was not masked
    if (framedata[index[0]].wholenight eq 0) then begin
    tp = where(framedata[index].photometric eq 0, n_no)
    print, "total frames", nfs, " non-photometric frames", n_no
    if (n_no gt 0) then begin
        if (n_no gt 0.5*nfs) then begin
            print, "!!! Mask the whole non-photometric night"
            framedata[index].wholenight = 1
            framedata[index].photometric = 0
            i1 = framedata[index[0]].i1
            i2 = framedata[index[nfs - 1]].i2
            id = lindgen(i2 - i1 + 1) + i1
            tparr = objs[id].flagq mod 317
            tpind = where(tparr gt 0, tpn)
            if (tpn gt 0) then begin
            objs[id[tpind]].flagq *= 317
            objs[id[tpind]].cal = 0
            objs[id[tpind]].stat = 0
            endif
            ; objs[id].remain = 1
            masknight += 1
        endif else begin
            for i = 0, nfs - 1 do begin
                if (framedata[index[i]].selfmask eq 1) then begin
                    print, "check selfmask frame", framedata[index[i]].frameindex
                    if (i gt 0) then begin
                        if (framedata[index[i-1]].photometric eq 1) then begin
                            print, "Mask former image"
                            i1 = framedata[index[i-1]].i1
                            i2 = framedata[index[i-1]].i2
                            id = lindgen(i2 - i1 + 1) + i1
                            framedata[index[i-1]].photometric = 0
                            idx = where(objs[id].cal eq 1, n_idx)
                            if (n_idx gt 0) then begin
                                objs[id[idx]].flagq *= 331
                                objs[id[idx]].cal = 0
                                objs[id[idx]].stat = 0
                                ; objs[id[idx]].remain = 1
                            endif
                        endif
                    endif
                    if (i lt nfs - 1) then begin
                        if (framedata[index[i+1]].photometric eq 1) then begin
                            print, "Mask latter image"
                            i1 = framedata[index[i+1]].i1
                            i2 = framedata[index[i+1]].i2
                            id = lindgen(i2 - i1 + 1) + i1
                            framedata[index[i+1]].photometric = 0
                            idx = where(objs[id].cal eq 1, n_idx)
                            if (n_idx gt 0) then begin
                                objs[id[idx]].flagq *= 331
                                objs[id[idx]].cal = 0
                                objs[id[idx]].stat = 0
                                ; objs[id[idx]].remain = 1
                            endif
                        endif
                    endif
                endif
            endfor
        endelse
    endif else begin
        print, "Very nice photometric night!"
    endelse
    endif else begin
        print, "Whole night already masked"
        masknight += 1
    endelse
endfor
; mwrfits, framedata, 'framedata_correction_2_nofit.fits'
sort_index = sort(objs.oind)
objs = objs[sort_index]
objs_out = objs
framedata_out = framedata
end
