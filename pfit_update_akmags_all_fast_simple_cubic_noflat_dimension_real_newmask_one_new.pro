function get_1d_cubic, x, y
   x = float(x)
   y = float(y)
   nall = n_elements(x)
   result = [transpose(x^3), transpose(y^3), transpose(x^2*y), transpose(x*y^2), $
   transpose(x^2), transpose(y^2), transpose(x*y), $
   transpose(x), transpose(y), $
   transpose(fltarr(nall)+1.0)]
   return, result
end

function get_atcinvb, npars, nthis, aind, kind, kdotind, thisb_wbvec, thisb_wbvec_air, thisb_wbvec_air_dt
    atcinvbone = dblarr(npars)
    for k=0L, nthis-1 do begin
      ; A term part
      atcinvbone[aind[k]] += thisb_wbvec[k]
      ; K term part
      atcinvbone[kind[k]] -= thisb_wbvec_air[k]
      ; Kdot term part
      atcinvbone[kdotind] -= thisb_wbvec_air_dt[k]
      ; flat term part
      ; fd = flatind[k]
      ; atcinvbone[fd: fd+9] += thisb_wbvec[k] * xy[*, k]
    endfor
    return, atcinvbone
end

function get_atcinva, npars, nthis, aind, kind, kdotind, wtprod, thisair, dt
   atcinvaone = dblarr(npars, npars)
   for i=0L, nthis-1 do begin
      for j=0L, nthis-1 do begin
         wtpd = wtprod[i, j]
         ai = aind[i]
         aj = aind[j]
         ki = kind[i]
         kj = kind[j]
         dij = kdotind
         ; fi = flatind[i]
         ; fj = flatind[j]
         ; xyi = xy[*, i]
         ; xyj = xy[*, j]
         airi = thisair[i]
         airj = thisair[j]
         dti = dt[i]
         dtj = dt[j]
         ;
         ; AA
         atcinvaone[ai, aj] += wtpd
         ; KK
         atcinvaone[ki,kj] += wtpd * airi * airj
         ; KA
         atcinvaone[ki,aj] -= wtpd * airi
         ; AK
         atcinvaone[ai,kj] -= wtpd * airj
         ; Ad
         atcinvaone[ai,dij] -= wtpd * airj * dtj
         ; dA
         atcinvaone[dij,aj] -= wtpd * airi * dti
         ; Kd
         atcinvaone[ki,dij] += wtpd * airi * airj * dtj
         ; dK
         atcinvaone[dij,kj] += wtpd * airi * dti * airj
         ; dd
         atcinvaone[dij,dij] += wtpd * airi * dti * airj * dtj
         ; ;---------- f -----------
         ; ; Af
         ; atcinvaone[ai, fj:fj+9] += wtpd * xyj
         ; ; fA
         ; atcinvaone[fi:fi+9, aj] += wtpd * xyi
         ; ; kf
         ; atcinvaone[ki, fj:fj+9] -= wtpd * airi * xyj
         ; ; fk
         ; atcinvaone[fi:fi+9, kj] -= wtpd * airj * xyi
         ; ; df
         ; atcinvaone[di, fj:fj+9] -= wtpd * airi * dti * xyj
         ; ; fd
         ; atcinvaone[fi:fi+9, dj] -= wtpd * airj * dtj * xyi
         ; ; ff
         ; atcinvaone[fi:fi+9, fj:fj+9] += wtpd * (xyi ## xyj)
            endfor
         endfor
return, atcinvaone
end

;+
; Reference:
;   learn from the pfit_update_akmags routine by Finkbeiner, Padmanabhan & Schlegel
; PURPOSE:
;   return best-fit aterms and kterms, given everything else
; CALLING SEQUENCE:
;   pfit_update_akmags, objs, magarr, uparam, filter=, catmagivar=
; INPUTS:
;   objs      - Array of OBJS structures, one per observation
;   uparam    - Structure with uber-calibration parameters
;-
;------------------------------------------------------------------------------
pro pfit_update_akmags_all_fast_simple_cubic_noflat_dimension_real_newmask_one_new, objs, uparam, uparam_out = uparam_out, objs_out = objs_out

   common com_pcalib_flist, flist
   eps = (machar()).eps
   ; Get dimensions
   nobj = n_elements(objs)
   naterms = long(max(objs.aindx))+1
   nkterms = long(max(objs.kindx))+1
   ; nkdoterms = nkterms

   ; SVD parameters
   ; svdmin = 1.e-5 ;;; This is used for constraining the nullspace
   if NOT keyword_set(svthresh1) then begin
      svdthresh=1.0d4
   endif else begin
      svdthresh = svdthresh1   ; also somewhat arbitrary...
   endelse

   npars = naterms + nkterms + 1;nkdoterms; + nflaterms
   alist = lindgen(naterms)
   klist = lindgen(nkterms)
   ; kdotlist = lindgen(nkdoterms)
   ;----------
   ; Construct the matrix terms...

   atcinvb = dblarr(npars)
   atcinva = dblarr(npars, npars)

   file = 'sdssRefStars_DR13_g.fits'
   dd = mrdfits(file, 1)
   nstar = long(n_elements(dd))
   objidarr = dd.objid

   t1 = systime(1)
   ; nbefore = naterms + nkterms + nkdoterms
   ; flatarr = uparam.flatarr


   for istar=0L, nstar-1 do begin
    print, objidarr[istar]
      ; Get all observations of this star
      i1 = dd[istar].i1
      i2 = dd[istar].i2
      idjj = lindgen(dd[istar].number) + i1
      objj = objs[i1:i2]
      ; indx = where((objj.bmask eq 0),nthis)
      ; indx = where((objj.bmask eq 0) and (objj.mask eq 1),nthis)
      indx = where((objj.cal eq 1),nthis)

      if (nthis GT 1) then begin ; Don't bother if only 1 good observation
         theseobjs = objj[indx]
         ;ind
         aind = alist[theseobjs.aindx]
         kind = klist[theseobjs.kindx]+naterms;*(ngoodk gt 0)
         kdotind = naterms + nkterms ;kdotlist[theseobjs.kindx] +
         ; flatind = nbefore + (theseobjs.ccdnum - 1) * 10
         ;
         thisair = double(theseobjs.airmass)
         dt = theseobjs.dt;theseobjs.mjd - floor(theseobjs.mjd)
         ;dt = theseobjs.dt
         ; flatindx = theseobjs.flatindx
         ; x = theseobjs.xall
         ; y = theseobjs.yall
         ; -------- Assert that dt is good
         if (max(abs(dt)) GT 12.0) then message, 'Weird dt values:',dt
         instmag = theseobjs.instmag + theseobjs.zeropoint
         wvec = theseobjs.instmagivar
         wvec = wvec / total(wvec, /double) ; Make these normalized weights w/happy properties
         meanmag = total(wvec * instmag)
        ;  objs[idjj].refmag = meanmag
         if (max(wvec)+10*eps) GE 1.0 then $
           message, 'wvec too large!' ; this should never happen!

         ; Construct "C^-1 * b" (from A*x=b)
           thisb = (meanmag - instmag)* theseobjs.instmagivar

         ; Construct "A^T * C^-1 * b"
         wbvec = total(thisb)*wvec
         thisb_wbvec = thisb - wbvec
         thisb_wbvec_air = (thisb - wbvec)*thisair
         thisb_wbvec_air_dt = (thisb - wbvec)*thisair*dt
         ; xy = get_1d_cubic(x, y)
         atcinvbone = get_atcinvb(npars, nthis, aind, kind, kdotind, thisb_wbvec, thisb_wbvec_air, thisb_wbvec_air_dt)
         atcinvb += atcinvbone

         ; Construct "A^T * C^-1 * A"
         ; subarray of A transpose
         wts = identity(nthis)-(fltarr(nthis)+1)#wvec
         ; subarray of A
         wtst = transpose(wts)

         for i = 0L, nthis-1L do wtst[*,i] *= theseobjs[i].instmagivar
         wtprod = wts##wtst

         atcinvaone = get_atcinva(npars, nthis, aind, kind, kdotind, wtprod, thisair, dt)
         atcinva += atcinvaone
      endif
   endfor

   t2 = systime(1)
   splog, 'Time to build matrices:   ', t2-t1, ' sec'
   ;----------
   ; Solve the matrix

   splog, 'Inverting matrix'
      svdc, atcinva, ww, uu, vv, /double, itmax=100
      print, 'WW:'
      print, ww
      print
      svdsave = {atcinva:atcinva, atcinvb:atcinvb, ww:ww, uu:uu, vv:vv, svdthresh:svdthresh}
      wdegen = where(ww LT svdthresh, ndegen)
      if (ndegen GT 0) then begin
         ww[wdegen] = 0.0d0
         splog, 'SVD: Setting ', ndegen, ' elements of ww to zero'
      endif
      par = svsol(uu, ww, vv, atcinvb, /double)

   ; Don't update parameters that had no stars to constrain the fit
   print, '------------ a ------------'
   print, par[0:naterms-1];uparam.aterm
   print, '------------ k ------------'
   print, par[naterms:naterms+nkterms-1];uparam.kterm
   print, '------------ kdot ------------'
   print, par[naterms+nkterms];uparam.kdoterm
   ; print, '------------ flat ------------'
   ; print, par[nbefore: npars - 1]
   print, '------------------------------'
    uparam.aterm = par[0:naterms-1]
    uparam.kterm = par[naterms:naterms+nkterms-1]
    uparam.kdoterm = par[naterms+nkterms]
    ; uparam.flaterm = par[nbefore: npars - 1]


   splog, 'Time for everything else: ', systime(1)-t2,  ' sec'
   uparam_out = uparam
   ; mwrfits, uparam_out, 'uparam_all_noflat.fits'
   ; objs_out = objs
   return
end
;------------------------------------------------------------------------------
