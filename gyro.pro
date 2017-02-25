PRO read_gyro, r, file=file

 IF ~keyword_set(file) THEN file='gyro.out'
 openr,unit,file,/get_lun
 IF ~eof(unit) THEN BEGIN 
    readf,unit,nfq
    
    fq=dblarr(nfq)
    jo=dblarr(nfq)
    jx=dblarr(nfq)
    ko=dblarr(nfq)
    kx=dblarr(nfq)

    i=0
    tmp=dblarr(5)
    WHILE ~eof(unit) DO BEGIN 
       readf,unit,tmp
       fq[i]=tmp[0]
       jo[i]=tmp[1]
       jx[i]=tmp[2]
       ko[i]=tmp[3]
       kx[i]=tmp[4]
       i++
    ENDWHILE 
 ENDIF 
 free_lun,unit
 
 r = {fq:fq,jo:jo,ko:ko,jx:jx,kx:kx}
 
END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO gyro,freq,flux $
          ,delta=delta $
          ,energy=energy $
          ,nel=nel $
          ,anor=anor $
          ,ntot=ntot $
          ,m=m $
          ,bmag=bmag $
          ,np=np $
          ,alpha=alpha $
          ,angle=angle $
          ,size=size $
          ,height=height $
          ,temperature=temperature $
          ,e1d=e1d,e2d=e2d,a1d=a1d,a2d=a2d $
          ,phi1=phi1, phi2=phi2 $
          ,electron=electron $
          ,plot=plot $ 
          ,oplot=oplot $
          ,QUIET=QUIET $
          ,ctime=ctime $ 
          ,single_core=single_core $
          ,polariz=polariz $
          ,struct=struct $
          ;,keep=keep,nulls=nulls,Anor_1Mev=Anor_1Mev $
          ,_EXTRA=_EXTRA

 ;; Constants:
 AU = 1.49597870d13                ; Astronomic Unit
 arc2cm = !dtor/3600d0 * AU        ; arcsec to cm in Sun
 m0 = 9.1094d-28                   ; electron mass [g]
 c = 2.998d10                      ; [cm/s] speed of light
 e = 4.803d-10                     ; electron charge
 E0 = m0 * c^2. / 1.6022d-12 / 1d3 ; electron rest energy [keV]

 if not keyword_set(delta) then delta=3.0
 if not keyword_set(energy) then energy=[10,5e3]
 if not keyword_set(m) then m=fix(0) ;; in radians
 if not keyword_set(bmag) then bmag=500.
 if not keyword_set(angle) then angle=45.
 if not keyword_set(np) then np=1e9
 if not keyword_set(size) then size=12.
 if not keyword_set(height) then height=5d8
 if not keyword_set(temperature) then temperature=0.0

 vol = Height * !dpi*(size/2.*arc2cm)^2

IF n_elements(delta) EQ 1 THEN delta=delta[0]

 IF (n_elements(Nel) eq 0) AND (n_elements(Ntot) eq 0) AND (n_elements(Anor) eq 0) THEN $
  nel=1d7

IF n_elements(ntot) NE 0 THEN BEGIN 
   ntot=double(ntot)
   Nel  = Ntot / vol
   anor  = Ntot * (delta-1) / ((Energy[0]/1d3)^(-delta+1) - (Energy[1]/1d3)^(-delta+1)) 
 ENDIF ELSE BEGIN 
    IF n_elements(nel) NE 0 THEN BEGIN  
       Ntot  = Nel * vol
       anor  = Ntot * (delta-1) / ((Energy[0]/1d3)^(-delta+1) - (Energy[1]/1d3)^(-delta+1)) 
    ENDIF ELSE BEGIN 
       IF n_elements(anor) NE 0 THEN BEGIN 
          IF n_elements(energy) GT 2 THEN $
           message,'Anor input only for single power-law. Input Ntot or Nel instead.'
          ;; Anor is the number of electrons at 1 MeV
          ;; since we use energy in keV, we divide by 1E3
          Ntot = anor * ((Energy[0]/1d3)^(-delta+1) - (Energy[1]/1d3)^(-delta+1)) / (delta-1) 
          nel = Ntot / vol
       ENDIF 
    ENDELSE 
 ENDELSE 

print,'Ntot='+string(Ntot)
print,'Anor='+string(anor)
print,'Nel='+string(Nel,form='(e0.6)')

 IF keyword_set(alpha) THEN BEGIN 
    vb = 0.5 * e / (!pi * c * m0) * Bmag ; gyrofrequency
    vp = 1.5*vb / alpha
    np = vp^2 * !dpi *m0 / e^2
 ENDIF 

;; gyro.cpp stuff:

 ned = n_elements(energy)
 nfreq = n_elements(freq)

 phi=0
 gphi=0

if (nel eq 0) and (nfreq gt 0) then begin
print,'User set Nel=0; FREE-FREE ONLY'
  r={fq:freq,jo:fltarr(nfreq),ko:fltarr(nfreq),jx:fltarr(nfreq),kx:fltarr(nfreq)}
endif else begin

 IF n_elements(m) EQ 2 THEN BEGIN 
    ;;print,'gaussian anisotropy'
    p0 = m[0];;center in rad
    p1 = m[1];;width in rad
    pitchangle=[p0,p1]
    mm = fix(1)
    npd = fix(2)
 ENDIF 
 IF n_elements(m) GT 2 THEN BEGIN 
    ;;print,'array anisotropy'
    size_m = size(m)
    IF size_m[2] EQ 2 THEN BEGIN 
       mm=fix(2)
       phi=m[*,0]
       gphi=m[*,1]
       npd=n_elements(phi)
       pitchangle=[[phi],[gphi]]
       ;; p0=!pi/2.*0. ;; temp
       ;; p1=!pi/6. ;; temp
       ;; npd = 100
       ;; phi=interpol([0.0,!pi],npd)
       ;; gphi=exp(-0.5*(phi-p0)^2/(p1)^2)
    ENDIF ELSE BEGIN 
       mm=fix(4)
       npd=fix(m[0])
       pitchangle=m[1:*]
    ENDELSE 
 ENDIF 
 IF n_elements(m) LE 1 THEN BEGIN
    ;;print,'isotropic'
    mm=0
    npd=2 ;; dummy values
 ENDIF 
    
 openw,unit,'gyro.in',/get_lun
 printf,unit, [bmag,angle,np,nel]
 printf,unit, [ned,mm,npd,nfreq]
 IF nfreq NE 0 THEN printf,unit, freq 
 printf,unit, [energy, delta]
 IF mm NE 0 THEN printf,unit,pitchangle
 
 free_lun,unit

 time0=systime(/sec)
 command = '~/solar/gyro/gyro'
 spawn, command, nulls, exit_status=exit_status,stderr=stderr
 ctime=systime(/sec)-time0

 IF ~keyword_set(QUIET) THEN begin
 print,'gyro...'
 print,'Elapsed time:',ctime,' seconds.'
 endif

 read_gyro, r

 ;; remove input and output ascii files after running the program (set
 ;; /keep to keep the files).
 if ~keyword_set(keep) then file_delete,'gyro.in','gyro.out',/ALLOW_NONEXISTENT,/NOEXPAND_PATH,/quiet
endelse
 freq=r.fq

 cgs2sfu = 1e19
 ;; circular area
 area = !pi * (size/2.0*arc2cm)^2.
 volume = area * height
 ;; Source solid angle:
 omega = area / AU^2.
 ;; square area
 ;;th = size * !dtor / 3600.
 ;;omega = th^2.

 gam = energy/E0
 gam1Mev = 1.0D+03/E0
 ; if n_elements(delta) gt 1 then begin
 ;   ien = max(where(energy le 1.0D+03))
 ;   if (ien lt n_elements(energy)) then idelta = ien else idelta = ien - 1
 ;   Anor_1Mev = Anor[ien] * gam1Mev^(-delta[idelta]) * Ntot * gam1Mev 
 ; endif  else  Anor_1Mev = Anor*Ntot*gam1Mev^(-delta) * gam1Mev

 IF KEYWORD_SET(TEMPERATURE) THEN BEGIN
    ;; free-free coefficient by Dulk (1985).
    ;; assuming fully ionized hidrogen isothermal plasma.
    kb = 1.38000e-016           ;; Boltzmann constant
    z2=1.4 ;; mean atomic number (corona)
    t1 = (TEMPERATURE LT 2e5) ? 17.9 + alog(TEMPERATURE^(1.5)) - alog(freq) : $
         24.5 + alog(TEMPERATURE) - alog(freq)
    kff = 9.78e-3 * np^2 / freq^2 / TEMPERATURE^(1.5) * t1 *z2
    jff = kff * kb * TEMPERATURE * FREQ^2 / C^2
    r.jo += jff
    r.jx += jff
    r.ko += kff
    r.kx += kff
;;  print,'FREE-FREE EMISSION/ABSORPTION INCLUDED.'
 ENDIF

;; RADIATIVE TRANSFER
 phi1=dblarr(n_elements(freq))
 phi2=dblarr(n_elements(freq))
 tau_o=r.ko*height
 tau_x=r.kx*height

 toLow = where(tau_o LT 0.001)
 toMed = where((tau_o GE 0.001) AND (tau_o LE 20.))
 toHig = where(tau_o GT 20.)
 txLow = where(tau_x LT 0.001)
 txMed = where((tau_x GE 0.001) AND (tau_x LE 20.))
 txHig = where(tau_x GT 20.)
 IF toLow[0] NE -1 THEN phi1[toLow] = r.jo[toLow] * omega * height
 IF toMed[0] NE -1 THEN phi1[toMed] = omega * r.jo[toMed]/r.ko[toMed]*(1d0-exp(-tau_o[toMed]))
 IF toHig[0] NE -1 THEN phi1[toHig] = r.jo[toHig]/r.ko[toHig] * omega
 IF txLow[0] NE -1 THEN phi2[txLow] = r.jx[txLow] * omega * height
 IF txMed[0] NE -1 THEN phi2[txMed] = omega * r.jx[txMed]/r.kx[txMed]*(1d0-exp(-tau_x[txMed]))
 IF txHig[0] NE -1 THEN phi2[txHig] = r.jx[txHig]/r.kx[txHig] * omega
 
;; phi1 = r.jo/r.ko*(1.0-exp(-r.ko*height))*omega*cgs2sfu
;; phi2 = r.jx/r.kx*(1.0-exp(-r.kx*height))*omega*cgs2sfu
 phi1*=cgs2sfu
 phi2*=cgs2sfu
 flux = phi1+phi2

 ;; phi1 = r.jo/r.ko*(1.0-exp(-r.ko*height))*omega*cgs2sfu
 ;; phi2 = r.jx/r.kx*(1.0-exp(-r.kx*height))*omega*cgs2sfu
 ;; flux = phi1+phi2

 e1d=r.jo
 e2d=r.jx
 a1d=r.ko
 a2d=r.kx

 struct=r
 electron={nel:nel,ntot:ntot,delta:delta,energy:energy}
 
 if keyword_set(plot) then BEGIN

  if keyword_set(polariz) then BEGIN
     yr=[min(flux[where(flux GT 0)]),max(flux[where(flux LT 1e10)])]
     if n_elements(freq) lt 10 then psym=-4
     plot,freq/1e9,flux,/xl,/yl,psym=psym,xs=3 ,_EXTRA=_EXTRA $
         ,ytit='Flux Density [sfu]',yr=yr,pos=[.1,.52,.9,.9]
     plot,freq/1e9,cos(angle*!dtor)/abs(cos(angle*!dtor))*(phi2-phi1)/flux,/xl,xs=3,xtit='Frequency [GHz]',ytit='Polariz. degree',/noerase,pos=[.1,.1,.9,.47]
     plots,10^!x.crange,0.,lines=1
      endif else begin 

    yr=[min(flux[where(flux GT 0)]),max(flux[where(flux LT 1e10)])]
    if n_elements(freq) lt 10 then psym=-4
    plot,freq/1e9,flux,/xl,/yl,psym=psym ,_EXTRA=_EXTRA $
         ,xtit='Frequency [GHz]',ytit='Flux Density [sfu]',yr=yr
    if n_elements(flux_old) gt 0 THEN $
     oplot,freq/1e9,flux_old,col=fsc_color('red'),lines=2 ,_EXTRA=_EXTRA
     endelse

 ENDIF

 IF keyword_set(oplot) THEN oplot,freq,flux,_EXTRA=_EXTRA

END 

FUNCTION gyro,freq $
          ,delta=delta $
          ,energy=energy $
          ,nel=nel $
          ,anor=anor $
          ,ntot=ntot $
          ,m=m $
          ,bmag=bmag $
          ,np=np $
          ,alpha=alpha $
          ,angle=angle $
          ,size=size $
          ,height=height $
          ,temperature=temperature $
          ,e1d=e1d,e2d=e2d,a1d=a1d,a2d=a2d $
          ,phi1=phi1, phi2=phi2 $
          ,electron=electron $
          ,plot=plot $ 
          ,oplot=oplot $
          ,QUIET=QUIET $
          ,ctime=ctime $ 
          ,single_core=single_core $
          ,keep=keep $
          ,_EXTRA=_EXTRA

 gyro,freq,flux $
          ,delta=delta $
          ,energy=energy $
          ,nel=nel $
          ,anor=anor $
          ,ntot=ntot $
          ,m=m $
          ,bmag=bmag $
          ,np=np $
          ,alpha=alpha $
          ,angle=angle $
          ,size=size $
          ,height=height $
          ,temperature=temperature $
          ,e1d=e1d,e2d=e2d,a1d=a1d,a2d=a2d $
          ,phi1=phi1, phi2=phi2 $
          ,electron=electron $
          ,plot=plot $ 
          ,oplot=oplot $
          ,QUIET=QUIET $
          ,ctime=ctime $ 
          ,single_core=single_core $
          ,keep=keep $
          ,_EXTRA=_EXTRA

return, flux

END 
