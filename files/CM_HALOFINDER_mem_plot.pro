pro CM_HALOFINDER_mem_plot, NPART, NFILES, BOUNDARY, BOXSIZE, INTSIZE, FLOATSIZE, MEM_LIMIT, MAXPROC, outfile1

; =========================================================================
; Set colors
red = reform([255, 0, 0], 1, 3)
green = reform([0, 185, 0], 1, 3)
blue = reform([0, 0, 255], 1, 3)
yellow = reform([255, 205, 0], 1, 3)
purple = reform([255, 0, 255], 1, 3)
black = reform([0, 0, 0],1 , 3)
tvlct, red, 1
tvlct, green, 2
tvlct, blue, 3
tvlct, yellow, 4
tvlct, purple, 5
tvlct, black, 6

NPART = double(NPART)
NFILES = double(NFILES)
BOXSIZE = double(BOXSIZE)
INTSIZE = double(INTSIZE)
BOUNDARY = double(BOUNDARY)
MEM_LIMIT = double(MEM_LIMIT)
FLOATSIZE = double(FLOATSIZE)

MEM = dblarr(MAXPROC+1)
NPROC = dindgen(MAXPROC+1);

for i = 0, MAXPROC do MEM[i] = 1.0e30; 

ratio = BOUNDARY/BOXSIZE
for nx = 2l, ceil(MAXPROC/4.0) do begin
  for ny = 2l, ceil(MAXPROC/4.0) do begin
    for nz = 2l, ceil(MAXPROC/4.0) do begin
      if nx*ny*nz le MAXPROC then begin
        extra = 2.0*ratio*((nx+ny+nz) + 2.0*ratio*(nx*ny+nx*nz+ny*nz) + 4.0*ratio*ratio*nx*ny*nz)
        MAXPARTICLES = (1.1*(1.0+extra)*NPART*NPART*NPART)/(nx*ny*nz)
        
        fac_1 = 6.0*FLOATSIZE*MAXPARTICLES                    ; Particle structure
        fac_2 = (6.0*FLOATSIZE*NPART*NPART*NPART)/NFILES      ; Opening and reading input data
        fac_3 = fac_1/NFILES                                  ; Temporary memory for sorted/transferred data
        fac_4 = 4.0*INTSIZE*MAXPARTICLES                      ; FOF algorithm

        MEM_FOF  = fac_1+fac_4
        MEM_READ = fac_1+fac_2+fac_3
        MEM_TRANSFER = fac_1+2.0*fac_3
      
        MEM_TOT = (max([MEM_FOF,MEM_READ,MEM_TRANSFER]))/(1024.0*1024.0*1024.0)
        if MEM_TOT lt MEM[nx*ny*nz] then MEM[nx*ny*nz] = MEM_TOT
      endif
    endfor
  endfor
endfor

index = where(MEM lt MEM_LIMIT)
outputstring = 'Minimum Number of Processors = ' + strtrim(fix(NPROC[index[0]]), 2)

index = where(MEM lt 1.0e30)

print, NPROC[index], MEM[index]

set_plot, 'ps'
device,filename=outfile1,/color

plot, NPROC[index], MEM[index], /nodata, background=6, xrange=[0.0, MAXPROC], yrange=[min(MEM[n_elements(MEM)-1])-0.2, MEM_LIMIT*2], xtitle = 'Number of Processors', ytitle='Memory Per Processor (GByte)', charthick=3, thick=3, xthick=3, ythick=3, charsize=1.5, xstyle=1, ystyle=1
oplot, NPROC[index], MEM[index], color=6, thick=3
oplot, NPROC, dblarr(MAXPROC)+MEM_LIMIT, color = 1, thick = 3, linestyle = 2
xyouts, 3900, 2500, outputstring, charsize = 1.0, charthick = 3, /device

device, /close
set_plot, 'x'

end  


