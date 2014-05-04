pro CM_HALOFINDER_mem, NPART, NFILES, BOUNDARY, BOXSIZE, INTSIZE, FLOATSIZE, NPROC

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
FLOATSIZE = double(FLOATSIZE)

MEM_BEST = 1.0e30
ratio = BOUNDARY/BOXSIZE
for nx = 2l, ceil(NPROC/4.0) do begin
  for ny = 2l, ceil(NPROC/4.0) do begin
    for nz = 2l, ceil(NPROC/4.0) do begin
      if nx*ny*nz eq NPROC then begin
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
        if MEM_TOT lt MEM_BEST then begin
          NX_BEST = nx
          NY_BEST = ny
          NZ_BEST = nz
          MEM_BEST = MEM_TOT
        endif
        print, nx, ny, nz, MEM_TOT, 'GB'
      endif
    endfor
  endfor
endfor

print, 'Lowest memory usage:'
print, '  Nx = ', NX_BEST
print, '  Ny = ', NY_BEST
print, '  Nz = ', NZ_BEST
print, '  Memory = ', MEM_BEST, 'GB'

end  


