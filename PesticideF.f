	program main
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50, nk = 101
	integer, parameter :: nc = 5, nb = 3, nl = 2
      integer, parameter :: iiniter = 1000, initer = 20000
      integer, parameter :: niter = 100000, niter1 = 20000
	      
	real*8 cancer, base(nc-1)
      real*8 durat(nj), ddurat
      real*8 zmean(nc), zstd(nc)
	
	real*8 yi(ni), zi(ni,nc), xij(ni,nj,nb)    
      real*8 wi(ni), taui(ni), vnu 
	real*8 betaa(nc), alphaa(nb,nl)
      real*8 azeta(nl), agamma
      real*8 gj(nj), weight(nl), theta(nl)
      
      real*8 sigma_betaa, sigma_theta, a0, b0
      
      integer ngj(nl), jk(nj)
      real*8 xa, temp, sumw
      real*8 xstar(nj,nk,nb)    

      real*8 cngj(nj,nl)
      real*8 mu(nj,nk)
      real*8 mul(nj,nk,nl)
      real*8 prob(nj,nl)
      
      real*8 mean_mu(nj,nk)
      real*8 mean_mul(nj,nk,nl)
      real*8 mean_prob(nj,nl)

      real*8 std_mu(nj,nk)
      real*8 std_prob(nj,nl)

      real*8 lower_mu(nj,nk)
      real*8 lower_prob(nj,nl)

      real*8 upper_mu(nj,nk)
      real*8 upper_prob(nj,nl)

      real*8 seq_mu(niter1,nj,nk)
      real*8 seq_prob(niter1,nj,nl)      
      
      real*8 ahpd(niter1), alow(2), aupp(2)
            
      real*8 delta(nj,nb)
	real*8 mean_betaa(nc)
      real*8 mean_delta(nj,nb)      
	real*8 bardic, dicbar, pd, dic      
      
      common /vyi/yi
      common /vzi/zi
      common /vxij/xij
      
      common /vwi/wi
      common /vtaui/taui
      common /vvnu/vnu

      common /vbetaa/betaa
      
      common /valphaa/alphaa
      common /vazeta/azeta
      common /vagamma/agamma
      
      common /vgj/gj
      common /vweight/weight
      common /vtheta/theta
            
      common /vsigma_betaa/sigma_betaa
      common /vsigma_theta/sigma_theta
      common /va0/a0
      common /vb0/b0
      
      common /vprob/prob
            
      external gen_dic
            
    	open(unit = 5, file = 'nAHSData.txt') 
     	open(unit = 6, file = 'nUniqueDuration.txt') 
     	open(unit = 7, file = 'PesticideF_Initial.txt') 
     
    	open(unit = 11, file = 'PesticideF_Output1.txt') 
    	open(unit = 12, file = 'PesticideF_Output2.txt') 
    	open(unit = 13, file = 'PesticideF_Output3.txt') 
    	open(unit = 14, file = 'PesticideF_Output4.txt') 
    	open(unit = 15, file = 'PesticideF_Output5.txt') 
    	open(unit = 16, file = 'PesticideF_Output6.txt') 
      
    	open(unit = 21, file = 'PesticideF_Output7.txt') 
    	open(unit = 22, file = 'PesticideF_Output8.txt') 
    	open(unit = 23, file = 'PesticideF_Output9.txt') 
    	open(unit = 24, file = 'PesticideF_Output10.txt') 
    	open(unit = 25, file = 'PesticideF_Output11.txt') 
    	open(unit = 26, file = 'PesticideF_Output12.txt') 

    	open(unit = 30, file = 'PesticideF_Output13.txt') 
      
      iseed = 9999999
                
      do jj = 1, nc
          zmean(jj) = 0.d0; zstd(jj) = 0.d0
      enddo   
	do ii = 1, ni
		             
	    read(5,*) i, cancer, base, durat
                    
	    yi(i) = cancer
          
          zi(i,1) = 1.d0
          do jj = 2, nc
              zi(i,jj) = base(jj-1)
              zmean(jj) = zmean(jj) + zi(i,jj)/dfloat(ni)
              zstd(jj) = zstd(jj) + zi(i,jj)**2
          enddo
          
          do j = 1, nj
              xij(i,j,1) = durat(j)
              xij(i,j,2) = durat(j)**2
              xij(i,j,3) = durat(j)**3
          enddo
                    
      enddo 
               
      do jj = 2, nc  
          temp = (zstd(jj) - dfloat(ni)*zmean(jj)**2)
     +            /dfloat(ni-1)
          zstd(jj) = dsqrt(temp)
      enddo

      do i = 1, ni
          do jj = 2, nc  
              zi(i,jj) = (zi(i,jj) - zmean(jj))/zstd(jj)
          enddo
      enddo

      do j = 1, nj
          do k = 1, nk
              do jj = 1, nb
                  xstar(j,k,jj) = -999.d0
              enddo
          enddo
      enddo      
      do ii = 1, 5050
          
          read(6,*) j, k, ddurat
                                 
          xstar(j,k,1) = ddurat
          xstar(j,k,2) = ddurat**2
          xstar(j,k,3) = ddurat**3
                                  
      enddo

	do j = 1, nj
	    jk(j) = 0
	    do k = 1, nk
              if (xstar(j,k,1) .ne. -999.d0) then
                  jk(j) = jk(j) + 1
              endif
          enddo
      enddo    
      
c     set hyper-parameters
                  
      sigma_betaa = 1000.d0
      sigma_theta = 1000.d0
      a0 = 1.0d0 ; b0 = 0.1d0
      vnu = 7.d0 
                       
c     set initial values    
                      
      read(7,*) betaa
      read(7,*) alphaa
      read(7,*) azeta
      read(7,*) agamma
      read(7,*) theta
      
      sumw = 0.d0
      do l = 1, nl
          sumw = sumw + dexp(theta(l))
      enddo
      do l = 1, nl
          weight(l) = dexp(theta(l))/sumw
      enddo
                 
      do i = 1, ni
          wi(i) = 0.d0
          taui(i) = 1.d0
      enddo                     
      
      do j = 1, nj
          gj(j) = 1.d0
      enddo

c     Run Gibbs
      call rnset(iseed)
                  
      do ir = 1, iiniter
          write(*,*) ir
          call gibbs_gj(iseed) 
          call gibbs_wi(iseed)          
          call gibbs_taui(iseed)          
      enddo
                         
      icount = 0
      do ir = 1, initer
              
          call gibbs(iseed)

          do l = 1, nl
              ngj(l) = 0
              do j = 1, nj
                  if (nint(gj(j)) .eq. l) then
                      ngj(l) = ngj(l) + 1
                  endif
              enddo
          enddo
                            
          write(*,1) ir, ngj
          write(*,7) ir, 1, betaa
          do l = 2, nl
              write(*,7) ir, l, (alphaa(jj,l),jj=1,nb)
          enddo
          write(*,7) ir, 1, weight
          write(*,7) ir, 1, azeta, agamma
          write(*,*)
                    
          if (ir - ir/5*5 .eq. 1) then
        
              icount = icount + 1
              
              write(11,1) icount, ngj
              write(12,3) icount, betaa
              write(13,3) icount, alphaa, azeta, agamma
              write(14,3) icount, weight
              write(15,3) icount, theta
              write(16,3) icount, gj
                                      
          endif
                                   
      enddo    	  
                  
      do j = 1, nj
          
          do k = 1, nk
              
              mean_mu(j,k) = 0.d0
              std_mu(j,k) = 0.d0
              
              do l = 1, nl
                  mean_mul(j,k,l) = 0.d0                  
              enddo
              
          enddo
          
          do l = 1, nl
              
              mean_prob(j,l) = 0.d0
              std_prob(j,l) = 0.d0
              
              cngj(j,l) = 0.d0

          enddo
          
      enddo                        
            
      do jj = 1, nc
	    mean_betaa(jj) = 0.d0          
      enddo
      do j = 1, nj
          do jj = 1, nb
	        mean_delta(j,jj) = 0.d0          
          enddo
      enddo
             
      bardic = 0.0d0      
      
      icount = 0
      do ir = 1, niter
              
          call gibbs(iseed)

          do l = 1, nl
              ngj(l) = 0
              do j = 1, nj
                  if (nint(gj(j)) .eq. l) then
                      ngj(l) = ngj(l) + 1
                  endif
              enddo
          enddo
                            
          write(*,1) ir, ngj
          write(*,7) ir, 1, betaa
          do l = 2, nl
              write(*,7) ir, l, (alphaa(jj,l),jj=1,nb)
          enddo
          write(*,7) ir, 1, weight
          write(*,7) ir, 1, azeta, agamma
          write(*,*)
                    
          if (ir - ir/5*5 .eq. 1) then
        
              icount = icount + 1
              
              write(11,1) icount, ngj
              write(12,3) icount, betaa
              write(13,3) icount, alphaa, azeta, agamma
              write(14,3) icount, weight
              write(15,3) icount, theta
              write(16,3) icount, gj
                        
              do j = 1, nj
                  
                  do k = 1, jk(j)
                      
                      do l = 1, nl
                          
                          xa = 0.d0
                          do jj = 1, nb
                              xa = xa + xstar(j,k,jj)*alphaa(jj,l)
                          enddo                          
                                                    
                          mean_mul(j,k,l) = mean_mul(j,k,l) + xa
                              
                          if (l .eq. nint(gj(j))) then
                              
                              mean_mu(j,k) = mean_mu(j,k) + xa
                              std_mu(j,k) = std_mu(j,k) + xa**2
                              seq_mu(icount,j,k) = xa
                              
                          endif
                          
                      enddo
                                            
                  enddo
                                         
                  do l = 1, nl
                      
                      mean_prob(j,l) = mean_prob(j,l) + prob(j,l)
                      std_prob(j,l) = std_prob(j,l) + prob(j,l)**2
                      seq_prob(icount,j,l) = prob(j,l)
                      
                      if (l .eq. nint(gj(j))) then
                          cngj(j,l) = cngj(j,l) + 1.d0/dfloat(niter1)
                      endif

                  enddo
                  
              enddo
              
              rewind (unit = 21)
              rewind (unit = 22)
              do j = 1, nj
                  
                  do k = 1, jk(j)
                         
                      ddurat = xstar(j,k,1)
                      
                      write(21,5) icount, j, k, ddurat, 
     +                 mean_mu(j,k)/dfloat(icount), 
     +                (mean_mul(j,k,l)/dfloat(icount),l = 1,nl)
                                            
                  enddo
                                
                  write(22,6) icount, j, 
     +                        (mean_prob(j,l)/dfloat(icount),l=1,nl)
                  
              enddo              
                        
              do jj = 1, nc
	            mean_betaa(jj) = mean_betaa(jj) 
     +                           + betaa(jj)/dfloat(niter1)
              enddo
              
              do j = 1, nj
                  
                  l = nint(gj(j))
                  
                  do jj = 1, nb
                      
                      delta(j,jj) = alphaa(jj,l)
                      
                      mean_delta(j,jj) = mean_delta(j,jj) 
     +                                 + delta(j,jj)/dfloat(niter1)
                      
                  enddo
                  
              enddo
            
              bardic = bardic 
     +               - 2.d0*gen_dic(betaa,delta)/dfloat(niter1)
              
          endif
                                   
      enddo    	  
      
      call rnget(iseed)
                        
      do j = 1, nj
      do k = 1, nk
              
          mean_mu(j,k) = mean_mu(j,k)/dfloat(icount)
              
          temp = (std_mu(j,k) 
     +            - dfloat(icount)*mean_mu(j,k)**2)
     +            /dfloat(icount-1)
		std_mu(j,k) = dsqrt(temp)
          
          do l = 1, nl
              mean_mul(j,k,l) = mean_mul(j,k,l)/dfloat(icount) 
          enddo
                                  
      enddo
      enddo

      do j = 1, nj
          do l = 1, nl
              
              mean_prob(j,l) = mean_prob(j,l)/dfloat(icount) 
              
              temp = (std_prob(j,l)
     +                - dfloat(icount)*mean_prob(j,l)**2)
     +                /dfloat(icount-1)
		    std_prob(j,l) = dsqrt(temp)
                            
          enddo
      enddo
      
      do j = 1, nj
      do k = 1, nk
              
          do ir = 1, icount
              ahpd(ir) = seq_mu(ir,j,k) 
          enddo
          call hpd(icount, 0.05d0, ahpd, alow, aupp)
          lower_mu(j,k) = alow(1)
          upper_mu(j,k) = aupp(1)
                    
      enddo
      enddo

      do j = 1, nj
          do l = 1, nl
              
              do ir = 1, icount
                  ahpd(ir) = seq_prob(ir,j,l)
              enddo
              call hpd(icount, 0.05d0, ahpd, alow, aupp)
              lower_prob(j,l) = alow(1)
              upper_prob(j,l) = aupp(1)
                                          
          enddo
      enddo
      
      do j = 1, nj
                  
          do k = 1, jk(j)
                         
              ddurat = xstar(j,k,1)
                      
              write(23,6) j, k, ddurat, 
     +                    mean_mu(j,k), std_mu(j,k), 
     +                    lower_mu(j,k), upper_mu(j,k)
              
              write(24,6) j, k, ddurat, 
     +                    (mean_mul(j,k,l), l = 1, nl)
              
          enddo

          write(25,3) j,(mean_prob(j,l), std_prob(j,l), 
     +                   lower_prob(j,l), upper_prob(j,l),
     +                   l = 1, nl)
              
          write(26,3) j, (cngj(j,l),l = 1, nl)
          
      enddo              
      
      dicbar = -2.d0*gen_dic(mean_betaa,mean_delta)
      pd = bardic - dicbar
      dic = dicbar + 2.0d0*pd 

      write(30,8) dicbar, pd, dic
                         
    1 format(1000i5)
    2 format(i5,1000f10.5)
    3 format(i5,1000f20.10)
    4 format(i5,1000f10.5)
    5 format(3i5,1000f20.10)
    6 format(2i5,1000f20.10)
    7 format(2i5,1000f10.5)
    8 format(1000f20.10)
            
      end program
               
      
      real*8 function gen_dic(betaa,delta)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2
		
	real*8 yi(ni), zi(ni,nc), xij(ni,nj,nb)    
      real*8 vnu, betaa(nc), delta(nj,nb)      
      
      real*8 zb, xa, fxij, cdf, dic
      
      common /vyi/yi
      common /vzi/zi
      common /vxij/xij
      
      common /vvnu/vnu
      
      external dtdf
                        
      dic = 0.d0
      do i = 1, ni
            
          zb = 0.d0
          do jj = 1, nc
              zb = zb + zi(i,jj)*betaa(jj)
          enddo
                    
          fxij = 0.d0               
          do j = 1, nj  
              if (xij(i,j,1) .ne. 0.d0) then
                                                  
                  xa = 0.d0
                  do jj = 1, nb
                      xa = xa + xij(i,j,jj)*delta(j,jj)
                  enddo
                      
                  fxij = fxij + xa
                  
              endif
          enddo
          
          cdf = dtdf(zb + fxij, vnu)
                  
          if (yi(i) .eq. 1.d0) then
                                            
              dic = dic + dlog(cdf)
                      
          else

              dic = dic + dlog(1.d0 - cdf)
                      
          endif
              
      enddo
                        
	gen_dic = dic

	end function          
      
 	include 'PesticideF_Gibbs.f'
	include 'tnorm.f'
	include 'optim1.f'
      include 'rGiGDist.f'                        
      include 'hpd.f'                        
            