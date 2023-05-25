      real*8 function rgenGIG(p,a,b,iseed)
      implicit real*8 (a-h, o-z)  
      
c     generate the random number from generalized inverse gaussian distribution 
c     f(x) = (a/b)^(p/2) / (2 K_p (sqrt(ab))) x^(p-1) exp( -(ax + b/x) /2 ), 
c     x > 0, a > 0, and b > 0    
c     We will only consider the two-parameter version, because
c     if X follows the three-parameter law with parameters (p, a, b), then
c     sqrt(a/b)X follows the two parameter law with tau = p and w = sqrt(ab)
c     so, gig = sqrt(b/a)X
      
      real*8 p, a, b
      
      real*8 tau, omega, alpha
      real*8 xx, tt, ss, temp1, temp2

      real*8 eta, zeta, theta, xi
      real*8 pp, rr, td, sd, qq
      real*8 uu, vv, ww, f1, f2, rv
      
c     Setup -- we sample from the two parameter version of the GIG(alpha,omega)
      tau = p
      omega = dsqrt(a*b)
      alpha = dsqrt(omega**2 + tau**2) - tau
            
c     Find t
      xx = -fgig(1.d0, alpha, tau)
      
      if ( (xx .ge. 0.5d0) .and. (xx .le. 2.d0) ) then
          
          tt = 1.d0
          
      else if (xx .gt. 2.d0) then
          
          tt = dsqrt(2.d0/(alpha + tau))
          
      else if (xx .lt. 0.5d0) then
          
          tt = dlog(4.d0/(alpha + 2.d0*tau))
          
      endif
      
c     Find s
      xx = -fgig(-1.d0, alpha, tau)
      
      if ( (xx .ge. 0.5d0) .and. (xx .le. 2.d0) ) then
          
          ss = 1.d0
          
      else if (xx .gt. 2.d0) then
          
          ss = dsqrt(4.d0/(alpha*dcosh(1.d0) + tau))
      
      else if (xx .lt. 0.5d0) then
          
          temp1 = 1.d0/tau
          temp2 = dlog( 1.d0 + 1.d0/alpha 
     +                  + dsqrt(1.d0/alpha**2 + 2.d0/alpha) )
          
          if (temp1 .le. temp2) then
              
              ss = temp1
          
          else
              
              ss = temp2
              
          endif
          
      endif
             
c     Generation
      
      eta = -fgig(tt, alpha, tau)
      
      zeta = -dfgig(tt, alpha, tau)
      
      theta = -fgig(-ss, alpha, tau)
      
      xi = dfgig(-ss, alpha, tau)
      
      pp = 1.d0/xi
      rr = 1.d0/zeta
      td = tt - rr*eta
      sd = ss - pp*theta
      qq = td + sd
                  
      iiflag = 1
      xx = 0.d0
      do while (iiflag .eq. 1) 
                
	    call rnset(iseed)
	    uu = drnunf()
	    call rnget(iseed)

	    call rnset(iseed)
	    vv = drnunf()
	    call rnget(iseed)

	    call rnset(iseed)
	    ww = drnunf()
	    call rnget(iseed)
      
          if ( uu .lt. qq/(pp + qq + rr) ) then
          
              xx = -sd + qq*vv
            
          else if ( uu .lt. (qq + rr)/(pp + qq + rr) ) then
          
              xx = td - rr*dlog(vv)
            
          else
          
              xx = -sd + pp*dlog(vv)
          
          endif

          f1 = dexp(-eta - zeta*(xx - tt))
          f2 = dexp(-theta + xi*(xx + ss))
          
          temp1 = ww*ffgig(xx, sd, td, f1, f2)
          temp2 = fgig(xx, alpha, tau)
      
          if ( temp1 .le. dexp(temp2) ) then          
              
              iiflag = 0
              
          endif     
      
      enddo
           
c     Transform X back to the three parameter GIG(p,a,b)
     
      rv = dexp(xx)*(tau/omega 
     +   + dsqrt(1.d0 + (tau/omega)**2))
      
      rgenGIG = rv*dsqrt(b/a)
            
      end function     
            
      
      real*8 function fgig(x, alpha, tau)
      implicit real*8 (a-h, o-z) 
      
      real*8 x, alpha, tau, val
      
      val = -alpha*(dcosh(x) - 1.d0) - tau*(dexp(x) - x - 1.d0)
      
      fgig = val
      
      end function          
      
      
      real*8 function dfgig(x, alpha, tau)
      implicit real*8 (a-h, o-z) 
      
      real*8 x, alpha, tau, val
      
      val = -alpha*dsinh(x) - tau*(dexp(x) - 1.d0)
      
      dfgig = val
      
      end function   
      
           
      real*8 function ffgig(x, sd, td, f1, f2)
      implicit real*8 (a-h, o-z) 
      
      real*8 x, sd, td, f1, f2
      real*8 a, b, c, val
      
      a = 0.d0
      b = 0.d0
      c = 0.d0
      
      if ( (x .ge. -sd) .and. (x .le. td) ) then
          
          a = 1.d0
          
      else if (x .gt. td) then
          
          b = f1
          
      else if (x .lt. -sd) then
          
          c = f2
          
      endif

      val = a + b + c
      
      ffgig = val
      
	end function          
