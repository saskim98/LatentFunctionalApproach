	real*8 function ytuvn(a,b,la,lb,iseed)
c  this function generates a n(0,1) random variable
c  subject to the constraint that it be in an interval
c  (a,b), where the endpoints may be finite or infinite.
c  references:  j. geweke, 1991.  "efficient simulation from the multivariate
c               normal and student-t distributions subject to linear
c               constraints," in e.m. keramidas (ed.), computing science and
c               statistics: proceedings of the 23rd symposium on the
c               interface, pp. 571-578.  fairfax, va: interface foundation
c               of north america, inc.
c  update history:  develped in gibbs/tnorm.  brought into ylib 1/23/92
c                   modified to allow singleton 2/20/92
c                   modified to prevent overflow before 10, 3/31/92
c  inputs:
c  a, b    endpoints of interval; a < b if la = lb = .false.
c  la      .true. if left endpoint is - infinity; in this
c          case a is ignored.
c  lb      .true. if right endpoint is + infinity; in this
c          case b is ignored.
c  output:
c  ggrnrm  random variable
      implicit real*8 (a-h,o-z)
      logical la,lb,lflip
      data eps,t1,t2,t3,t4/2.0d0,.375d0,2.18d0,.725d0,.45d0/
      
      external drnunf
      
      f(x)=dexp(-0.5d0*x**2)
      if(la.and.lb)go to 160
      lflip=.false.
      if(la.or.lb)go to 100
      if(b.le.a)go to 170
c ******* finite interval
      c1=a
      c2=b
      if((c1*c2).gt.0.0d0)go to 30
c ++++ (a,b) includes 0
      if((c1.gt.-t1).and.(c2.lt.t1))go to 20
c -- f(a) or f(b) small: full normal with rejection
   10 x=drnnof()
      if(x.lt.c1)go to 10
      if(x.gt.c2)go to 10
      go to 150
c -- f(a) and f(b) large: uniform importance sampling
   20 cdel=c2-c1
   25 x=c1+cdel*drnunf()
      if(drnunf().gt.f(x))go to 25
      go to 150
c ++++ (a,b) excludes 0
c -- transform to both positive
   30 if(c1.gt.0.0d0)go to 40
      c=c1
      c1=-c2
      c2=-c
      lflip=.true.
   40 f1=f(c1)
      f2=f(c2)
      if(f2.lt.eps)go to 60
      if((f1/f2).gt.t2)go to 60
c  -- f(a)/f(b) not large: uniform importance sampling
   50 cdel=c2-c1
   55 x=c1+cdel*drnunf()
      if(drnunf().gt.(f(x)/f1))go to 55
      go to 140
   60 if(c1.gt.t3)go to 80
c -- p(x>a) and f(a)/f(b) large: half-normal with rejection
   70 x=dabs(drnnof())
      if(x.lt.c1)go to 70
      if(x.gt.c2)go to 70
      go to 140
c -- p(x>a) small, f(a)/f(b) large: exponential importance
c    sampling with rejection
   80 c=c2-c1
   90 z=-dlog(drnunf())/c1
      if(z.gt.c)go to 90
      if(drnunf().gt.f(z))go to 90
      x=c1+z
      go to 140
c ****** half-line interval
  100 c1=a
c -- transform to bound from below if a = -infinity
      if(lb)go to 110
      c1=-b
      lflip=.true.
  110 if(c1.gt.t4)go to 130
c -- a not large: full normal with rejection
  120 x=drnnof()
      if(x.lt.c1)go to 120
      go to 140
c -- a small: exponential importance sampling
  130 z=-dlog(drnunf())/c1
      if(drnunf().gt.f(z))go to 130
      x=c1+z
  140 if(lflip)x=-x
  150 ytuvn=x
      return
c ****** whole interval
  160 ytuvn=drnnof()
      return
c  ***** singleton
  170 ytuvn=a

      end function