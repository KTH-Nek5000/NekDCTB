!> @file sprj_tool.f
!! @ingroup nekamr
!! @brief Set of routines for projection of point on a surface. I use
!!  routines from Numerical Recipes ch 10 to minimalise distance between
!!  given point in the space and the point on the surface. The distance to
!!  surface has to be defined in the sfunc function and the derivatives
!!  in dsfunc. Two different methods are implemented: Powell' (does not require
!!  function gradient) and Fletcher-Reeves (uses function gradient).
!!  I use both as the gradient is not always well defined.
!!  Few modifications added with respect to the original code. 
!! @author Adam Peplinski
!! @date MAR 23, 2017
!=======================================================================
!> @brief Project face onto 2D surface
!! @param[inout] xe       point position/position shift
!! @param[in]    ndim     space dimension (must be defined by parameter)
!! @param[in]    np       number of points (must be defined by parameter)
!! @param[in]    pmin     starting points list
!! @param[in]    nmin     number of starting points
!! @param[in]    ftol     convergance tolerance  
!! @param[in]    ifderiv  do we use derivatives
!! @param[inout] fret     minimal distance
!! @param[inout] ortho    orthogonality check
!! @param[in]    rpos     subroutine returning surface pos. for given parameter
!! @param[in]    rdpos    subroutine returning derivatives for given parameter
!! @param[in]    rbnd     subroutine returning parameter boundaries
      subroutine fcs_prj(xe,ndim,np,pmin,nmin,ftol,ifderiv,fret,ortho,
     $     rpos,rdpos,rbnd)
      implicit none

      ! argument list
      integer np, ndim, nmin
      real xe(ndim,np,np)
      real pmin(2,nmin)
      real ftol, fret, ortho
      logical ifderiv
      external rpos,rdpos,rbnd

      ! local variables
      integer nm
      parameter (nm=2)
      integer max_it
      parameter (max_it=20)
      integer il, jl, kl
      real xs(ndim), xstmp(ndim), xel(ndim,np,np)
      real pl(nm), fretl, orthol
      real ptmp(nm)
!-----------------------------------------------------------------------
      ! test initial parameter value
      ! point in the middle of the edge
      do il=1,ndim
         xstmp(il) = 0.25*(xe(il,1,1)+xe(il,1,np)+xe(il,np,1)+
     $        xe(il,np,np))
      enddo
      ! initil distance
      fret = 1.0e50
      ! set of test to  pick up correct local minima
      do jl = 1,nmin
         ! copy initial position
         do il =1,nm
            pl(il) = pmin(il,jl)
         enddo
         do il=1,ndim
            xs(il) = xstmp(il)
         enddo
         ! get local minima
         call pnt_prj(pl,nm,xs,ndim,ftol,ifderiv,fretl,orthol,
     $        rpos,rdpos,rbnd)
         ! is it global minima
         if (fretl.lt.fret) then
            fret = fretl
            ortho = orthol
            do il =1,nm
               ptmp(il) = pl(il)
            enddo
         endif
      enddo

      ! get initial points position
      fret = 0.0
      ortho = 0.0
      do kl=1,np
         do il =1,nm
            pl(il) = ptmp(il)
         enddo
         do jl=1,np
            do il=1,ndim
               xs(il) = xe(il,jl,kl)
            enddo
            call pnt_prj(pl,nm,xs,ndim,ftol,ifderiv,fretl,orthol,
     $           rpos,rdpos,rbnd)
            ! store projected point position and parameter value
            do il=1,ndim
               xel(il,jl,kl) = xs(il)
            enddo
            ! check max distance and inner product
            fret = max(fret,abs(fretl))
            ortho = max(ortho,abs(orthol))
            if (jl.eq.1) then
               do il =1,nm
                  ptmp(il) = pl(il)
               enddo
            endif
         enddo
      enddo

      ! get shift
      do kl=1,np
         do jl=1,np
            do il=1,ndim
               xe(il,jl,kl) = xel(il,jl,kl) - xe(il,jl,kl)
            enddo
         enddo
      enddo
      
      return
      end subroutine
!=======================================================================
!> @brief Project edge onto 2D surface
!! @param[inout] xe       point position/position shift
!! @param[in]    ndim     space dimension (must be defined by parameter)
!! @param[in]    np       number of points (must be defined by parameter)
!! @param[in]    pmin     starting points list
!! @param[in]    nmin     number of starting points
!! @param[in]    ftol     convergance tolerance  
!! @param[in]    ifderiv  do we use derivatives
!! @param[inout] fret     minimal distance
!! @param[inout] ortho    orthogonality check
!! @param[in]    rpos     subroutine returning surface pos. for given parameter
!! @param[in]    rdpos    subroutine returning derivatives for given parameter
!! @param[in]    rbnd     subroutine returning parameter boundaries
      subroutine edg2D_prj(xe,ndim,np,pmin,nmin,ftol,ifderiv,fret,ortho,
     $     rpos,rdpos,rbnd)
      implicit none

      ! argument list
      integer np, ndim, nmin
      real xe(ndim,np)
      real pmin(2,nmin)
      real ftol, fret, ortho
      logical ifderiv
      external rpos,rdpos,rbnd
      ! local variables
      integer nm
      parameter (nm=2)
      integer max_it
      parameter (max_it=20)
      integer il, jl, kl
      real xs(ndim), xstmp(ndim), xel(ndim,np)
      real pl(nm), fretl, orthol
      real ptmp(nm)
!-----------------------------------------------------------------------
      ! test initial parameter value
      ! point in the middle of the edge
      do il=1,ndim
         xstmp(il) = 0.5*(xe(il,1)+xe(il,np))
      enddo
      ! initil distance
      fret = 1.0e50
      ! set of test to  pick up correct local minima
      do jl = 1,nmin
         ! copy initial position
         do il =1,nm
            pl(il) = pmin(il,jl)
         enddo
         do il=1,ndim
            xs(il) = xstmp(il)
         enddo
         ! get local minima
         call pnt_prj(pl,nm,xs,ndim,ftol,ifderiv,fretl,orthol,
     $        rpos,rdpos,rbnd)
         ! is it global minima
         if (fretl.lt.fret) then
            fret = fretl
            ortho = orthol
            do il =1,nm
               ptmp(il) = pl(il)
            enddo
         endif
      enddo

      ! get initial points position
      fret = 0.0
      ortho = 0.0
      do il =1,nm
         pl(il) = ptmp(il)
      enddo
      do jl=1,np
         do il=1,ndim
            xs(il) = xe(il,jl)
         enddo
         call pnt_prj(pl,nm,xs,ndim,ftol,ifderiv,fretl,orthol,
     $        rpos,rdpos,rbnd)
         ! store projected point position and parameter value
         do il=1,ndim
            xel(il,jl) = xs(il)
         enddo
         ! check max distance and inner product
         fret = max(fret,abs(fretl))
         ortho = max(ortho,abs(orthol))
      enddo

      ! get shift
      do il=1,np
         do jl=1,ndim
            xe(jl,il) = xel(jl,il) - xe(jl,il)
         enddo
      enddo
      
      return
      end subroutine
!=======================================================================
!> @brief Project edge onto 1D curve
!! @param[inout] xe       point position/position shift
!! @param[in]    xv       point distribution (along the arc)
!! @param[in]    ndim     space dimension (must be defined by parameter)
!! @param[in]    np       number of points (must be defined by parameter)
!! @param[in]    pmin     starting points list
!! @param[in]    nmin     number of starting points
!! @param[in]    ftol     convergance tolerance  
!! @param[in]    ifderiv  do we use derivatives
!! @param[inout] fret     minimal distance
!! @param[inout] ortho    orthogonality check
!! @param[in]    rpos     subroutine returning surface pos. for given parameter
!! @param[in]    rdpos    subroutine returning derivatives for given parameter
!! @param[in]    rbnd     subroutine returning parameter boundaries
      subroutine edg1D_prj(xe,xv,ndim,np,pmin,nmin,ftol,ifderiv,fret,
     $     ortho,rpos,rdpos,rbnd)
      implicit none

      ! argument list
      integer np, ndim, nmin
      real xe(ndim,np), xv(np)
      real pmin(nmin)
      real ftol, fret, ortho
      logical ifderiv
      external rpos,rdpos,rbnd

      ! local variables
      integer nm
      parameter (nm=1)
      integer max_it
      parameter (max_it=20)
      integer il, jl, kl
      real xs(ndim), xstmp(ndim), xel(ndim,np)
      real pval(np)
      real pl(nm), fretl, orthol
      real ptmp(nm)
      integer nit, itest
      real arcl, arclo, dp, da
!-----------------------------------------------------------------------  
      ! test initial parameter value
      ! point in the middle of the edge
      do il=1,ndim
         xstmp(il) = 0.5*(xe(il,1)+xe(il,np))
      enddo
      ! initil distance
      fret = 1.0e50
      !set of test to  pick up correct local minima
      do jl = 1,nmin
         ! copy initial position
         pl(1) = pmin(jl)
         do il=1,ndim
            xs(il) = xstmp(il)
         enddo
         ! get local minima
         call pnt_prj(pl,nm,xs,ndim,ftol,ifderiv,fretl,orthol,
     $        rpos,rdpos,rbnd)
         ! is it global minima
         if (fretl.lt.fret) then
            fret = fretl
            ptmp(1) = pl(1)
         endif
      enddo

      ! get first vertex position
      fret = 0.0
      ortho = 0.0
      pl(1) = ptmp(1)
      do il=1,ndim
         xs(il) = xe(il,1)
      enddo
      call pnt_prj(pl,nm,xs,ndim,ftol,ifderiv,fretl,orthol,
     $             rpos,rdpos,rbnd)
      ! store projected point position and parameter value
      pval(1) = pl(1)
      do il=1,ndim
         xel(il,1) = xs(il)
      enddo
      ! check max distance and inner product
      fret = max(fret,abs(fretl))
      ortho = max(ortho,abs(orthol))

      ! get second vertex position
      pl(1) = ptmp(1)
      do il=1,ndim
         xs(il) = xe(il,np)
      enddo
      call pnt_prj(pl,nm,xs,ndim,ftol,ifderiv,fretl,orthol,
     $             rpos,rdpos,rbnd)
      ! store projected point position and parameter value
      pval(np) = pl(1)
      do il=1,ndim
         xel(il,np) = xs(il)
      enddo
      ! check max distance and inner product
      fret = max(fret,abs(fretl))
      ortho = max(ortho,abs(orthol))

      if (np.eq.2) then
         ! get shift
         do il=1,np
            do jl=1,ndim
               xe(jl,il) = xel(jl,il) - xe(jl,il)
            enddo
         enddo
         return
      endif

      ! get arc length
      arclo = 0.0
      nit = 10
      do il=1,max_it
         nit = nit*2
         dp = (pval(np)-pval(1))/real(nit)
         arcl = 0.0
         pl(1) = pval(1)
         do jl=1,ndim
            xstmp(jl) = xel(jl,1)
         enddo
         do jl=1,nit
            pl(1) = pl(1) + dp
            call rpos(xs,ndim,pl,nm)
            da = 0.0
            do kl=1,ndim
               da = da +(xstmp(kl)-xs(kl))**2
               xstmp(kl) = xs(kl)
            enddo
            arcl = arcl + sqrt(da)
         enddo
         if (abs(arcl-arclo).le.ftol) exit
         arclo = arcl
      enddo

      ! get parameter and point position
      arclo = 0.0
      pl(1) = pval(1)
      do jl=1,ndim
         xstmp(jl) = xel(jl,1)
      enddo
      do il = 2,np
         do jl=1,nit
            pl(1) = pl(1) + dp
            call rpos(xs,ndim,pl,nm)
            da = 0.0
            do kl=1,ndim
               da = da +(xstmp(kl)-xs(kl))**2
               xstmp(kl) = xs(kl)
            enddo
            arclo = arclo + sqrt(da)
            if (arclo.gt.(arcl*xv(il))) then
               ptmp(1) = pl(1) - dp*(arclo-arcl*xv(il))/sqrt(da)
               call rpos(xs,ndim,ptmp,nm)
               do kl=1,ndim
                  xel(kl,il) = xs(kl)
               enddo
               exit
            endif
         enddo
      enddo

      ! get shift
      do il=1,np
         do jl=1,ndim
            xe(jl,il) = xel(jl,il) - xe(jl,il)
         enddo
      enddo
      
      return
      end subroutine
!=======================================================================
!> @brief Project point onto the surface
!! @param[inout] p        starting point/min. location
!! @param[in]    n        number of variables
!! @param[inout] xs       point position
!! @param[in]    ndim     space dimension (must be defined by parameter)
!! @param[in]    ftol     convergance tolerance  
!! @param[in]    ifderiv  do we use derivatives
!! @param[out]   fret     minimal distance
!! @param[out]   ortho    orthogonality check
!! @param[in]    rpos     subroutine returning surface pos. for given parameter
!! @param[in]    rdpos    subroutine returning derivatives for given parameter
!! @param[in]    rbnd     subroutine returning parameter boundaries
      subroutine pnt_prj(p,n,xs,ndim,ftol,ifderiv,fret,ortho,
     $     rpos,rdpos,rbnd)
      implicit none

      ! argument list
      integer n, ndim
      real p(n), xs(ndim)
      real ftol, fret, ortho
      logical ifderiv
      external rpos,rdpos,rbnd

      ! local variables
      real x0(ndim)
      integer il, jl, kl
      integer iter
      real xi(ndim,ndim), dx(n)
!-----------------------------------------------------------------------
      ! set point for calculation
      do il=1,ndim
         x0(il) = xs(il)
      enddo
      ! find minimum
      if (ifderiv) then
         call frprmn(p,n,ftol,iter,fret,x0,ndim,rpos,rdpos,rbnd)
      else
         do jl=1,ndim
            do kl=1,ndim
               xi(jl,kl)=0.0
            enddo
            xi(jl,jl)=1.0
         enddo
         call powell(p,xi,n,ndim,ftol,iter,fret,x0,ndim,rpos,rdpos,rbnd)
      endif
      ! get position
      call rpos(xs,ndim,p,n)
      ! check orthogonality
      call dsfunc(p,dx,n,x0,ndim,rpos,rdpos)
      ortho = 0.0
      do il=1,n
         ortho = ortho + dx(il)
      enddo
      ! place to correct orthogonality
      return
      end subroutine
!=======================================================================
!> @brief Powell minimalisation of n-dimensional function; Numerical Recipes 10.5
!! @param[inout] p     starting point/min. location
!! @param[inout] xi    matrix of directions
!! @param[in]    n     number of variables
!! @param[in]    np    array sizes
!! @param[in]    ftol  convergance tolerance  
!! @param[out]   iter  iteration count
!! @param[out]   fret  minimal value of the function    
!! @param[in]    x0       reference point
!! @param[in]    ndim     space dimension
!! @param[in]    rpos     subroutine returning surface pos. for given parameter
!! @param[in]    rdpos    subroutine returning derivatives for given parameter
!! @param[in]    rbnd     subroutine returning parameter boundaries
      subroutine powell(p,xi,n,np,ftol,iter,fret,x0,ndim,
     $                  rpos,rdpos,rbnd)
      implicit none

      ! argument list
      integer iter,n,np,ndim
      real p(n),xi(np,np),fret,ftol,x0(ndim)
      external rpos,rdpos,rbnd

      ! local vairables
      integer nmax,itmax
      parameter (nmax=10,itmax=200)
      real tiny
      parameter (tiny=1.d-30)
      integer i,ibig,j
      real del,fp,fptt,t,pt(nmax),ptt(nmax),xit(nmax)
      logical ifderiv
      parameter (ifderiv=.FALSE.)

      ! functions
      real sfunc
!-----------------------------------------------------------------------
      fret=sfunc(p,n,x0,ndim,rpos)
      do j=1,n
         pt(j)=p(j)
      enddo
      iter=0
 1    iter=iter+1
      fp=fret
      ibig=0
      del=0.
      do i=1,n
         do j=1,n
            xit(j)=xi(j,i)
         enddo
         fptt=fret
         call linmin(p,xit,n,ftol,fret,ifderiv,x0,ndim,rpos,rdpos,rbnd)
         if(fptt-fret.gt.del)then
            del=fptt-fret
            ibig=i
         endif
      enddo
      if(2.*(fp-fret).le.ftol*(abs(fp)+abs(fret))+tiny) return
      if(iter.eq.itmax)
     $     write(*,*) 'Error: powell; exceeding maximum iterations'
      do j=1,n
         ptt(j)=2.*p(j)-pt(j)
         xit(j)=p(j)-pt(j)
         pt(j)=p(j)
      enddo
      fptt=sfunc(ptt,n,x0,ndim,rpos)
      if(fptt.ge.fp)goto 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0.)goto 1
      call linmin(p,xit,n,ftol,fret,ifderiv,x0,ndim,rpos,rdpos,rbnd)
      do j=1,n
         xi(j,ibig)=xi(j,n)
         xi(j,n)=xit(j)
      enddo
      goto 1
      end subroutine
!=======================================================================
!> @brief Fletcher-Reeves minimalisation of n-dimensional function; Numerical Recipes 10.6
!! @param[inout] p     starting point/min. location
!! @param[in]    n     number of variables
!! @param[in]    ftol  convergance tolerance  
!! @param[out]   iter  iteration count
!! @param[out]   fret  minimal value of the function
!! @param[in]    x0       reference point
!! @param[in]    ndim     space dimension
!! @param[in]    rpos     subroutine returning surface pos. for given parameter
!! @param[in]    rdpos    subroutine returning derivatives for given parameter
!! @param[in]    rbnd     subroutine returning parameter boundaries
      subroutine frprmn(p,n,ftol,iter,fret,x0,ndim,rpos,rdpos,rbnd)
      implicit none

      ! argument list
      integer n,iter,ndim
      real fret,ftol,p(n),x0(ndim)
      external rpos,rdpos,rbnd

      ! local variables
      integer nmax, itmax
      parameter (nmax=10,itmax=200)
      real eps
      parameter (eps=1.d-14)
      integer its,j
      real dgg,fp,gam,gg,g(nmax),h(nmax),xi(nmax)
      logical ifderiv
      parameter (ifderiv=.TRUE.)

      ! functions
      real sfunc
!     this for testing
!      real xs(ndim)
!-----------------------------------------------------------------------
      ! initialisation
      fp=sfunc(p,n,x0,ndim,rpos)
      call dsfunc(p,xi,n,x0,ndim,rpos,rdpos)
      do j=1,n
         g(j)=-xi(j)
         h(j)=g(j)
         xi(j)=h(j)
      enddo
      ! main iteration loop
      do its=1,itmax
         iter=its
         call linmin(p,xi,n,ftol,fret,ifderiv,x0,ndim,rpos,rdpos,rbnd)
         if(2.*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+eps)) return
         fp=fret
         call dsfunc(p,xi,n,x0,ndim,rpos,rdpos)
         gg=0.
         dgg=0.
         do j=1,n
            gg=gg+g(j)**2
            dgg=dgg+(xi(j)+g(j))*xi(j)
         enddo
         if(gg.eq.0.)return
         gam=dgg/gg
         do j=1,n
            g(j)=-xi(j)
            h(j)=g(j)+gam*h(j)
            xi(j)=h(j)
         enddo
      enddo
      write(*,*) 'Warning: frprmn; maximum iterations exceeded'
!     this for testing
!      call rpos(xs,ndim,p,n)
!      write(*,*) p(:),xs(:)
!     testing end
      return
      end subroutine
!=======================================================================
!> @brief 1D function min; Numerical Recipes 10.5
!! @param[inout] p        starting point/min. location
!! @param[inout] xi       n-dimensional direction/point shift
!! @param[in]    n        number of variables
!! @param[in]    ftol     convergance tolerance
!! @param[out]   fret     minimal value of the function
!! @param[in]    ifderiv  do we use function derivatives
!! @param[in]    x0       reference point
!! @param[in]    ndim     space dimension
!! @param[in]    rpos     subroutine returning surface pos. for given parameter
!! @param[in]    rdpos    subroutine returning derivatives for given parameter
!! @param[in]    rbnd     subroutine returning parameter boundaries
      subroutine linmin(p,xi,n,ftol,fret,ifderiv,x0,ndim,
     $                  rpos,rdpos,rbnd)
      implicit none

      ! argument list
      integer n, ndim
      real ftol,fret,p(n),xi(n), x0(ndim)
      logical ifderiv
      external rpos,rdpos,rbnd

      ! local variables
      integer j,iter
      real ax,bx,fa,fb,fx,xmin,xx,ps(n),pi(n)
      real x_min, x_max
!-----------------------------------------------------------------------
      do j=1,n
         ps(j)=p(j)
         pi(j)=xi(j)
      enddo
      ax=0.
      xx=1.
      call bnd_par(p,xi,n,x_min,x_max,rbnd)
      ax=max(min(ax,x_max),x_min)
      xx=max(min(xx,x_max),x_min)
      call mnbrak(ax,xx,bx,fa,fx,fb,x_min,x_max,ps,pi,n,x0,ndim,rpos)
      if (ifderiv) then
         call dbrent(ax,xx,bx,ftol,xmin,fret,iter,ps,pi,n,x0,ndim,
     $               rpos,rdpos)
      else
         call brent(ax,xx,bx,ftol,xmin,fret,iter,ps,pi,n,x0,ndim,rpos)
      endif
      do j=1,n    
         xi(j)=xmin*xi(j)
         p(j)=p(j)+xi(j)
      enddo
      return
      end subroutine
!=======================================================================
!> @brief Bracket function min; Numerical Recipes 10.1
!! @param[inout] ax, bx, cx    minimum brackets
!! @param[out]   fa, fb, fc    function values   
!! @param[in]    x_min,x_max   min and max bracket boundaries
!! @param[in]    ps            starting point
!! @param[in]    pi            displacement vector
!! @param[in]    np            number of parameters (manifold dimension)
!! @param[in]    x0            reference point
!! @param[in]    ndim          space dimension
!! @param[in]    rpos          subroutine returning surface pos. for given parameter
      subroutine mnbrak(ax,bx,cx,fa,fb,fc,x_min,x_max,ps,pi,np,x0,ndim,
     $                  rpos)
      implicit none

      ! argument list
      real ax,bx,cx,fa,fb,fc,x_min,x_max
      integer np, ndim
      real ps(np), pi(np), x0(ndim)
      external rpos

      ! local variables
      real gold,glimit,tiny
      parameter (gold=1.618034, glimit=100., tiny=1.d-30)
      real dum,fu,q,r,u,ulim,itest

      ! functions
      real f1dim
!-----------------------------------------------------------------------
      itest=0.0
      fa=f1dim(ax,ps,pi,np,x0,ndim,rpos)
      fb=f1dim(bx,ps,pi,np,x0,ndim,rpos)
      if(fb.gt.fa)then
         dum=ax
         ax=bx
         bx=dum
         dum=fb
         fb=fa
         fa=dum
      endif
      cx=bx+gold*(bx-ax)
      cx=max(min(cx,x_max),x_min)
      if (abs(cx-bx).lt.tiny) then
         bx = 0.5*(ax+cx)
      endif
      fc=f1dim(cx,ps,pi,np,x0,ndim,rpos)
 1    if(fb.ge.fc)then
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),tiny),q-r))
         u=max(min(u,x_max),x_min)
         ulim=bx+glimit*(cx-bx)
         if((bx-u)*(u-cx).gt.0.)then
            fu=f1dim(u,ps,pi,np,x0,ndim,rpos)
            if(fu.lt.fc)then
               ax=bx
               fa=fb
               bx=u
               fb=fu
               return
            else if(fu.gt.fb)then
               cx=u
               fc=fu
               return
            endif
            u=cx+gold*(cx-bx)
            u=max(min(u,x_max),x_min)
            fu=f1dim(u,ps,pi,np,x0,ndim,rpos)
         else if((cx-u)*(u-ulim).gt.0.)then
            fu=f1dim(u,ps,pi,np,x0,ndim,rpos)
            if(fu.lt.fc)then
               bx=cx
               cx=u
               u=cx+gold*(cx-bx)
               u=max(min(u,x_max),x_min)
               fb=fc
               fc=fu
               fu=f1dim(u,ps,pi,np,x0,ndim,rpos)
            endif
         else if((u-ulim)*(ulim-cx).ge.0.)then
            u=ulim
            u=max(min(u,x_max),x_min)
            fu=f1dim(u,ps,pi,np,x0,ndim,rpos)
         else
            u=cx+gold*(cx-bx)
            u=max(min(u,x_max),x_min)
            fu=f1dim(u,ps,pi,np,x0,ndim,rpos)
         endif
         ax=bx
         bx=cx
         cx=u
         fa=fb
         fb=fc
         fc=fu
         itest=itest+1
         if (ax.eq.cx) then
            return
         endif
         goto 1
      endif
      return
      end subroutine
!=======================================================================
!> @brief 1D function min, Brent's method; Numerical Recipes 10.2
!! @param[in]   ax, bx, cx    minimum brackets
!! @param[in]   tol           tolerance
!! @param[out]  xmin          minimum position
!! @param[out]  fmin          function min
!! @param[out]  iterd         iteration count
!! @param[in]   ps            starting point
!! @param[in]   pi            displacement vector
!! @param[in]   np            number of parameters (manifold dimension)
!! @param[in]   x0            reference point
!! @param[in]   ndim          space dimension
!! @param[in]   rpos          subroutine returning surface pos. for given parameter
      subroutine brent(ax,bx,cx,tol,xmin,fmin,iterd,ps,pi,np,x0,ndim,
     $                 rpos)
      implicit none

      ! argument list
      integer iterd
      real ax,bx,cx,tol,xmin,fmin
      integer np, ndim
      real ps(np), pi(np), x0(ndim)
      external rpos

      ! local variables
      integer itmax
      parameter (itmax=100)
      real cgold,zeps
      parameter (cgold=.3819660,zeps=1.0d-14)
      integer iter
      real a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm

      ! functions
      real f1dim
!-----------------------------------------------------------------------
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f1dim(x,ps,pi,np,x0,ndim,rpos)
      fv=fx
      fw=fx
      do iter=1,itmax
         xm=0.5*(a+b)
         tol1=tol*abs(x)+zeps
         tol2=2.*tol1
         if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
         if(abs(e).gt.tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.*(q-r)
            if(q.gt.0.) p=-p
            q=abs(q)
            etemp=e
            e=d
            if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.
     $           p.ge.q*(b-x)) goto 1
            d=p/q
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2
         endif
 1       if(x.ge.xm) then
            e=a-x
         else
            e=b-x
         endif
         d=cgold*e
 2       if(abs(d).ge.tol1) then
            u=x+d
         else
            u=x+sign(tol1,d)
         endif
         fu=f1dim(u,ps,pi,np,x0,ndim,rpos)
         if(fu.le.fx) then
            if(u.ge.x) then
               a=x
            else
               b=x
            endif
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
         else
            if(u.lt.x) then
               a=u
            else
               b=u
            endif
            if(fu.le.fw .or. w.eq.x) then
               v=w
               fv=fw
               w=u
               fw=fu
            else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
               v=u
               fv=fu
            endif
         endif      
      enddo
      write(*,*) 'brent exceed maximum iterations'
 3    xmin=x
      fmin=fx
      iterd = iter
      return
      end subroutine
!=======================================================================
!> @brief 1D function min using derivative; Numerical Recipes 10.3
!! @param[in]   ax, bx, cx    minimum brackets
!! @param[in]   tol           tolerance
!! @param[out]  xmin          minimum position
!! @param[out]  fmin          function min
!! @param[out]  iterd         iteration count
!! @param[in]   ps            starting point
!! @param[in]   pi            displacement vector
!! @param[in]   np            number of parameters (manifold dimension)
!! @param[in]   x0            reference point
!! @param[in]   ndim          space dimension
!! @param[in]   rpos          subroutine returning surface pos. for given parameter
!! @param[in]   rdpos         subroutine returning derivatives for given parameter
      subroutine dbrent(ax,bx,cx,tol,xmin,fmin,iterd,ps,pi,np,x0,ndim,
     $     rpos,rdpos)
      implicit none

      ! argument list
      integer iterd
      real ax,bx,cx,tol,xmin,fmin
      integer np, ndim
      real ps(np), pi(np), x0(ndim)
      external rpos, rdpos

      ! local variables
      integer itmax
      parameter (itmax=100)
      real zeps
      parameter (zeps=1.0d-14)
      integer iter
      real a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,
     $     u,u1,u2,v,w,x,xm
      real dsign
      logical ok1,ok2   

      ! functions
      real f1dim, df1dim
!-----------------------------------------------------------------------
      a=min(ax,cx)
      b=max(ax,cx)
      if (ax.eq.a) then
         dsign = 1.0
      else
         dsign = -1.0
      endif
      v=bx
      w=v
      x=v
      e=0.
      fx=f1dim(x,ps,pi,np,x0,ndim,rpos)
      fv=fx
      fw=fx
      dx=dsign*df1dim(x,ps,pi,np,x0,ndim,rpos,rdpos)
      dv=dx
      dw=dx
      do iter=1,itmax
         xm=0.5*(a+b)
         tol1=tol*abs(x)+zeps
         tol2=2.*tol1
         if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
         if(abs(e).gt.tol1) then
            d1=2.*(b-a)
            d2=d1
            if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
            if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
            u1=x+d1
            u2=x+d2
            ok1=((a-u1)*(u1-b).gt.0.).and.(dx*d1.le.0.)
            ok2=((a-u2)*(u2-b).gt.0.).and.(dx*d2.le.0.)
            olde=e
            e=d
            if(.not.(ok1.or.ok2))then
               goto 1
            else if (ok1.and.ok2)then
               if(abs(d1).lt.abs(d2))then
                  d=d1
               else
                  d=d2
               endif
            else if (ok1)then
               d=d1
            else
               d=d2
            endif
            if(abs(d).gt.abs(0.5*olde))goto 1
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2
         endif
 1       if(dx.ge.0.) then
            e=a-x
         else
            e=b-x
         endif
         d=0.5*e
 2       if(abs(d).ge.tol1) then
            u=x+d
            fu=f1dim(u,ps,pi,np,x0,ndim,rpos)
         else
            u=x+sign(tol1,d)
            fu=f1dim(u,ps,pi,np,x0,ndim,rpos)
            if(fu.gt.fx)goto 3
         endif
         du=dsign*df1dim(u,ps,pi,np,x0,ndim,rpos,rdpos)
         if(fu.le.fx) then
            if(u.ge.x) then
               a=x
            else
               b=x
            endif
            v=w
            fv=fw
            dv=dw
            w=x
            fw=fx
            dw=dx
            x=u
            fx=fu
            dx=du
         else
            if(u.lt.x) then
               a=u
            else
               b=u
            endif
            if(fu.le.fw .or. w.eq.x) then
               v=w
               fv=fw
               dv=dw
               w=u
               fw=fu
               dw=du
            else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
               v=u
               fv=fu
               dv=du
            endif
         endif
      enddo
      write(*,*) 'ERROR; dbrent exceeded maximum iterations'
 3    xmin=x
      fmin=fx
      iterd=iter
      return
      end subroutine
!=======================================================================
!> @brief 1D function value wrapper; Numerical Recipes 10.5
!! @param[in] xp   position along the direction
!! @param[in] ps    starting point
!! @param[in] pi    displacement vector
!! @param[in] np    number of parameters (manifold dimension)
!! @param[in] x0    reference point
!! @param[in] ndim  space dimension
!! @param[in] rpos subroutine returning surface pos. for given parameter
!! @return function value
      function f1dim(xp,ps,pi,np,x0,ndim,rpos)
      implicit none

      ! argument list
      real xp
      integer np, ndim
      real ps(np), pi(np), x0(ndim)
      external rpos

      ! returns
      real f1dim

      ! local variables
      integer jl
      real xt(np)

      ! functions
      real sfunc
!-----------------------------------------------------------------------
      do jl=1,np
         xt(jl)=ps(jl)+xp*pi(jl)
      enddo
      f1dim=sfunc(xt,np,x0,ndim,rpos)
      return
      end function
!=======================================================================
!> @brief 1D function derivative wrapper; Numerical Recipes 10.6
!! @param[in] xp    position along the direction
!! @param[in] ps    starting point
!! @param[in] pi    displacement vector
!! @param[in] np    number of parameters (manifold dimension)
!! @param[in] x0    reference point
!! @param[in] ndim  space dimension
!! @param[in] rpos  subroutine returning surface pos. for given parameter
!! @param[in] rdpos subroutine returning derivatives for given parameter
!! @return inner product of direction and derivative     
      function df1dim(xp,ps,pi,np,x0,ndim,rpos,rdpos)
      implicit none

      ! argument list
      real xp
      integer np, ndim
      real ps(np), pi(np), x0(ndim)
      external rpos,rdpos

      ! returns
      real df1dim

      ! local variables
      integer jl
      real df(np),xt(np)
!-----------------------------------------------------------------------
      do jl=1,np
         xt(jl)=ps(jl)+xp*pi(jl)
      enddo
      call dsfunc(xt,df,np,x0,ndim,rpos,rdpos)
      df1dim=0.
      do jl=1,np
         df1dim=df1dim+df(jl)*pi(jl)
      enddo
      return
      end function
!=======================================================================
!> @brief Distance function
!! @param[in]    xp    parameter set
!! @param[in]    np    number of parameters (manifold dimension)
!! @param[in]    x0    reference point
!! @param[in]    ndim  space dimension
!! @param[in]    rpos  subroutine returning surface pos. for given parameter
      function sfunc(xp,np,x0,ndim,rpos)
      implicit none

      ! argument list
      integer np, ndim
      real xp(np)
      real x0(ndim)
      external rpos

      ! returns
      real sfunc

      ! local variables
      integer il
      real xs(ndim)
!-----------------------------------------------------------------------   
      call rpos(xs,ndim,xp,np)
      sfunc = 0.0
      do il=1,ndim
         sfunc = sfunc + (x0(il) - xs(il))**2 
      enddo
      end function
!=======================================================================
!> @brief Partial derivatives of the distance function
!! @param[in]    xp    parameter set
!! @param[out]   df    partial derivatives
!! @param[in]    np    number of parameters (manifold dimension)
!! @param[in]    x0    reference point
!! @param[in]    ndim  space dimension
!! @param[in]    rpos  subroutine returning surface pos. for given parameter
!! @param[in]    rdpos subroutine returning derivatives for given parameter
      subroutine dsfunc(xp,df,np,x0,ndim,rpos,rdpos)
      implicit none

      ! argument list
      integer np, ndim
      real xp(np), df(np)
      external rpos,rdpos
      real x0(ndim)

      ! local variables
      integer il, jl
      real xs(ndim), dxs(ndim,np)
!-----------------------------------------------------------------------
      call rpos(xs,ndim,xp,np)
      call rdpos(dxs,ndim,xp,np)
      do il=1,ndim
         xs(il) = 2.0*(x0(il) - xs(il))
      enddo
      do il=1,np
         df(il) = 0.0
         do jl=1,ndim
            df(il) = df(il) + xs(jl)*dxs(jl,il)
         enddo
      enddo
      end subroutine
!=======================================================================
!> @brief Get parameter boundaries for given direction
!! @param[inout] ps    starting point
!! @param[inout] xi    n-dimensional direction
!! @param[in]    np    number of variables
!! @param[out]   x_min lower scaling limit
!! @param[out]   x_max upper scaling limit
!! @param[in]    rbnd  subroutine returning parameter boundaries
      subroutine bnd_par(ps,xi,np,x_min,x_max,rbnd)
      implicit none

      ! argument list
      integer np
      real ps(np),xi(np),x_min,x_max
      external rbnd

      ! local variables
      integer il,jl
      real bnd(2,np), btmp(2), rtmp
!-----------------------------------------------------------------------
      x_min=-1.0d50
      x_max= 1.0d50
      ! get parameter boundaries for given surface
      call rbnd(bnd,np)
      do il=1,np
         if (xi(il).ne.0.0) then
            do jl=1,2
               btmp(jl) =(bnd(jl,il)-ps(il))/xi(il)
            enddo
            if (btmp(1).gt.btmp(2)) then
               rtmp = btmp(1)
               btmp(1) = btmp(2)
               btmp(2) = rtmp
            endif
            x_min=max(x_min,btmp(1))
            x_max=min(x_max,btmp(2))
         endif
      enddo
      return
      end subroutine
!=======================================================================
