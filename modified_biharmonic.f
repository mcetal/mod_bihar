      program MODBIHARMONIC
c
c Multiply connected with new kernels
c use Bremer's trick
c
      implicit real*8 (a-h,o-z)
      parameter (ndmax = 512, kmax = 4, nmax = kmax*ndmax)
c
c Geometry of domain
      dimension ak(kmax), bk(kmax), phi_k(kmax), ncyc(kmax)
      dimension rkappa(nmax), dsdth(nmax), ds_k(kmax)
      complex*16 z(nmax), dz(nmax), zn(nmax), zs(nmax), zk(kmax)
c
c Boundary conditions
      dimension bc_n(nmax), bc_s(nmax)
c
c Density
      dimension sigma1(nmax), sigma2(nmax), alpha1(nmax), alpha2(nmax),
     1          Palpha2(nmax)
      dimension Ealpha1(nmax), Ealpha2(nmax)
c
c  Matrix equation variables for GMRES
      parameter (maxl = 200)
      parameter (lrwork=10+2*nmax*(maxl+6)+maxl*(maxl+3), liwork=30)
      dimension gmwork(lrwork), iwork(2*nmax)
      dimension rhs(2*nmax), soln(2*nmax)
c
c Special function arrays
      dimension psi(100), factorial(0:100)
c
c  Quadrature weights and nodes
      dimension uweights(16), vnodes(16), ds_vec1(nmax*15),
     *          ds_vec2(nmax*15)
      complex*16 z_vec1(nmax*15), z_vec2(nmax*15), zn_vec1(nmax*15),
     *           zn_vec2(nmax*15), zs_vec1(nmax*15), zs_vec2(nmax*15)
c
c bessel functions
      dimension besselK0(nmax*nmax), besselK1(nmax*nmax)
c
c log singularities in multiply connected case
      dimension aK0_k(kmax)
c
c series expansions for terms in kernel
      dimension c0(nmax*nmax), c1(nmax*nmax), c2(nmax*nmax), 
     1          c3(nmax*nmax), c4(nmax*nmax), c5(nmax*nmax)
c
c Solution on grid
      parameter (nrmax = 300, nthmax = 300, ngmax = nrmax*nthmax)
      dimension sol_grid(ngmax), x_grid(ngmax), y_grid(ngmax), 
     1          igrid(ngmax)
      complex*16 z_tar(100), zgrad_pnt(100)
c
c Other arrays
      dimension w1(2*nmax), w2(2*nmax), amat(4*nmax**2)
      complex*16 zf(nmax), wsave(20*nmax), zf2(nmax)
      REAL*4 TIMEP(2), ETIME
c
c Common Blocks
      common /geometry_info/ z, dz, zn, zs, zk, rkappa, dsdth, ds_k
      common /ellipse_info/ ak, bk, phi_k, ncyc
      common /size_info/ k0, k, nd, nbk
      common /other_parameters/ alambda
      common /alpert/ uweights, vnodes, numpts, ishift, iorder, iquad
      common /geo_nodes/ ds_vec1, ds_vec2, z_vec1, z_vec2
      common /geo_nodes2/ zn_vec1, zn_vec2, zs_vec1, zs_vec2
      common /bessel_fns/ besselK0, besselK1
      common /coefficients/ c0, c1, c2, c3, c4, c5
      common /cutoff/ rcut
c
c Initialize and read in geometry data
         call PRINI (6,13)
         rcut = 0.02d0
         call READ_DATA (k0, k, nd, nbk, ak, bk, phi_k, ncyc, zk,
     1                   alambda,iquad)
         call DCFFTI (nd, wsave)
c
c Get quadrature weights and nodes for Alpert's rules
         call QUAD2 (vnodes, uweights, numpts, ishift, iorder, iquad)
         call prinf (' ORDER OF QUADRATURE = *', iorder, 1)
c
c Construct hole geometry 
         call MAKE_GEO (k0, k, nd, nbk, ak, bk, phi_k, ncyc, z, zk,  
     1                  dz, zn, zs, dsdth, rkappa, ds_k)
c
c Construct grid geometry
         nx = 200
         ny = 200
ccc         call GET_GRID (k0, k, nd, nx, ny, ak, bk, phi_k, ncyc, zk, 
ccc     1                  x_grid, y_grid, igrid)
         call GET_GRID2 (k0, k, nd, nx, ny, ak, bk, phi_k, ncyc, zk, 
     1                     x_grid, y_grid, igrid)
ccc         call GET_POLAR_GRID (k0, k, nd, nx, ny, ak, bk, phi_k, ncyc,  
ccc     1                  zk, x_grid, y_grid, igrid)
c
c Get geometry data at nodes
         call REAL_AT_NODES (k0, k, nd, nbk, vnodes, numpts, dsdth,  
     1                      ds_vec1, ds_vec2, zf, zf2, wsave)
         call COMPLEX_AT_NODES (k0, k, nd, nbk, vnodes, numpts, z,   
     1                          z_vec1, z_vec2, zf, zf2, wsave)
         call COMPLEX_AT_NODES (k0, k, nd, nbk, vnodes, numpts, zn,   
     1                          zn_vec1, zn_vec2, zf, zf2, wsave)
         call COMPLEX_AT_NODES (k0, k, nd, nbk, vnodes, numpts, zs,   
     1                          zs_vec1, zs_vec2, zf, zf2, wsave)
c
c Look up tables of Bessel Functions
         call GET_BESSEL (nbk, alambda, z, besselK0, besselK1,  
     1                    c0, c1, c2, c3, c4, c5)
c
c Construct boundary conditions
         call GET_BCS (k0, k, nd, nbk, alambda, zk, z, zn, zs, 
     1                 bc_n, bc_s, rhs, dsdth, zf, wsave, ds_k, ak)
ccc         call GET_BCS_REORD (k0, k, nd, nbk, alambda, zk, z, zn, zs, 
ccc     1                       bc_n, bc_s, rhs, dsdth, zf, wsave)
ccc         call prin2 (' bc_n exact = *', rhs, nbk)
ccc         call prin2 (' bc_s exact = *', rhs(nbk+1), nbk)
ccc         call OPERATE_E (k0, k, nd, nbk, rhs, rhs(nbk+1), dsdth,  
ccc     1                      rkappa, Ealpha1, Ealpha2, Palpha2, zf, 
ccc     2                      wsave)
ccc         call prin2 (' ealpha1 = *', Ealpha1, nbk)
ccc         call prin2 (' ealpha2 = *', Ealpha2, nbk)
ccc         call DIAG (k0, k, nd, nbk, Ealpha1, Ealpha2, dsdth,  
ccc     1                      rkappa, alpha1, alpha2, Palpha2, zf, 
ccc     2                      wsave)
ccc         call prin2 (' alpha1 = *', alpha1, nbk)
ccc         call prin2 (' alpha2 = *', alpha2, nbk)
ccc         stop
cccc
ccc         call MATVEC (2*nbk,rhs,w1,nelt,ia,ja,dummy_mat,isym)
ccc         open (unit = 22, file = 'Dbcs.dat')
ccc         if = 22
ccc         do i = 1, 2*nbk
ccc            write(if,'(e20.13,$)')(w1(i))
ccc            write (if,'(a)')  ''
ccc         end do
ccc         close(22)
ccc         stop
ccc         call PRIn2 (' new alpha1 = *', w1, nbk)
ccc         call prin2 (' new alpha2 = *', w1(nbk+1), nbk)
ccc         stop
c
c Solve
         call SOLVE (maxl, k0, k, nd, nbk, rhs, soln, alpha1, alpha2, 
     1               Palpha2, sigma1, sigma2, dsdth, rkappa, zf,  
     1               wsave, gmwork, lrwork, iwork, liwork, aK0_k, ds_k)
ccc         call SOLVE_REORD (maxl, k0, k, nd, nbk, rhs, soln, alpha1,  
ccc     1               alpha2, Palpha2, sigma1, sigma2, dsdth, rkappa, zf,  
ccc     1               wsave, gmwork, lrwork, iwork, liwork, aK0_k)
ccc         call prin2 (' rhs after solve = *', rhs, 2*nbk)
ccc         do i = 1, nbk
ccc            w1(i) = alpha1(i)
ccc            w1(nbk+i) = alpha2(i)
ccc         end do
ccc         do kbod = 1, k
ccc            w1(2*nbk+kbod) = aK0_k(kbod)
ccc         end do
ccc         call MATVEC (2*nbk,w1,w2,nelt,ia,ja,dummy_mat,isym)
ccc         call prin2 (' w2 = *', w2, 2*nbk)
ccc         call prin2 (' bc_n = *', w2, nbk)
ccc         call prin2 (' bc_s = *', w2(nbk+1), nbk)
ccc         stop
c
c construct matrix and dump it out so it can be read in by matlab
ccc         t0 = etime(timep)
ccc         call GETMAT1 (k, nd, nbk, w1, w2, amat, gmwork, iwork)
ccc         t1 = etime(timep)
ccc         call PRIN2 (' time in GETMAT1 = *', t1-t0, 1)
ccc         call GETMAT2 (k, nd, nbk, w1, w2, zf, wsave)
c
ccc         t0 = etime(timep)
ccc         call BUILD_MAT (2*nbk, amat)
ccc         t1 = etime(timep)
ccc         call PRIN2 (' time in BUILD_MAT = *', t1-t0, 1)
c
c Construct solution on grid
         call GET_SOL_GRID (k0, k, nd, nbk, nx, ny, alambda, sigma1,  
     1                      sigma2, z, zn, dsdth, x_grid, y_grid, 
     2                      igrid, sol_grid, zf, zf2, wsave, aK0_k)
ccc         call GET_SOL_POLAR_GRID (k0, k, nd, nbk, nx, ny, alambda,   
ccc     1                      sigma1, sigma2, z, zn, dsdth, x_grid,  
ccc     2                      y_grid, igrid, sol_grid, zf, zf2, wsave, 
ccc     3                      aK0_k)
c
c Get gradient at test points
         call GET_GRAD_PNT (k0, k, nd, nbk, alambda, sigma1, sigma2, 
     1                      z, zn, dsdth, zf, zf2, wsave, z_tar, 
     2                      zgrad_pnt, zk, aK0_k)
c
c Check error in case of a known analytic solution
         call CHECK_ERROR (k0, k, zk, nx, ny, alambda, x_grid, y_grid, 
     1                     igrid, sol_grid, z_tar, zgrad_pnt, ak)
c
c Dump out solution on grid
         open (unit = 31, file = 'xgrid.dat')
         open (unit = 32, file = 'ygrid.dat')
         open (unit = 33, file = 'igrid.dat')
         open (unit = 34, file = 'ugrid.dat')
         call DUMP (nx, ny, x_grid, igrid, 1, 31)
         call DUMP (nx, ny, y_grid, igrid, 1, 32)
         call DUMP (nx, ny, y_grid, igrid, 0, 33)
         call DUMP (nx, ny, sol_grid, igrid, 1, 34)
c
      stop
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine READ_DATA (k0, k, nd, nbk, ak, bk, phi_k, ncyc, zk,
     1                      alambda, iquad)
c---------------
c
c  k0 = 0 if interior problem, 1 if exterior
c  k  = degree of connectivity (i.e. k is the number of holes)
c  nd = points per hole
c
      implicit real*8 (a-h,o-z)
      dimension ak(*), bk(*), phi_k(*), ncyc(*)
      complex*16 zk(*)
         open (unit = 11, file = 'input.dat')
c
         read (11,*) k0, k, nd
         read (11,*) alambda
         read (11,*) iquad
         nbk = (k-k0+1)*nd
         call PRINF (' k0 = *', k0, 1)
         call PRINF (' k = *', k, 1)
         call PRINF (' nbk = *', nbk, 1)
         call PRIN2 (' alambda = *', alambda, 1)
         call PRINF (' iquad = *', iquad, 1)
         do kbod = k0, k
            i = kbod-k0+1
            read(11,*) ak(i), bk(i), phi_k(i), ncyc(i), cx, cy
            zk(i) = dcmplx(cx,cy)
         end do
         call PRIN2 (' zk = *', zk, 2*(k-k0+1))
c
         close(11)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine ELLIPSE (alpha, a, b, phi_k, zk, x, y, dx, dy, 
     1                    d2x, d2y, NCYC) 
c---------------
      implicit real*8 (a-h,o-z)
      complex*16 zk
c
c returns parametrization of a rotated ellipse:
c  x = cx + a cos(phi) cos(alpha) - b*sin(phi)*sin(alpha), 
c  y = cy + b cos(phi) sin(alpha) + a*sin(phi)*cos(alpha)
c
         iflag = 1
         if (iflag.eq.0) then 
            x = dreal(zk) +
     1          a*dcos(phi)*dcos(alpha) - b*dsin(phi)*dsin(alpha)
               dx = -a*dcos(phi)*dsin(alpha) - b*dsin(phi)*dcos(alpha)
               d2x = -a*dcos(phi)*dcos(alpha) + b*dsin(phi)*dsin(alpha)
            y = dimag(zk) + 
     1          b*dcos(phi)*dsin(alpha) + a*dsin(phi)*dcos(alpha)
               dy = b*dcos(phi)*dcos(alpha) - a*dsin(phi)*dsin(alpha)
               d2y = -b*dcos(phi)*dsin(alpha) - a*dsin(phi)*dcos(alpha)
          else
ccc            NCYC = 6
            eps = b
            cs = dcos(alpha)
            sn = dsin(alpha)
            snn = dsin(NCYC*alpha)
            csn = dcos(NCYC*alpha)
            r = dsqrt(a**2 + b**2 + 2.d0*a*b*dcos((NCYC+1)*alpha))
            rdot = -a*b*(NCYC+1)*dsin((NCYC+1)*alpha)/r
            r2dot = -a*b*(NCYC+1)**2*dcos((NCYC+1)*alpha) - rdot**2
            r2dot = r2dot/r
            x = dreal(zk) + r*dcos(alpha)
            y = dimag(zk) + r*dsin(alpha)
               dx = rdot*dcos(alpha) - r*dsin(alpha)
               d2x = r2dot*dcos(alpha) -2.d0*rdot*dsin(alpha) 
     1                 - r*dcos(alpha)
               dy = rdot*dsin(alpha) + r*dcos(alpha)
               d2y = r2dot*dsin(alpha) + 2.d0*rdot*dcos(alpha) 
     1                 - r*dsin(alpha)
         end if

      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine MAKE_GEO (k0, k, nd, nbk, ak, bk, phi_k, ncyc, z, zk, 
     1                     dz, zn, zs, dsdth, rkappa, ds_k)
c---------------
      implicit real*8 (a-h,o-z)
      dimension ak(k0:k), bk(k0:k), phi_k(k0:k), ncyc(k0:k), ds_k(k0:k)
      dimension rkappa(nbk), dsdth(nbk)
      complex*16 z(nbk), dz(nbk), zn(nbk), zs(nbk), zk(k0:k)
c
         open (unit = 15, file = 'geometry.m')
c
         pi = 4.d0*datan(1.d0)
c
         dalph = 2.d0*pi/nd
         ist = 0
         do kbod = k0, k
            dsmin = 1.d10
            dsmax = 0.d0
            do i = 1, nd
               call PRINF('+++ i = *', i, 1)
               alpha = dalph*(i-1.d0)
               call ELLIPSE (alpha, ak(kbod), bk(kbod), phi_k(kbod), 
     1                       zk(kbod), x, y, dx, dy, d2x, d2y,
     2                       NCYC(kbod))
               call PRIN2 (' dx = *', dx, 1)
               call PRIN2 (' dy = *', dy, 1)
               call PRIN2 (' d2x = *', d2x, 1)
               call PRIN2 (' d2y = *', d2y, 1)
               z(ist+i) = dcmplx(x,y)
               dsdth(ist+i) = dsqrt(dx**2 + dy**2)
               dsmin = min(dsmin,dsdth(ist+i))
               dsmax = max(dsmax,dsdth(ist+i))
               if (kbod.eq.0) then 
                  dz(ist+i) = dcmplx(dx,dy)
                  zs(ist+i) = dz(ist+i)/dsdth(ist+i)
                  zn(ist+i) = dcmplx(dy,-dx)/dsdth(ist+i)
                  rkappa(ist+i) = (dx*d2y-dy*d2x)/dsdth(ist+i)**3
                 else
                  dz(ist+i) = -dcmplx(dx,dy)
                  zs(ist+i) = dz(ist+i)/dsdth(ist+i)
                  zn(ist+i) = -dcmplx(dy,-dx)/dsdth(ist+i)
                  rkappa(ist+i) = -(dx*d2y-dy*d2x)/dsdth(ist+i)**3
               end if
            end do
            dsmin = dsmin*2.d0*pi/nd
            dsmax = dsmax*2.d0*pi/nd
ccc            ds_k(kbod) = 1.d0
            ds_k(kbod) = dsmin
            call RSCPLOT (z(ist+1),nd,1,15)
            ist = ist + nd
         end do
         call PRIN2 (' ds_k = *', ds_k(k0), k-k0+1)
         call PRIN2 (' dsdth = *', dsdth, nbk)
         call PRIN2 (' rkappa = *', rkappa, nbk)
         close(15)
ccc         stop
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GET_GRID (k0, k, nd, nx, ny, ak, bk, phi_k, ncyc, zk, 
     1                     x_grid, y_grid, igrid)
c---------------
c Lays down a grid; assumes a single hole
c
      implicit real*8 (a-h,o-z)
      dimension ak(k0:k), bk(k0:k), phi_k(k0:k), ncyc(k0:k)
      dimension x_grid(nx,ny), y_grid(nx,ny), igrid(nx,ny)
      complex*16 zk(k0:k)
c
         pi = 4.d0*datan(1.d0)
         dth = 2.d0*pi/nd
         nfac = 2
c
c find dimensions of box         
         if (k0.eq.0) then
            xmax = dreal(zk(0)) + ak(0)
            xmin = dreal(zk(0)) - ak(0)
            ymax = dimag(zk(0)) + bk(0)
            ymin = dimag(zk(0)) - bk(0)
           else
            xmax = 4.d0
            xmin = -4.d0
            ymax = 4.d0
            ymin = -4.d0
         end if
         call PRIN2 (' xmin = *', xmin, 1)
         call PRIN2 (' xmax = *', xmax, 1)
         call PRIN2 (' ymin = *', ymin, 1)
         call PRIN2 (' ymax = *', ymax, 1)
c
         dx = (xmax - xmin)/(nx - 1.d0)
         dy = (ymax - ymin)/(ny - 1.d0)
c
         do ix = 1, nx
            x = xmin + (ix-1.d0)*dx
            do iy = 1, ny
               y = ymin + (iy-1.d0)*dy
ccc               call prin2 (' x = *', x, 1)
ccc               call prin2 (' y = *', y, 1)
               x_grid(ix,iy) = x
               y_grid(ix,iy) = y
               isafe = 0
               do kbod = k0, k
                  dsmax = max(ak(kbod),bk(kbod))*dth
                  eps = nfac*dsmax
                  xk = dreal(zk(kbod))
                  yk = dimag(zk(kbod))
                  if (kbod.eq.0) then
                     asafe = ak(kbod) - eps
                     bsafe = bk(kbod) - eps
                    else
                     asafe = ak(kbod) + eps
                     bsafe = bk(kbod) + eps
                  end if
                  check = ((x-xk)/asafe)**2 + ((y-yk)/bsafe)**2
                  if ((kbod.ge.1).and.(check.lt.1.d0)) then
                     isafe = isafe + 1
                    elseif((kbod.eq.0).and.(check.gt.1.d0)) then 
                     isafe = isafe + 1
                  end if
               end do
               if (isafe.ne.0) then
                  igrid(ix,iy) = 0
                 else
                  igrid(ix,iy) = 1
               end if
            end do
         end do
c
         close(21)
         close(22)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GET_POLAR_GRID (k0, k, nd, nx, ny, ak, bk, phi_k, ncyc,  
     1                     zk, x_grid, y_grid, igrid)
c---------------
c Lays down a grid; x = r y = th
c for concentric circle example only!
c
      implicit real*8 (a-h,o-z)
      dimension ak(k0:k), bk(k0:k), phi_k(k0:k), ncyc(k0:k)
      dimension x_grid(nx,ny), y_grid(nx,ny), igrid(nx,ny)
      complex*16 zk(k0:k)
c
         pi = 4.d0*datan(1.d0)
         dth = 2.d0*pi/nd
         nfac = 5
         eps = nfac*ak(1)*dth
c
c find dimensions of box  
         r0 = ak(1)+nfac*ak(1)*dth
         r1 = ak(0)-nfac*ak(1)*dth
         call PRIN2 ('MIN R = *', r0, 1)
         call PRIN2 ('MAX R = *', r1, 1)
c
         dx = (r1 - r0)/(nx - 1.d0)
         dy = 2.d0*pi/(ny - 1.d0)
c
         do ix = 1, nx
            x = r0 + (ix-1.d0)*dx
            do iy = 1, ny
               y = ymin + (iy-1.d0)*dy
ccc               call prin2 (' x = *', x, 1)
ccc               call prin2 (' y = *', y, 1)
               x_grid(ix,iy) = x
               y_grid(ix,iy) = y
               igrid(ix,iy) = 1
            end do
         end do
c
         close(21)
         close(22)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GET_GRID2 (k0, k, nd, nx, ny, ak, bk, phi_k, ncyc, zk, 
     1                     x_grid, y_grid, igrid)
c---------------
c Lays down a grid; assumes a single hole
c
      implicit real*8 (a-h,o-z)
      dimension ak(k0:k), bk(k0:k), phi_k(k0:k), ncyc(k0:k)
      dimension x_grid(nx,ny), y_grid(nx,ny), igrid(nx,ny)
      complex*16 zk(k0:k)
c
         pi = 4.d0*datan(1.d0)
         dth = 2.d0*pi/nd
         nfac = 8
ccc         N = 6
c
c find dimensions of box         
         if (k0.eq.0) then
            xmax = dreal(zk(0)) + ak(0) + bk(0)
            xmin = dreal(zk(0)) - ak(0) - bk(0)
            ymax = dimag(zk(0)) + bk(0) + ak(0)
            ymin = dimag(zk(0)) - bk(0) - ak(0)
           else
            xmax = 3.5d0
            xmin = -3.5d0
            ymax = 3.5d0
            ymin = -3.5d0
         end if
         call PRIN2 (' xmin = *', xmin, 1)
         call PRIN2 (' xmax = *', xmax, 1)
         call PRIN2 (' ymin = *', ymin, 1)
         call PRIN2 (' ymax = *', ymax, 1)
c
         dx = (xmax - xmin)/(nx - 1.d0)
         dy = (ymax - ymin)/(ny - 1.d0)
c
         do ix = 1, nx
            x = xmin + (ix-1.d0)*dx
            do iy = 1, ny
               y = ymin + (iy-1.d0)*dy
ccc               call prin2 (' x = *', x, 1)
ccc               call prin2 (' y = *', y, 1)
               x_grid(ix,iy) = x
               y_grid(ix,iy) = y
               isafe = 0
               do kbod = k0, k
                  N = NCYC(kbod)
                  dsmax = ak(kbod)*dth
                  eps = nfac*dsmax
                  xk = dreal(zk(kbod))
                  yk = dimag(zk(kbod))
                  r = dsqrt((x-xk)**2 + (y-yk)**2)
                  if (yk.ge.0.d0) then
                     alpha = dacos((x-xk)/r)
                    else
                     alpha = 2.d0*pi - dacos((x-xk)/r) 
                  end if
                  if (k0.eq.0) then
                     rc = dsqrt(ak(kbod)**2 + bk(kbod)**2 
     1                + 2.d0*ak(kbod)*bk(kbod)*dcos((N+1)*alpha)) - eps
                    else
                     rc = dsqrt(ak(kbod)**2 + bk(kbod)**2 
     1                + 2.d0*ak(kbod)*bk(kbod)*dcos((N+1)*alpha)) + eps
                  end if
                  check = r/rc
                  if ((kbod.ge.1).and.(check.lt.1.d0)) then
                     isafe = isafe + 1
                    elseif((kbod.eq.0).and.(check.gt.1.d0)) then 
                     isafe = isafe + 1
                  end if
               end do
               if (isafe.ne.0) then
                  igrid(ix,iy) = 0
                 else
                  igrid(ix,iy) = 1
               end if
            end do
         end do
c
         close(21)
         close(22)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine COMPLEX_AT_NODES (k0, k, nd, nbk, v, numpts, z,   
     1                             z_vec1, z_vec2, zf, zf2, wsave)
c---------------
      implicit real*8 (a-h,o-z)
      dimension v(16)
      complex*16 z(nbk), z_vec1(nbk,15), z_vec2(nbk,15), eye
      complex*16 zf(nbk), zf2(nbk), wsave(*)
c
      pi = 4.d0*datan(1.d0)
      eye = dcmplx(0.d0,1.d0)
      dalph = 2.d0*pi/nd
      n = nd/2
      dth = 2.d0*pi/nd
c
         ist = 0
         do kbod = k0, k
c
c Fourier interpolate z to nodes
            do i = 1, nd
               zf(i) = z(ist+i)
            end do
            call DCFFTF (nd, zf, wsave)
            do i = 1, nd
               zf(i) = zf(i)/nd
            end do
            do inode = 1, numpts
               zf2(1) = zf(1)
               do kmode = 1, n-1
                  zf2(kmode+1) 
     1                 = zf(kmode+1)*cdexp(eye*kmode*v(inode)*dth)
                  zf2(nd-kmode+1) 
     1                 = zf(nd-kmode+1)*cdexp(-eye*kmode*v(inode)*dth)
               end do
               zf2(n+1) = 0.d0
               call DCFFTB (nd, zf2, wsave)
               do i = 1, nd
                  z_vec1(ist+i,inode) = zf2(i)
               end do
c
               zf2(1) = zf(1)
               do kmode = 1, n-1
                  zf2(kmode+1) 
     1                  = zf(kmode+1)*cdexp(-eye*kmode*v(inode)*dth)
                  zf2(nd-kmode+1) 
     1                  = zf(nd-kmode+1)*cdexp(eye*kmode*v(inode)*dth)
               end do
               zf2(n+1) = 0.d0
               call DCFFTB (nd, zf2, wsave)
               do i = 1, nd
                  z_vec2(ist+i,inode) = zf2(i)
               end do
            end do
            ist = ist+nd
         end do
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine MATRIX_MULT (nbk, inode, x, f_vec, fx, E_vec)
c---------------
      implicit real*8 (a-h,o-z)
      dimension f_vec(nbk,nbk,15), x(nbk), fx(nbk), E_vec(nbk,15)
      dimension uweights(16), vnodes(16)
      common /alpert/ uweights, vnodes, numpts, ishift, iorder, iquad
c
         do i = 1, nbk
            fx(i) = 0.d0
            do j = 1, nbk
               fx(i) = fx(i) + f_vec(i,j,inode)*x(j)
            end do
            call PRIN2 (' fx = *', fx(i), 1)
            call PRIn2 ('    E_vec = *', E_vec(i,inode), 1)
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine NODE_MAT (nd, numpts, v,    
     1                     f_vec1, f_vec2, zf, zf2, wsave)
c---------------
      implicit real*8 (a-h,o-z)
      dimension v(16), f_vec1(nd,nd,15), f_vec2(nd,nd,15)
      complex*16 eye
      complex*16 zf(*), zf2(*), wsave(*)
c
      pi = 4.d0*datan(1.d0)
      eye = dcmplx(0.d0,1.d0)
      dalph = 2.d0*pi/nd
      n = nd/2
      dth = 2.d0*pi/nd
c
         call DCFFTI (nd, wsave)
c
c Fourier interpolate z to nodes
            do i = 1, nd
               do j = 1, i-1
                  zf(j) = 0.d0
               end do
               zf(i) = 1.d0
               do j = i+1,nd
                  zf(j) = 0.d0
               end do
               call DCFFTF (nd, zf, wsave)
               do j = 1, nd
                  zf(j) = zf(j)/nd
               end do
               do inode = 1, numpts
                  zf2(1) = zf(1)
                  do kmode = 1, n-1
                     zf2(kmode+1) 
     1                 = zf(kmode+1)*cdexp(eye*kmode*v(inode)*dth)
                     zf2(nd-kmode+1) 
     1                 = zf(nd-kmode+1)*cdexp(-eye*kmode*v(inode)*dth)
                  end do
                  zf2(n+1) = 0.d0
                  call DCFFTB (nd, zf2, wsave)
                  do j = 1, nd
                     f_vec1(j,i,inode) = dreal(zf2(j))
                  end do
c
                  zf2(1) = zf(1)
                  do kmode = 1, n-1
                     zf2(kmode+1) 
     1                  = zf(kmode+1)*cdexp(-eye*kmode*v(inode)*dth)
                     zf2(nd-kmode+1) 
     1                  = zf(nd-kmode+1)*cdexp(eye*kmode*v(inode)*dth)
                  end do
                  zf2(n+1) = 0.d0
                  call DCFFTB (nd, zf2, wsave)
                  do j = 1, nd
                     f_vec2(j,i,inode) = dreal(zf2(j))
                  end do
               end do
            end do
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine REAL_AT_NODES (k0, k, nd, nbk, v, numpts, f,   
     1                          f_vec1, f_vec2, zf, zf2, wsave)
c---------------
      implicit real*8 (a-h,o-z)
      dimension v(16), f(nbk), f_vec1(nbk,15), f_vec2(nbk,15)
      complex*16 eye
      complex*16 zf(nbk), zf2(nbk), wsave(*)
c
      pi = 4.d0*datan(1.d0)
      eye = dcmplx(0.d0,1.d0)
      dalph = 2.d0*pi/nd
      n = nd/2
      dth = 2.d0*pi/nd
c
         ist = 0
         call DCFFTI (nd, wsave)
         do kbod = k0, k
c
c Fourier interpolate z to nodes
            do i = 1, nd
               zf(i) = f(ist+i)
            end do
            call DCFFTF (nd, zf, wsave)
            do i = 1, nd
               zf(i) = zf(i)/nd
            end do
            do inode = 1, numpts
               zf2(1) = zf(1)
               do kmode = 1, n-1
                  zf2(kmode+1) 
     1                 = zf(kmode+1)*cdexp(eye*kmode*v(inode)*dth)
                  zf2(nd-kmode+1) 
     1                 = zf(nd-kmode+1)*cdexp(-eye*kmode*v(inode)*dth)
               end do
               zf2(n+1) = 0.d0
               call DCFFTB (nd, zf2, wsave)
               do i = 1, nd
                  f_vec1(ist+i,inode) = dreal(zf2(i))
               end do
c
               zf2(1) = zf(1)
               do kmode = 1, n-1
                  zf2(kmode+1) 
     1                  = zf(kmode+1)*cdexp(-eye*kmode*v(inode)*dth)
                  zf2(nd-kmode+1) 
     1                  = zf(nd-kmode+1)*cdexp(eye*kmode*v(inode)*dth)
               end do
               zf2(n+1) = 0.d0
               call DCFFTB (nd, zf2, wsave)
               do i = 1, nd
                  f_vec2(ist+i,inode) = dreal(zf2(i))
               end do
            end do
            ist = ist+nd
         end do
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GET_BESSEL (nbk, alambda, z, besselK0, besselK1, 
     1                       c0, c1, c2, c3, c4, c5)
c---------------
      implicit real*8 (a-h,o-z)
      complex*16 z(nbk)
      dimension besselK0(nbk,nbk), besselK1(nbk,nbk), c0(nbk,nbk),
     1          c1(nbk,nbk), c2(nbk,nbk), c3(nbk,nbk), c4(nbk,nbk),
     2          c5(nbk,nbk)
      dimension psi(100), factorial(0:100), coeff(3)
      dimension bk(0:1)
      common /spec_fun/ psi, factorial, nterms
      common /cutoff/ rcut
c
      pi = 4.d0*datan(1.d0)
      dalph = 2.d0*pi/nd
c
         az_min = 1.d10
         do i = 1, nbk
            do j = 1, i-1
               r = cdabs(z(i) - z(j))
               az = alambda*r
               az_min = min(az_min,az)
               if (az.lt.100.d0) then
                  call RKBESL (alambda*r, 0.d0, 2, 1, bk, ncalc)
                 else
                  call RKBESL (alambda*r, 0.d0, 2, 2, bk, ncalc)
                  do ik = 0, 1
                     bk(ik) = dexp(-az)*bk(ik)
                  end do
               end if
               if (ncalc.ne.2) then 
                  write (6,*) 'ERROR IN GET_BESSEL! ', ncalc
                  stop
               end if
               besselK0(i,j) = bk(0)
               besselK1(i,j) = bk(1)
               c0(i,j) = bk(0) + 2.d0*bk(1)/az - 2.d0/az**2
               c1(i,j) = 3.d0*c0(i,j) + az*bk(1) + 0.5d0
               c2(i,j) = 4.d0*c0(i,j) + az*bk(1)
               c3(i,j) = 12.d0*c0(i,j) + 5.d0*az*bk(1) + az**2*bk(0)
     1                   + 1.d0
               c4(i,j) = 24.d0*c0(i,j) + 8.d0*az*bk(1) + az**2*bk(0)
               c5(i,j) = 2.d0*c0(i,j) + az*bk(1)
            end do
            besselK0(i,i) = 0.d0
            besselK1(i,i) = 0.d0
            c0(i,i) = 0.d0
            c1(i,i) = 0.d0
            c2(i,i) = 0.d0
            c3(i,i) = 0.d0
            c4(i,i) = 0.d0
            c5(i,i) = 0.d0
            do j = i+1,nbk
               r = cdabs(z(i) - z(j))
               az = alambda*r
               if (az.lt.100.d0) then
                  call RKBESL (alambda*r, 0.d0, 2, 1, bk, ncalc)
                 else
                  call RKBESL (alambda*r, 0.d0, 2, 2, bk, ncalc)
                  do ik = 0, 1
                     bk(ik) = dexp(-az)*bk(ik)
                  end do
               end if
               if (ncalc.ne.2) then 
                  write (6,*) 'ERROR IN GET_BESSEL! ', ncalc
                  stop
               end if
               besselK0(i,j) = bk(0)
               besselK1(i,j) = bk(1)
               c0(i,j) = bk(0) + 2.d0*bk(1)/az - 2.d0/az**2
               c1(i,j) = 3.d0*c0(i,j) + az*bk(1) + 0.5d0
               c2(i,j) = 4.d0*c0(i,j) + az*bk(1)
               c3(i,j) = 12.d0*c0(i,j) + 5.d0*az*bk(1) + az**2*bk(0)
     1                   + 1.d0
               c4(i,j) = 24.d0*c0(i,j) + 8.d0*az*bk(1) + az**2*bk(0)
               c5(i,j) = 2.d0*c0(i,j) + az*bk(1)
            end do
         end do
         call PRIN2 (' az min in GET BESSEL = *', az_min, 1)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GET_BESSEL_PNT (az, c0, c1, c2, c3, c4, c5)
c---------------
      implicit real*8 (a-h,o-z)
      dimension bk(0:1)
      data euler /0.57721566490153286d0/
      common /cutoff/ rcut
c
         pi = 4.d0*datan(1.d0)
c
               if (az.lt.100.d0) then
                  call RKBESL (az, 0.d0, 2, 1, bk, ncalc)
                 else
                  call RKBESL (az, 0.d0, 2, 2, bk, ncalc)
                  do ik = 0, 1
                     bk(ik) = dexp(-az)*bk(ik)
                  end do
               end if
         if (ncalc.ne.2) then 
                  write (6,*) 'ERROR IN GET_BESSEL_PNT! ', ncalc
                  stop
         end if
         if (az.gt.rcut) then
            c0 = bk(0) + 2.d0*bk(1)/az - 2.d0/az**2
            c1 = 3.d0*c0 + az*bk(1) + 0.5d0
            c2 = 4.d0*c0 + az*bk(1)
            c3 = 12.d0*c0 + 5.d0*az*bk(1) + az**2*bk(0)
     1                   + 1.d0
            c4 = 24.d0*c0 + 8.d0*az*bk(1) + az**2*bk(0)
            c5 = 2.d0*c0 + az*bk(1)
          else
c
c I've added additional terms from the paper
            alnz=dlog(2.d0/az)-euler
            az2=az**2
            az4=az2**2
            az6=az2**3
            az8=az4**2
            c0=-0.5d0 + (alnz+0.75d0)*az2/8.d0 
     1             + (alnz+17.d0/12.d0)*az4/96.d0
     2             + (alnz+43.d0/24.d0)*az6/3072.d0
     3             + (alnz+247.d0/120.d0)*az6/184320.d0
            c1=-(alnz-0.25d0)*az2/8.d0-(alnz+13.d0/12.d0)*az4/32.d0
     1         -5.d0*(alnz + 191.d0/120.d0)/3072.d0
     2         -7.d0*(alnz + 1609.d0/840.d0)/184320.d0
            c2=-1.d0+az2/8.d0-(alnz+11.d0/12.d0)*az4/48.d0
     1         -(alnz + 37.d0/24.d0)/768.d0
     2         -(alnz + 227.d0/120.d0)/30720.d0
            c3=-az2/8.d0+(alnz+7.d0/12.d0)*az4/16.d0
     1         +5.d0*(alnz + 161.d0/120.d0)/768.d0
     2         +7.d0*(alnz + 1469.d0/840.d0)/30720.d0
            c4=-4.d0+az2/4.d0-az4/48.d0
     1         +(alnz + 25.d0/24.d0)/384.d0
     2         +(alnz + 197.d0/120.d0)/7680.d0
            c5=-(alnz+0.25d0)*az2/4.d0-(alnz+7.d0/6.d0)*az4/24.d0
     1         -(alnz + 13.d0/8.d0)/512.d0
     2         -(alnz + 29.d0/15.d0)/23040.d0
         end if
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GET_BCS (k0, k, nd, nbk, alambda, zk, z, zn, zs, 
     1                    bc_n, bc_s, rhs, dsdth, zf, wsave, ds_k, ak)
c---------------
      implicit real*8 (a-h,o-z)
      complex*16 z(nbk), zn(nbk), zs(nbk), zk(k0:k), zgrad_psi
      complex*16 zf(nd), wsave(*)
      dimension bc_n(nbk), bc_s(nbk), rhs(2*nbk), dsdth(nbk)
      dimension bi(0:1), ds_k(k0:k), ak(k0:k)
      REAL*4 TIMEP(2), ETIME
c
      pi = 4.d0*datan(1.d0)
      dalph = 2.d0*pi/nd
      call DCFFTI (nd, wsave)
c
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               call EXACT_SOL (k0, k, z(ist+i), zk, alambda, psi, 
     1                         zgrad_psi, ak(k0))
               bc_n(ist+i) = dreal(zgrad_psi*dconjg(zn(ist+i)))
               bc_s(ist+i) = dreal(zgrad_psi*dconjg(zs(ist+i)))
               zf(i) = bc_s(ist+i)*dsdth(ist+i)
ccc               bc_n(ist+i) = kbod+1.d0
ccc               bc_s(ist+i) = 0.d0
               rhs(ist+i) = ds_k(kbod)*bc_n(ist+i)
               rhs(nbk+ist+i) = ds_k(kbod)*bc_s(ist+i)
            end do
            call DCFFTF (nd, zf, wsave)
            flux = 2.d0*pi*dreal(zf(1))/nd
            call PRINF (' body = *', kbod, 1)
            call PRIN2 ('   flux = *', flux, 1)
            ist = ist + nd
         end do
ccc         do kbod = 1, k
ccc            rhs(2*nbk+kbod) = 0.d0
ccc         end do
c
         open (unit = 21, file = 'bcs.dat')
         if = 21
         do i = 1, nbk
            write(if,'(e20.13,$)')(bc_n(i))
            write (if,'(a)')  ''
         end do
         do i = 1, nbk
            write(if,'(e20.13,$)')(bc_s(i))
            write (if,'(a)')  ''
         end do
         close(21)
c
         call PRIN2 (' bc_n = *', bc_n, nbk)
         call PRIN2 (' bc_s = *', bc_s, nbk)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GET_BCS_REORD (k0, k, nd, nbk, alambda, zk, z, zn, zs, 
     1                          bc_n, bc_s, rhs, dsdth, zf, wsave)
c---------------
      implicit real*8 (a-h,o-z)
      complex*16 z(nbk), zn(nbk), zs(nbk), zk(k0:k), zgrad_psi
      complex*16 zf(nd), wsave(*)
      dimension bc_n(nbk), bc_s(nbk), rhs(2*nbk), dsdth(nbk)
      dimension bi(0:1)
      REAL*4 TIMEP(2), ETIME
c
      pi = 4.d0*datan(1.d0)
      dalph = 2.d0*pi/nd
      call DCFFTI (nd, wsave)
c
         ist = 0
         ist2 = 0
         do kbod = k0, k
            do i = 1, nd
               call EXACT_SOL (k0, k, z(ist+i), zk, alambda, psi, 
     1                         zgrad_psi)
               bc_n(ist+i) = dreal(zgrad_psi*dconjg(zn(ist+i)))
               bc_s(ist+i) = dreal(zgrad_psi*dconjg(zs(ist+i)))
               zf(i) = bc_s(ist+i)*dsdth(ist+i)
ccc               bc_n(ist+i) = kbod+1.d0
ccc               bc_s(ist+i) = 0.d0
               rhs(ist2+i) = bc_n(ist+i)
               rhs(nd+ist2+i) = bc_s(ist+i)
            end do
            call DCFFTF (nd, zf, wsave)
            flux = 2.d0*pi*dreal(zf(1))/nd
            call PRINF (' body = *', kbod, 1)
            call PRIN2 ('   flux = *', flux, 1)
            ist = ist + nd
            ist2 = ist2 + 2*nd
         end do
ccc         do kbod = 1, k
ccc            rhs(2*nbk+kbod) = 0.d0
ccc         end do
c
         open (unit = 21, file = 'bcs.dat')
         if = 21
         do i = 1, nbk
            write(if,'(e20.13,$)')(bc_n(i))
            write (if,'(a)')  ''
         end do
         do i = 1, nbk
            write(if,'(e20.13,$)')(bc_s(i))
            write (if,'(a)')  ''
         end do
         close(21)
c
         call PRIN2 (' bc_n = *', bc_n, nbk)
         call PRIN2 (' bc_s = *', bc_s, nbk)
ccc         call PRIN2 (' rhs = *', rhs, 2*nbk)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine EXACT_SOL (k0, k, z, zk, alambda, psi, zgrad_psi, ak)
c---------------
c
      implicit real*8 (a-h,o-z)
      complex*16 zgrad_psi, zk(k0:k), z, z0, eye
      dimension bi(0:1), bk(0:1), ak(k0:k)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
         psi = 0.d0
         zgrad_psi = 0.d0
         itest = 0
c
         if (itest.eq.0) then 
c
c First example in paper
            z0 = -0.6d0+0.5d0*eye
            r = cdabs(z - z0) 
            az = alambda*r
                  if (az.lt.1000.d0) then
                     call RKBESL (alambda*r, 0.d0, 2, 1, bi, ncalc)
                    else
                     call RKBESL (alambda*r, 0.d0, 2, 2, bi, ncalc)
                     do ik = 0, 1
                        bi(ik) = dexp(-az)*bi(ik)
                     end do
                  end if
                  zgrad_psi = zgrad_psi
     1              - 0.5d0*(1.d0/r - alambda*bi(1))
     2                             *(z - z0)/(r)
                  psi = psi -0.5d0*(dlog(r)+bi(0))
           elseif (itest.eq.1) then
c
c Second example in paper
ccc            z0 = -0.625d0+0.625d0*eye
ccc            z0 = -2.5d0+2.5d0*eye
ccc            z0 = -.77d0+.7d0*eye
            z0 = -1.d0+1.d0*eye
            r = cdabs(z - z0) 
            az = alambda*r
                  if (az.lt.1000.d0) then
                     call RKBESL (alambda*r, 0.d0, 2, 1, bi, ncalc)
                    else
                     call RKBESL (alambda*r, 0.d0, 2, 2, bi, ncalc)
                     do ik = 0, 1
                        bi(ik) = dexp(-az)*bi(ik)
                     end do
                  end if
                  zgrad_psi = zgrad_psi
     1              - 0.5d0*(1.d0/r - alambda*bi(1))
     2                             *(z - z0)/(r)
                  psi = psi -0.5d0*(dlog(r)+bi(0))
            do kbod = 1, k
               r = cdabs(z - zk(kbod))
                  az = alambda*r
                  if (az.lt.1000.d0) then
                     call RKBESL (alambda*r, 0.d0, 2, 1, bi, ncalc)
ccc                     call RIBESL (alambda*r, 0.d0, 2, 1, bi, ncalc)
                    else
                     call RKBESL (alambda*r, 0.d0, 2, 2, bi, ncalc)
ccc                     call RIBESL (alambda*r, 0.d0, 2, 2, bi, ncalc)
                     do ik = 0, 1
                        bi(ik) = dexp(-az)*bi(ik)
ccc                        bi(ik) = dexp(az)*bi(ik)
                     end do
                  end if
                  zgrad_psi = zgrad_psi - alambda*bi(1)
     2                          *(z - zk(kbod))/(r)
ccc                  zgrad_psi = zgrad_psi + alambda*bi(1)
ccc     2                          *(z - zk(kbod))/(r)
                  psi = psi + bi(0)
            end do 
           elseif (itest.eq.2) then
c
c An example of rotating concentric cylinders
            r0 = ak(1)
            r1 = ak(0)
            V0 = -2.d0
            V1 = 2.d0
c
c   set coefficient matrix
            az = alambda*r0
               if (az.lt.1000.d0) then
                  call RKBESL (az, 0.d0, 2, 1, bk, ncalc)
                  call RIBESL (az, 0.d0, 2, 1, bi, ncalc)
                 else
                  call RKBESL (az, 0.d0, 2, 2, bk, ncalc)
                  call RIBESL (az, 0.d0, 2, 2, bi, ncalc)
                  do ik = 0, 1
                     bk(ik) = dexp(-az)*bk(ik)
                     bi(ik) = dexp(-az)*bi(ik)
                  end do
               end if
            a11 = -bk(1)
            a12 = bi(1)
            az = alambda*r1
               if (az.lt.1000.d0) then
                  call RKBESL (az, 0.d0, 2, 1, bk, ncalc)
                  call RIBESL (az, 0.d0, 2, 1, bi, ncalc)
                 else
                  call RKBESL (az, 0.d0, 2, 2, bk, ncalc)
                  call RIBESL (az, 0.d0, 2, 2, bi, ncalc)
                  do ik = 0, 1
                     bk(ik) = dexp(-az)*bk(ik)
                     bi(ik) = dexp(-az)*bi(ik)
                  end do
               end if
            a21 = -bk(1)
            a22 = bi(1)
c
c   set rhs
            b1 = -V0/alambda
            b2 = -V1/alambda
c
c   solve, find A and B where psi = A K_0 (lambda r) + B I_0 (lambda r)
            det = a11*a22 - a12*a21
            B = -(b1*a21-b2*a11)/det
            A = (b1 - a12*B)/a11
ccc            call prin2 (' a11 = *', a11, 1)
ccc            call prin2 (' a12 = *', a12, 1)
ccc            call prin2 (' a21 = *', a21, 1)
ccc            call prin2 (' a22 = *', a22, 1)
ccc            call prin2 (' b1 = *', b1, 1)
ccc            call prin2 (' b2 = *', b2, 1)
ccc            call prin2 ('  A = *', A, 1)
ccc            call prin2 ('  B = *', B, 1)
ccc            stop
c
c   Now construct solution at specified point
            r = cdabs(z)
            az = alambda*r
               if (az.lt.1000.d0) then
                  call RKBESL (az, 0.d0, 2, 1, bk, ncalc)
                  call RIBESL (az, 0.d0, 2, 1, bi, ncalc)
                 else
                  call RKBESL (az, 0.d0, 2, 2, bk, ncalc)
                  call RIBESL (az, 0.d0, 2, 2, bi, ncalc)
                  do ik = 0, 1
                     bk(ik) = dexp(-az)*bk(ik)
                     bi(ik) = dexp(-az)*bi(ik)
                  end do
               end if
             psi = A*bk(0) + B*bi(0)
             zgrad_psi = alambda*(-A*bk(1) + B*bi(1))
     1                          *z/r
           elseif (itest.eq.3) then
             zgrad_psi = dcmplx(0.d0,-1.d0)  
          end if
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine SOLVE (maxl, k0, k, nd, nbk, rhs, soln, alpha1, alpha2, 
     1                  Palpha2, sigma1, sigma2, dsdth, rkappa, zf,  
     1                  wsave, rwork, lrwork, iwork, liwork, aK0_k,
     3                  ds_k)
c---------------
c
      implicit real*8 (a-h,o-z)
      external MATVEC, MSOLVE, MATVEC_NOPREC
c
c  System
      dimension soln(*), rhs(*), sigma1(nbk), sigma2(nbk), aK0_k(k)
c
c  DGMRES work arrays
      dimension rwork(lrwork), iwork(liwork)
c
c  Arrays to calculate E o alpha
      dimension alpha1(nbk), alpha2(nbk), Palpha2(nbk), dsdth(nbk),
     1          rkappa(nbk), ds_k(k0:k)
      complex*16 zf(nbk), wsave(*)
c
c  Timings
      real*4 timep(2), etime
c
         pi = 4.d0*datan(1.d0)
c
c  solve linear system using GMRES.
c     parameters for DGMRES
         itol = 0
         tol = 1.0d-14
         isym = 0
         iwork(1) = maxl
ccc         iwork(1) = 50
         do i=2,liwork
            iwork(i) = 0
         enddo
c
c  Preconditioner flag
         iwork(4) = 0
c
c  Restart flag
         iwork(5) = 100      
c
c     provide initial guess soln
         norder = 2*nbk 
ccc         norder = 2*nbk + k
         call PRINF (' SYSTEM SIZE = *', norder, 1)
ccc         norder = nbk
         do i=1,norder
            soln(i) = rhs(i)
         enddo
c
         t0 = etime(timep)
         call DGMRES (norder, rhs, soln, nelt, ia, ja, a, isym,
     1               MATVEC, MSOLVE, itol, tol, itmax, iter,   
     1               err, ierr, 6, sb, sx, rwork, lrwork, iwork, 
     1               liwork, rw, iw)
         call Prin2 (' after  solve, err = *', err, 1)
         call PrinF (' after  solve, ierr = *', ierr, 1)
         call PRINI (6,13)
         call PRINF ('  # GMRES ITERATIONS = *',iter,1)
         if (ierr.ne.0) then
            call PRINF ('  SOMETHING WRONG IN GMRES, IERR = *',ierr,1)
            call PRINF ('  iwork = *',iwork,10)
            call PRINF ('  lrwork = *', lrwork, 1)
            stop
           else
            t1 = etime(timep)
            tsec = t1 - t0
            call PRIn2 (' time in solve = *', tsec, 1)
c
c  unpack soln into U
            istart = 0
            do i = 1, nbk
               alpha1(i) = soln(i)
               alpha2(i) = soln(nbk+i)
            end do  
            call PRIN2 (' alpha1 = *', alpha1, nbk)
            call prin2 (' alpha2 = *', alpha2, nbk)
c
c check int density
            ist = 0
            do kbod = k0, k
               call PRINF (' kbod = *', kbod, 1)
ccc            call prin2 ('     alpha1 = *', alpha1(ist+1), nd)
ccc            call prin2 ('     alpha2 = *', alpha2(ist+1), nd)
               do i = 1, nd
                 zf(i) = alpha1(ist+i)*dsdth(ist+i)
               end do
               call DCFFTF (nd, zf, wsave)
               if (kbod.ne.0) then
                  aK0_k(kbod) = 2.d0*pi*dreal(zf(1))/nd
                  aK0_k(kbod) = aK0_k(kbod)/ds_k(kbod)
               end if
               ai = 2.d0*pi*dreal(zf(1))/nd
               call prin2 (' int alpha1 ds = *', ai, 1)
               do i = 1, nd
                 zf(i) = alpha2(ist+i)*dsdth(ist+i)
               end do
               call DCFFTF (nd, zf, wsave)
               ai = 2.d0*pi*dreal(zf(1))/nd
               call prin2 (' int alpha2 ds = *', ai, 1)
               do i = 1, nd
                  alpha1(ist+i) = alpha1(ist+i)/ds_k(kbod)
                  alpha2(ist+i) = alpha2(ist+i)/ds_k(kbod)
               end do
               ist = ist+nd
            end do
ccc            call prin2 (' alpha1 rescaled = *', alpha1, nbk)
ccc            call prin2 (' alpha2 rescaled = *', alpha2, nbk)
ccc            do kbod = 1, k
ccc               aK0_k(kbod) = soln(2*nbk+kbod)
ccc            end do
ccc            call prin2 ('ak0_k in solve = *', aK0_k, k)
c
c  Calculate E o alpha
            call OPERATE_E (k0, k, nd, nbk, alpha1, alpha2, dsdth,  
     1                      rkappa, sigma1, sigma2, Palpha2, zf, 
     2                      wsave)
c
            open (unit = 22, file = 'alpha.dat')
            do i = 1, nbk
               write(22,'(e20.13,$)')(alpha1(i))
               write (22,'(a)')  ''
            end do
            do i = 1, nbk
               write(22,'(e20.13,$)')(alpha2(i))
               write (22,'(a)')  ''
            end do
            close (22)
            open (unit = 23, file = 'sigma.dat')
            do i = 1, nbk
               write(23,'(e20.13,$)')(sigma1(i))
               write (23,'(a)')  ''
            end do
            do i = 1, nbk
               write(23,'(e20.13,$)')(sigma2(i))
               write (23,'(a)')  ''
            end do
            close (23)
ccc            do i = 1, nbk
ccc               zf(i) = sigma1(i)*dsdth(i)
ccc            end do
ccc            call DCFFTF (nd, zf, wsave)
ccc            far_field = dreal(zf(1))/nd
ccc            call prin2 (' far-field log growth = *', far_field, 1)
ccc            ist = 0
ccc            do kbod = k0, k
ccc               call PRINF (' kbod = *', kbod, 1)
ccc               call PRIN2 ('    sigma1 = *', sigma1(ist+1), nd)
ccc               call PRIN2 ('    sigma2 = *', sigma2(ist+1), nd)
ccc               ist = ist + nd
ccc            end do
            call prin2 (' sigma1 = *', sigma1, nbk)
            call prin2 (' sigma2 = *', sigma2, nbk)
ccc            call prin2 (' alog_k = *', alog_k, k)
         end if
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GET_SOL_GRID (k0, k, nd, nbk, nx, ny, alambda,  
     1                         sigma1, sigma2, z, zn, dsdth, x_grid, 
     2                         y_grid, igrid, u, zf, zf2, wsave)
c---------------
      implicit real*8 (a-h,o-z)
      dimension u(nx,ny), sigma1(nbk), sigma2(nbk), dsdth(nbk),
     1          igrid(nx,ny), G(2)
      dimension x_grid(nx,ny), y_grid(nx,ny)
      complex*16 zn(nbk), z_tar, z(nbk), eye, zf(nbk), zf2(nbk), 
     1           wsave(*), zdis, z_src, zn_src
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
         do i = 1, nx
            do j = 1, ny
               if (igrid(i,j).eq.1) then 
                  z_tar = x_grid(i,j) + eye*y_grid(i,j)
                  psi = 0.d0
                  ist = 0
                  do kbod = k0, k
                     do l = 1, nd
                        z_src = z(ist+l)
                        zn_src = zn(ist+l)
                        call KERNEL_PSI (alambda, z_tar, z_src, zn_src, 
     1                                   G)
ccc                        call prin2 (' G = *', G, 2)
                        zf(l) = sigma1(ist+l)*G(1) + 
     1                             sigma2(ist+l)*G(2)
                        zf(l) = zf(l)*dsdth(ist+l)
                     end do
                     call DCFFTF (nd, zf, wsave)
                     psi = psi + 2.d0*pi*(zf(1) )/nd
                     ist = ist + nd
                  end do
                  u(i,j) = psi
ccc
cccc for uniform flow
ccc                  u(i,j) = psi + y_grid(i,j)
                 else
                  u(i,j) = 0.d0
               end if
            end do
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GET_SOL_POLAR_GRID (k0, k, nd, nbk, nx, ny, alambda,  
     1                         sigma1, sigma2, z, zn, dsdth, x_grid, 
     2                         y_grid, igrid, u, zf, zf2, wsave)
c---------------
      implicit real*8 (a-h,o-z)
      dimension u(nx,ny), sigma1(nbk), sigma2(nbk), dsdth(nbk),
     1          igrid(nx,ny), G(2)
      dimension x_grid(nx,ny), y_grid(nx,ny)
      complex*16 zn(nbk), z_tar, z(nbk), eye, zf(nbk), zf2(nbk), 
     1           wsave(*), zdis, z_src, zn_src
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
         do i = 1, nx
            do j = 1, ny
                  r = x_grid(i,j)
                  th = y_grid(i,j)
                  x = r*dcos(th)
                  y = r*dsin(th)
                  z_tar = x + eye*y
                  psi = 0.d0
                  ist = 0
                  do kbod = k0, k
                     do l = 1, nd
                        z_src = z(ist+l)
                        zn_src = zn(ist+l)
                        call KERNEL_PSI (alambda, z_tar, z_src, zn_src, 
     1                                   G)
ccc                        call prin2 (' G = *', G, 2)
                        zf(l) = sigma1(ist+l)*G(1) + 
     1                             sigma2(ist+l)*G(2)
                        zf(l) = zf(l)*dsdth(ist+l)
                     end do
                     call DCFFTF (nd, zf, wsave)
                     psi = psi + 2.d0*pi*(zf(1) )/nd
                     ist = ist + nd
                  end do
                  u(i,j) = psi 
            end do
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GET_GRAD_PNT (k0, k, nd, nbk, alambda, sigma1, sigma2, 
     1                         z, zn, dsdth, zf1, zf2, wsave, z_tar, 
     2                         zgrad_pnt, zk, aK0_k)
c---------------
      implicit real*8 (a-h,o-z)
      dimension grad_kern(2,2)
      dimension sigma1(nbk), sigma2(nbk), dsdth(nbk), aK0_k(k), bi(0:2)
      complex*16 zn(nbk), z_tar(100), z(nbk), eye, zf1(nbk), zf2(nbk), 
     1           wsave(*), zgrad_K0, zgrad_pnt(100), z_src, zn_src, 
     2           zk(k0:k)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
         open (unit = 15, file = 'target_points.m')
c
ccc         z_tar = -0.1d0 + 0.05d0*eye
ccc         z_tar = -0.225d0 + 0.1875d0*eye
ccc         z_tar = -0.3d0 + 0.225d0*eye
ccc         z_tar = -0.39d0 + 0.33d0*eye
         call DCFFTI (nd, wsave)
         call PRIN2 (' ak = *', aK0_k, k)
         call prin2 (' zk = *', zk(1), 2*k)
c
c  mesh spacing
         dth = 2.d0*pi/nd
c
c np = number of test points
         np = 20
c
c Integrate
         do ip = 1, np
            zgrad_pnt(ip) = 0.d0
c
c error sample points for example 1
ccc         dphi = 2.d0*pi/np
ccc         rad = 0.2d0
ccc            z_tar(ip) = rad*dcmplx(dcos((ip-1)*dphi),dsin((ip-1)*dphi))
c
c error sample points for example 2
            dphi = 2.d0*pi/np
            phi = (ip-1)*dphi
            phi = 7.d0*pi/8.d0 + (ip-1)*dphi
            z_tar(ip) = dcmplx(0.01*dcos(phi), 0.01*dsin(phi))
ccc            z_tar(ip) = dcmplx(-.15+.05*dcos(phi), .11+0.05*dsin(phi))
ccc            z_tar(ip) = dcmplx(0.6*dcos(phi), 0.6*dsin(phi))
ccc            z_tar(ip) = dcmplx(.625+.1*dcos(phi), 0.1*dsin(phi))
ccc            call prin2 (' z_tar = *', z_tar, 2)
            ist = 0
            do kbod = k0, k
               do i = 1, nd
                  z_src = z(ist+i)
                  zn_src = zn(ist+i)
                  call GRAD_KERNEL (alambda, z_tar(ip), z_src, zn_src, 
     1                              grad_kern)
                  zf1(i) = grad_kern(1,1)*sigma1(ist+i)
     1                     + eye*grad_kern(2,1)*sigma1(ist+i)
                  zf1(i) = zf1(i)*dsdth(ist+i)
                  zf2(i) = grad_kern(1,2)*sigma2(ist+i)
     1                     + eye*grad_kern(2,2)*sigma2(ist+i)
                  zf2(i) = zf2(i)*dsdth(ist+i)
               end do
               call DCFFTF (nd, zf1, wsave)
               call DCFFTF (nd, zf2, wsave)
               zgrad_pnt(ip) = zgrad_pnt(ip)
     1                          + 2.d0*pi*(zf1(1) + zf2(1))/nd
               ist = ist + nd
            end do
c
c add on K0 sources
            do lbod = 1, k
               r = cdabs(z_tar(ip)- zk(lbod))
ccc               call prin2 (' r = *', r, 1)
               az = alambda*r
               if (az.lt.1000.d0) then
                  call RKBESL (alambda*r, 0.d0, 2, 1, bi, ncalc)
                 else
                  call RKBESL (alambda*r, 0.d0, 2, 2, bi, ncalc)
                  do ik = 0, 1
                     bi(ik) = dexp(-az)*bi(ik)
                  end do
               end if
               zgrad_pnt(ip) = zgrad_pnt(ip) - aK0_k(lbod)*alambda*bi(1)
     2                          *(z_tar(ip) - zk(lbod))/(r)
            end do
         end do
c
c  Calculate dpsi/dx1
         call PRIN2 (' ztar = *', z_tar, 2*np)
         call PRIN2 (' zgrad_pnt = *', zgrad_pnt, 2*np)
         call RSCPLOT (z_tar,np,1,15)
         close (15)
c
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine CHECK_ERROR (k0, k, zk, nx, ny, alambda, x_grid,  
     1                        y_grid, igrid, u, z_tar, zgrad_pnt, ak)
c---------------
      implicit real*8 (a-h,o-z)
      complex*16 zgrad_psi
      dimension x_grid(nx,ny), y_grid(nx,ny), u(nx,ny), bi(0:1)
      dimension igrid(nx,ny), ak(k0:k)
      complex*16 z, eye, zk(k0:k), z_tar(100), zgrad_pnt(100)
c
         eye = dcmplx(0.d0,1.d0)
c
         err = 0.d0
         cmax = -1.d10
         cmin = 1.d10
         umax = 0.d0
         do i = 1, nx
            do j = 1, ny
               if (igrid(i,j).eq.1) then
                  z = x_grid(i,j) + eye*y_grid(i,j)
                  r = cdabs(z - zk(k0))
                  call EXACT_SOL (k0, k, z, zk, alambda, psi, zgrad_psi,
     1                            ak(k0))
                  u_exact = psi
                  umax = max(umax,dabs(u(i,j)))
                  diff = u_exact-u(i,j)
                  cmax = max (cmax,diff)
                  cmin = min (cmin,diff)
ccc                  call prin2 (' u_exact = *', u_exact,1)
ccc                  call PRIN2 ('    u  = *', u(i,j), 1)
ccc                  call PRIN2 ('    diff  = *', diff, 1)
                  err = max(err,dabs(u_exact-u(i,j)))
               end if
            end do
         end do
         call prin2 (' minimum difference = *', cmin, 1)
         call prin2 (' maximum difference = *', cmax, 1)
         call PRIN2 (' ERROR IN SOLUTION = *', dabs(cmin-cmax), 1)
         relerr = dabs(cmin-cmax)/umax
         call PRIN2 (' relative error = *', relerr, 1)
c
c Check gradient at point
         rel_err_max = 0.d0
         abs_err_max = 0.d0
         np = 20
         do ip = 1, np
            call EXACT_SOL (k0, k, z_tar(ip), zk, alambda, psi, 
     1                      zgrad_psi, ak(k0))
ccc            call PRIN2 (' exact grad at point = *', zgrad_psi, 2)
ccc            call PRIN2 (' calc grad at point = *', zgrad_pnt(ip), 2)
            diff = cdabs(zgrad_psi-zgrad_pnt(ip))
ccc            call prin2 (' diff = *', diff, 1)
            rel_err_max = max(rel_err_max,diff/cdabs(zgrad_psi))
            abs_err_max = max(abs_err_max,diff)
         end do
         call PRIN2 (' max absolute error = *', abs_err_max, 1)
         call PRIN2 (' max relative error = *', rel_err_max, 1)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine MSOLVE (N,x,u,nelt,ia,ja,a,isym,rwork,iwork)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension s(N),u(N)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GETMAT1 (k, nd, nbk, u, w, amat, rwork, iwork)
c---------------
      implicit real*8 (a-h,o-z)
      character norm
      dimension u(2*nbk+k), w(2*nbk+k), amat(2*nbk,2*nbk), 
     1          rwork(*), iwork(*)
      dimension a(10,10), gwork(100000), ipiv(100000), iwk(100000)
      complex*16 eye
c
        write (6,*) 'in getmat1'
        call prini (6,13)
         open (unit = 24, file = 'matrix.dat')
         do i = 1, 2*nbk
            u(i) = 0.d0
         end do
         do i = 1, 2*nbk
            u(i) = 1.d0
            call MATVEC (2*nbk,u,w,nelt,ia,ja,dummy_mat,isym)
            do j = 1, 2*nbk
               write(24,'(e20.13,$)')(w(j))
               write (24,'(a)')  ''
               amat(j,i) = w(j)
            end do
            u(i) = 0.d0
            gwork(i) = 1.d0
         end do
         close (24) 
c
c estimate the condition number
         norm = '1' 
         lda = 2*nbk
ccc         call DGECO (amat, lda, lda, ipiv, RCOND, gwork)
ccc         rcond = 1.d0/rcond
ccc            call PRIN2 (' RCOND = *', rcond, 1)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GETMAT2 (k, nd, nbk, u, w, zf, wsave)
c---------------
      implicit real*8 (a-h,o-z)
      dimension u(2*nbk), w(2*nbk)
      complex*16 eye, zf(nd), wsave(*)
c
         pi = 4.d0*datan(1.d0)
c
         open (unit = 25, file = 'G11.dat')
         open (unit = 26, file = 'G12.dat')
         open (unit = 27, file = 'G21.dat')
         open (unit = 28, file = 'G22.dat')
c
         call DCFFTI (nd, wsave)
c
         kmode = 32
         do i = 1, nbk
            th = 2.d0*pi*(i-1.d0)/nd
            u(i) = dcos(kmode*th)
         end do
         call prin2 (' u = *', u, nbk)
         do i = nbk+1, 2*nbk
            u(i) = 0.d0
         end do
         call MATVEC (2*nbk,u,w,nelt,ia,ja,dummy_mat,isym)
         do j = 1, nbk
            write(25,'(e20.13,$)')(w(j))
            write (25,'(a)')  ''
         end do
         do j = nbk+1, 2*nbk
            write(27,'(e20.13,$)')(w(j))
            write (27,'(a)')  ''
         end do
c
         do i = 1, nbk
            u(i) = 0.d0
         end do
         do i = nbk+1, 2*nbk
            th = 2.d0*pi*(i-1.d0)/nd
            u(i) = dcos(kmode*th)
         end do
         call MATVEC (2*nbk,u,w,nelt,ia,ja,dummy_mat,isym)
         do j = 1, nbk
            write(26,'(e20.13,$)')(w(j))
            write (26,'(a)')  ''
         end do
         do j = nbk+1, 2*nbk
            write(28,'(e20.13,$)')(w(j))
            write (28,'(a)')  ''
         end do
         close (26)
         close (28)    
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine ASSEMBLE_MAT (k0, k, nd, nbk, pmat, G11, G12, G21, G22,
     1                         emat, amat, rkappa, wmat)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension pmat(nbk,nbk), G11(nbk,nbk), G12(nbk,nbk), G21(nbk,nbk),
     1          G22(nbk,nbk), amat(2*nbk,2*nbk), emat(2*nbk,2*nbk),
     2          wmat(2*nbk,2*nbk), rkappa(*)
c
c
         do i = 1, nbk
            do j = 1, nbk
               amat(i,j) = G11(i,j)
               amat(i,nbk+j) = G12(i,j)
               amat(nbk+i,j) = G21(i,j)
               amat(nbk+i,nbk+j) = G22(i,j)
               if (i.eq.j) then
                  emat(i,j) = 2.d0
                 else
                  emat(i,j) = 0.d0
               end if
               emat(i,nbk+j) = 4.d0*rkappa(i)*pmat(i,j)
               emat(nbk+i,j) = 0.d0
               emat(nbk+i,nbk+j) = 2.d0*pmat(i,j)
            end do
         end do
c
c  assemble AE
         do i = 1, 2*nbk
            do j = 1, 2*nbk
               wmat(i,j) = 0.d0
               do kk = 1,2*nbk
                  wmat(i,j) = wmat(i,j) + amat(i,kk)*emat(kk,j)
               end do
            end do
         end do
c
c add in +/- I depending on exterior or interior
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               wmat(ist+i,ist+i) = 1.d0 + wmat(ist+i,ist+i)
               if (kbod.eq.0) then
                  wmat(nbk+ist+i,nbk+ist+i) 
     1                = 1.d0 + wmat(nbk+ist+i,nbk+ist+i)
                 else
                  wmat(nbk+ist+i,nbk+ist+i) 
     1                = -1.d0 + wmat(nbk+ist+i,nbk+ist+i)
               end if
            end do
            ist = ist+nd
         end do
c
         open (unit = 24, file = 'matrix2.dat')
         do i = 1, 2*nbk
            do j = 1, 2*nbk
               write(24,'(e20.13,$)')(wmat(i,j))
               write (24,'(a)')  ''
            end do
         end do
         close (24) 
c
      return
      end
c
c
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine BUILD_MAT (N,wmat)
c---------------
c
c hopefully speedier, some of the computational redundancies are 
c removed
      implicit real*8 (a-h,o-z)
      parameter (ndmax = 512, kmax = 4, nmax = kmax*ndmax)
c
      dimension wmat(N**2)
      dimension pmat(nmax**2), emat(4*nmax**2), amat(4*nmax**2)
      dimension rkappa(nmax), dsdth(nmax), ds_k(kmax)
      complex*16 z(nmax), dz(nmax), zn(nmax), zs(nmax), zk(kmax)
c
c  kernels and integral operators
      dimension G11(nmax**2), G12(nmax**2), G21(nmax**2), 
     1          G22(nmax**2)
c
c  work arrays for Alpert's rule
      complex*16 zf2(nmax)
      dimension f_vec1(15*ndmax**2), f_vec2(15*ndmax**2)
      dimension uweights(16), vnodes(16)
      dimension ds_vec1(15*nmax), ds_vec2(15*nmax)
      complex*16 z_vec1(15*nmax), z_vec2(15*nmax)
      complex*16 zn_vec1(15*nmax), zn_vec2(15*nmax)
      complex*16 zs_vec1(15*nmax), zs_vec2(15*nmax)
c
c log singularities for multiply-connected domains
      dimension aK0_k(kmax), bi(0:2)
c
c  table of Bessel functions and other special functions
      dimension besselK0(nmax*nmax), besselK1(nmax*nmax)
      dimension c0(nmax*nmax), c1(nmax*nmax), c2(nmax*nmax),
     1          c3(nmax*nmax), c4(nmax*nmax), c5(nmax*nmax)
c
c Other arrays
      complex*16 zf(nmax), wsave(20*nmax)
      REAL*4 TIMEP(2), ETIME
c
c complex work variables
      complex*16 z_tar, z_src, zn_tar, zn_src, zs_tar, zs_src, zdis, 
     1           zgrad_K0
c
c Common Blocks
      common /geometry_info/ z, dz, zn, zs, zk, rkappa, dsdth, ds_k
      common /size_info/ k0, k, nd, nbk
      common /other_parameters/ alambda
      common /alpert/ uweights, vnodes, numpts, ishift, iorder, iquad
      common /geo_nodes/ ds_vec1, ds_vec2, z_vec1, z_vec2
      common /geo_nodes2/ zn_vec1, zn_vec2, zs_vec1, zs_vec2
      common /bessel_fns/ besselK0, besselK1
      common /coefficients/ c0, c1, c2, c3, c4, c5
c
         pi = 4.0d0*datan(1.0d0)
c
c Precondition
         call BUILD_FOR_MAT (k0, k, nd, nbk, dsdth, pmat, zf, wsave)
         call NODE_MAT (nd, numpts, vnodes, f_vec1,    
     1                  f_vec2, zf, zf2, wsave) 
         call ALPERT_MAT (k0, k, nd, nbk, alambda, z, z_vec1,  
     1               z_vec2, dsdth, ds_vec1, ds_vec2, zn, zn_vec1, 
     2               zn_vec2, zs, zs_vec1, zs_vec2, f_vec1, f_vec2,
     3               c0, c1, c2, c3, c4, c5, G11, G12, G21, 
     6               G22, ds_k)
         call ASSEMBLE_MAT (k0, k, nd, nbk, pmat, G11, G12, G21, G22,
     1                      emat, amat, rkappa, wmat)
c
ccc        call prin2 (' G12 = *', G12, nbk)
ccc         err = 0.d0
ccc         do i = 1, nbk
ccc            err = max(err,dabs(G12(i)-G12_O(i)))
ccc         end do
ccc         call prin2 (' err g12 = *', err, 1)
ccc         stop
ccc         call prin2 (' G11 = *', G11, nbk)
ccc         call prin2 (' G12 = *', G12, nbk)
ccc         call prin2 (' G21 = *', G21, nbk)
ccc         call prin2 (' G22 = *', G22, nbk)
ccc         ist = 0
ccc         do kbod = k0, k
ccc            do i = 1, nd
ccc               w(ist+i) = w(ist+i) + G11(ist+i) + G12(ist+i)
ccc               if (kbod.eq.0) then
ccc                  w(nbk+ist+i) = w(nbk+ist+i) + G21(ist+i) + G22(ist+i)
ccc                else
ccc                  w(nbk+ist+i) = -w(nbk+ist+i) + G21(ist+i) + G22(ist+i)
ccc               endif
ccc            end do
ccc            ist = ist + nd
ccc         end do
c            
ccc         do i = 1, nbk
ccc            w(i) = w(i) + G11(i) + G12(i)
ccc            w(nbk+i) = w(nbk+i) + G21(i) + G22(i)
ccc         end do
ccc         call prin2 (' wn before sources = *', w, nbk)
ccc         call prin2 (' ws = *', w(nbk+1), nbk)
c
c Add on singularities from holes
ccc         ist = 0
ccc         do kbod = k0, k
ccc            kk = kbod+1-k0
ccc            do i = 1, nd
ccc               zgrad_K0 = 0.d0
ccc               do lbod = 1, k
ccc                  ll = lbod+1-k0
ccc                  ds_fact = ds_k(kk)/ds_k(ll)
ccc                  r = cdabs(z(ist+i) - zk(ll))
ccc                  az = alambda*r
ccc                  if (az.lt.1000.d0) then
ccc                     call RKBESL (alambda*r, 0.d0, 2, 1, bi, ncalc)
ccc                    else
ccc                     call RKBESL (alambda*r, 0.d0, 2, 2, bi, ncalc)
ccc                     do ik = 0, 1
ccc                        bi(ik) = dexp(-az)*bi(ik)
ccc                     end do
ccc                  end if
ccc                  call prin2 (' ds_fact = *', ds_fact, 1)
ccc                  zgrad_K0 = zgrad_K0 - aK0_k(lbod)*alambda*bi(1)
ccc     2                             *(z(ist+i) - zk(ll))*ds_fact/(r)
ccc               end do
ccc               w(ist+i) = w(ist+i) + dreal(zgrad_K0*dconjg(zn(ist+i)))
ccc               w(nbk+ist+i) = w(nbk+ist+i) 
ccc     1                         + dreal(zgrad_K0*dconjg(zs(ist+i)))
ccc            end do
ccc            ist = ist + nd
ccc         end do
ccc         call prin2 (' wn = *', w, nbk)
ccc         call prin2 (' ws = *', w(nbk+1), nbk)
c
c  Constraints
ccc         do kbod = 1, k
ccc            w(2*nbk+kbod) = 0.d0
ccc            do i = 1,nd
ccc               index = (kbod-k0)*nd + i
ccc               w(2*nbk+kbod) = w(2*kbk+kbod) + alpha1(index)
ccc            end do
ccc         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine MATVEC (N,u,w,nelt,ia,ja,dummy_mat,isym)
c---------------
c
c hopefully speedier, some of the computational redundancies are 
c removed
      implicit real*8 (a-h,o-z)
      parameter (ndmax = 512, kmax = 4, nmax = kmax*ndmax)
c
      dimension rkappa(nmax), dsdth(nmax), ds_k(kmax)
      complex*16 z(nmax), dz(nmax), zn(nmax), zs(nmax), zk(kmax)
c
      dimension u(N),w(N)
      dimension alpha1(nmax), alpha2(nmax), Palpha2(nmax)
      dimension Ealpha1(nmax), Ealpha2(nmax)
      dimension ak(kmax), bk(kmax), phi_k(kmax), ncyc(kmax)
c
c  kernels and integral operators
      dimension G11(nmax), G12(nmax), G21(nmax), G22(nmax)
c
c  work arrays for Alpert's rule
      complex*16 zf2(nmax)
      dimension Ealph1_vec1(15*nmax), Ealph1_vec2(15*nmax), 
     1          Ealph2_vec1(15*nmax), Ealph2_vec2(15*nmax)
      dimension uweights(16), vnodes(16)
      dimension ds_vec1(15*nmax), ds_vec2(15*nmax)
      complex*16 z_vec1(15*nmax), z_vec2(15*nmax)
      complex*16 zn_vec1(15*nmax), zn_vec2(15*nmax)
      complex*16 zs_vec1(15*nmax), zs_vec2(15*nmax)
c
c log singularities for multiply-connected domains
      dimension aK0_k(kmax), bi(0:2)
c
c  table of Bessel functions and other special functions
      dimension besselK0(nmax*nmax), besselK1(nmax*nmax)
      dimension c0(nmax*nmax), c1(nmax*nmax), c2(nmax*nmax),
     1          c3(nmax*nmax), c4(nmax*nmax), c5(nmax*nmax)
c
c Other arrays
      complex*16 zf(nmax), wsave(20*nmax)
      REAL*4 TIMEP(2), ETIME
c
c complex work variables
      complex*16 z_tar, z_src, zn_tar, zn_src, zs_tar, zs_src, zdis, 
     1           zgrad_K0
c
c Common Blocks
      common /geometry_info/ z, dz, zn, zs, zk, rkappa, dsdth, ds_k
      common /size_info/ k0, k, nd, nbk
      common /other_parameters/ alambda
      common /ellipse_info/ ak, bk, phi_k, ncyc
      common /alpert/ uweights, vnodes, numpts, ishift, iorder, iquad
      common /geo_nodes/ ds_vec1, ds_vec2, z_vec1, z_vec2
      common /geo_nodes2/ zn_vec1, zn_vec2, zs_vec1, zs_vec2
      common /bessel_fns/ besselK0, besselK1
      common /coefficients/ c0, c1, c2, c3, c4, c5
c
         pi = 4.0d0*datan(1.0d0)
c
c  Identity operator
         do i = 1, 2*nbk
            w(i) = u(i)
         end do
c
c  Calculate E o alpha
         do i = 1, nbk
            alpha1(i) = u(i)
            alpha2(i) = u(nbk+i)
         enddo
ccc         do kbod = 1, k
ccc            aK0_k(kbod) = u(2*nbk+kbod)
ccc         end do
ccc         call PRIN2 (' in matvec, zk = *', zk, 4)
c
c  Calculate int_Gamma_k alpha_1 ds
         ist = nd
         do kbod = 1, k
            aK0_k(kbod) = 0.d0
ccc            call prin2 (' alpha1 = *', alpha1(ist+1), nd)
            do i = 1, nd
               aK0_k(kbod) = aK0_k(kbod) + alpha1(ist+i)*dsdth(ist+i)
            end do
            aK0_k(kbod) = aK0_k(kbod)*2.d0*pi/nd
            ist = ist + nd
         end do
ccc         call prin2 ('aK0_k = *', aK0_k, k)
c
c Precondition
         call OPERATE_E (k0, k, nd, nbk, alpha1, alpha2, dsdth,  
     1                   rkappa, Ealpha1, Ealpha2, Palpha2, zf, 
     2                   wsave)
ccc         call prin2 (' Ealpha1 = *', Ealpha1, nbk)
ccc         call prin2 (' Ealpha2 = *', Ealpha2, nbk)
c
c  Interpolate layer densities to nodes
         call REAL_AT_NODES (k0, k, nd, nbk, vnodes, numpts, Ealpha1,   
     1                       Ealph1_vec1, Ealph1_vec2, zf, zf2, wsave)
         call REAL_AT_NODES (k0, k, nd, nbk, vnodes, numpts, Ealpha2,   
     1                       Ealph2_vec1, Ealph2_vec2, zf, zf2, wsave)
c
         call ALPERT_QUAD_GEN (k0, k, nd, nbk, alambda, z, z_vec1,    
     2            z_vec2, dsdth, ds_vec1, ds_vec2, zn, zn_vec1, zn_vec2, 
     3            zs, zs_vec1, zs_vec2, Ealpha1, Ealph1_vec1,  
     4            Ealph1_vec2, Ealpha2, Ealph2_vec1, Ealph2_vec2, c0,  
     5            c1, c2, c3, c4, c5, G11, G12, G21, G22, ds_k)
ccc        call prin2 (' G12 = *', G12, nbk)
ccc         err = 0.d0
ccc         do i = 1, nbk
ccc            err = max(err,dabs(G12(i)-G12_O(i)))
ccc         end do
ccc         call prin2 (' err g12 = *', err, 1)
ccc         stop
ccc         call prin2 (' G11 = *', G11, nbk)
ccc         call prin2 (' G12 = *', G12, nbk)
ccc         call prin2 (' G21 = *', G21, nbk)
ccc         call prin2 (' G22 = *', G22, nbk)
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               w(ist+i) = w(ist+i) + G11(ist+i) + G12(ist+i)
               if (kbod.eq.0) then
                  w(nbk+ist+i) = w(nbk+ist+i) + G21(ist+i) + G22(ist+i)
                else
                  w(nbk+ist+i) = -w(nbk+ist+i) + G21(ist+i) + G22(ist+i)
               endif
            end do
            ist = ist + nd
         end do
c            
ccc         do i = 1, nbk
ccc            w(i) = w(i) + G11(i) + G12(i)
ccc            w(nbk+i) = w(nbk+i) + G21(i) + G22(i)
ccc         end do
ccc         call prin2 (' wn before sources = *', w, nbk)
ccc         call prin2 (' ws = *', w(nbk+1), nbk)
c
c Add on singularities from holes
ccc         ist = 0
ccc         do kbod = k0, k
ccc            kk = kbod+1-k0
ccc            do i = 1, nd
ccc               zgrad_K0 = 0.d0
ccc               do lbod = 1, k
ccc                  ll = lbod+1-k0
ccc                  ds_fact = ds_k(kk)/ds_k(ll)
ccc                  r = cdabs(z(ist+i) - zk(ll))
ccc                  az = alambda*r
ccc                  if (az.lt.1000.d0) then
ccc                     call RKBESL (alambda*r, 0.d0, 2, 1, bi, ncalc)
ccc                    else
ccc                     call RKBESL (alambda*r, 0.d0, 2, 2, bi, ncalc)
ccc                     do ik = 0, 1
ccc                        bi(ik) = dexp(-az)*bi(ik)
ccc                     end do
ccc                  end if
ccc                  call prin2 (' ds_fact = *', ds_fact, 1)
ccc                  zgrad_K0 = zgrad_K0 - aK0_k(lbod)*alambda*bi(1)
ccc     2                             *(z(ist+i) - zk(ll))*ds_fact/(r)
ccc               end do
ccc               w(ist+i) = w(ist+i) + dreal(zgrad_K0*dconjg(zn(ist+i)))
ccc               w(nbk+ist+i) = w(nbk+ist+i) 
ccc     1                         + dreal(zgrad_K0*dconjg(zs(ist+i)))
ccc            end do
ccc            ist = ist + nd
ccc         end do
ccc         call prin2 (' wn = *', w, nbk)
ccc         call prin2 (' ws = *', w(nbk+1), nbk)
c
c  Constraints
ccc         do kbod = 1, k
ccc            w(2*nbk+kbod) = 0.d0
ccc            do i = 1,nd
ccc               index = (kbod-k0)*nd + i
ccc               w(2*nbk+kbod) = w(2*kbk+kbod) + alpha1(index)
ccc            end do
ccc         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine MATVEC_CHECK (N,u,w,nelt,ia,ja,dummy_mat,isym)
c---------------
c
c hopefully speedier, some of the computational redundancies are 
c removed
      implicit real*8 (a-h,o-z)
      parameter (ndmax = 512, kmax = 4, nmax = kmax*ndmax)
c
      dimension rkappa(nmax), dsdth(nmax), ds_k(kmax)
      complex*16 z(nmax), dz(nmax), zn(nmax), zs(nmax), zk(kmax)
c
      dimension u(N),w(N)
      dimension alpha1(nmax), alpha2(nmax), Palpha2(nmax)
      dimension Ealpha1(nmax), Ealpha2(nmax)
      dimension ak(kmax), bk(kmax), phi_k(kmax), ncyc(kmax)
c
c  kernels and integral operators
      dimension G11(nmax), G12(nmax), G21(nmax), G22(nmax)
c
c  work arrays for Alpert's rule
      complex*16 zf2(nmax)
      dimension Ealph1_vec1(15*nmax), Ealph1_vec2(15*nmax), 
     1          Ealph2_vec1(15*nmax), Ealph2_vec2(15*nmax)
      dimension uweights(16), vnodes(16)
      dimension ds_vec1(15*nmax), ds_vec2(15*nmax)
      complex*16 z_vec1(15*nmax), z_vec2(15*nmax)
      complex*16 zn_vec1(15*nmax), zn_vec2(15*nmax)
      complex*16 zs_vec1(15*nmax), zs_vec2(15*nmax)
c
c log singularities for multiply-connected domains
      dimension aK0_k(kmax), bi(0:2)
c
c  table of Bessel functions and other special functions
      dimension besselK0(nmax*nmax), besselK1(nmax*nmax)
      dimension c0(nmax*nmax), c1(nmax*nmax), c2(nmax*nmax),
     1          c3(nmax*nmax), c4(nmax*nmax), c5(nmax*nmax)
c
c Other arrays
      complex*16 zf(nmax), wsave(20*nmax)
      REAL*4 TIMEP(2), ETIME
c
c complex work variables
      complex*16 z_tar, z_src, zn_tar, zn_src, zs_tar, zs_src, zdis, 
     1           zgrad_K0
c
c Common Blocks
      common /geometry_info/ z, dz, zn, zs, zk, rkappa, dsdth, ds_k
      common /size_info/ k0, k, nd, nbk
      common /other_parameters/ alambda
      common /ellipse_info/ ak, bk, phi_k, ncyc
      common /alpert/ uweights, vnodes, numpts, ishift, iorder, iquad
      common /geo_nodes/ ds_vec1, ds_vec2, z_vec1, z_vec2
      common /geo_nodes2/ zn_vec1, zn_vec2, zs_vec1, zs_vec2
      common /bessel_fns/ besselK0, besselK1
      common /coefficients/ c0, c1, c2, c3, c4, c5
c
         pi = 4.0d0*datan(1.0d0)
c
c  Calculate E o alpha
         do i = 1, nbk
            alpha1(i) = u(i)
            alpha2(i) = u(nbk+i)
         enddo
ccc         do kbod = 1, k
ccc            aK0_k(kbod) = u(2*nbk+kbod)
ccc         end do
ccc         call PRIN2 (' in matvec, zk = *', zk, 4)
ccc         call prin2 ('aK0_k = *', aK0_k, k)
c
c Precondition
         call OPERATE_E (k0, k, nd, nbk, alpha1, alpha2, dsdth,  
     1                      rkappa, Ealpha1, Ealpha2, Palpha2, zf, 
     2                      wsave)
ccc         call OPERATE_E_CHK (k0, k, nd, nbk, alpha1, alpha2, dsdth,  
ccc     1                   rkappa, Ealpha1, Ealpha2, Palpha2, zf, 
ccc     2                   wsave)
ccc         call prin2 (' Ealpha1 = *', Ealpha1, nbk)
ccc         call prin2 (' Ealpha2 = *', Ealpha2, nbk)
ccc         do i = 1, nbk
ccc            Ealpha1(i) = alpha1(i)
ccc            Ealpha2(i) = alpha2(i)
ccc         end do
c
c  Interpolate layer densities to nodes
         call REAL_AT_NODES (k0, k, nd, nbk, vnodes, numpts, Ealpha1,   
     1                       Ealph1_vec1, Ealph1_vec2, zf, zf2, wsave)
         call REAL_AT_NODES (k0, k, nd, nbk, vnodes, numpts, Ealpha2,   
     1                       Ealph2_vec1, Ealph2_vec2, zf, zf2, wsave)
c
         call ALPERT_QUAD_GEN (k0, k, nd, nbk, alambda, z, z_vec1,    
     2            z_vec2, dsdth, ds_vec1, ds_vec2, zn, zn_vec1, zn_vec2, 
     3            zs, zs_vec1, zs_vec2, Ealpha1, Ealph1_vec1,  
     4            Ealph1_vec2, Ealpha2, Ealph2_vec1, Ealph2_vec2, c0,  
     5            c1, c2, c3, c4, c5, G11, G12, G21, G22, ds_k)
ccc        call prin2 (' G12 = *', G12, nbk)
ccc         err = 0.d0
ccc         do i = 1, nbk
ccc            err = max(err,dabs(G12(i)-G12_O(i)))
ccc         end do
ccc         call prin2 (' err g12 = *', err, 1)
ccc         stop
ccc         call prin2 (' G11 = *', G11, nbk)
ccc         call prin2 (' G12 = *', G12, nbk)
ccc         call prin2 (' G21 = *', G21, nbk)
ccc         call prin2 (' G22 = *', G22, nbk)
         ist = 0
         do kbod = k0, k
            do i = 1, nd
ccc               w(ist+i) = w(ist+i) + G11(ist+i) + G12(ist+i)
               w(ist+i) =  G11(ist+i) + G12(ist+i)
               if (kbod.eq.0) then
ccc                  w(nbk+ist+i) = w(nbk+ist+i) + G21(ist+i) + G22(ist+i)
                  w(nbk+ist+i) = G21(ist+i) + G22(ist+i)
                else
ccc                  w(nbk+ist+i) = -w(nbk+ist+i) + G21(ist+i) + G22(ist+i)
                  w(nbk+ist+i) =  G21(ist+i) + G22(ist+i)
               endif
            end do
            ist = ist + nd
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine MATVEC_NOPREC (N,u,w,nelt,ia,ja,dummy_mat,isym)
c---------------
c
c hopefully speedier, some of the computational redundancies are 
c removed
      implicit real*8 (a-h,o-z)
      parameter (ndmax = 512, kmax = 4, nmax = kmax*ndmax)
c
      dimension rkappa(nmax), dsdth(nmax), ds_k(kmax)
      complex*16 z(nmax), dz(nmax), zn(nmax), zs(nmax), zk(kmax)
c
      dimension u(N),w(N)
      dimension alpha1(nmax), alpha2(nmax), Palpha2(nmax)
      dimension sigma1(nmax), sigma2(nmax)
      dimension Dalpha1(nmax), Dalpha2(nmax)
      dimension ak(kmax), bk(kmax), phi_k(kmax), ncyc(kmax)
c
c  kernels and integral operators
      dimension G11(nmax), G12(nmax), G21(nmax), G22(nmax)
c
c  work arrays for Alpert's rule
      complex*16 zf2(nmax)
      dimension Ealph1_vec1(15*nmax), Ealph1_vec2(15*nmax), 
     1          Ealph2_vec1(15*nmax), Ealph2_vec2(15*nmax)
      dimension Esig1_vec1(15*nmax), Esig1_vec2(15*nmax), 
     1          Esig2_vec1(15*nmax), Esig2_vec2(15*nmax)
      dimension uweights(16), vnodes(16)
      dimension ds_vec1(15*nmax), ds_vec2(15*nmax)
      complex*16 z_vec1(15*nmax), z_vec2(15*nmax)
      complex*16 zn_vec1(15*nmax), zn_vec2(15*nmax)
      complex*16 zs_vec1(15*nmax), zs_vec2(15*nmax)
c
c log singularities for multiply-connected domains
      dimension aK0_k(kmax)
c
c  table of Bessel functions and other special functions
      dimension besselK0(nmax*nmax), besselK1(nmax*nmax)
      dimension c0(nmax*nmax), c1(nmax*nmax), c2(nmax*nmax),
     1          c3(nmax*nmax), c4(nmax*nmax), c5(nmax*nmax)
c
c Other arrays
      complex*16 zf(nmax), wsave(20*nmax)
      REAL*4 TIMEP(2), ETIME
c
c complex work variables
      complex*16 z_tar, z_src, zn_tar, zn_src, zs_tar, zs_src, zdis, 
     1           zgrad_log
c
c Common Blocks
      common /geometry_info/ z, dz, zn, zs, zk, rkappa, dsdth, ds_k
      common /size_info/ k0, k, nd, nbk
      common /other_parameters/ alambda
      common /ellipse_info/ ak, bk, phi_k, ncyc
      common /alpert/ uweights, vnodes, numpts, ishift, iorder, iquad
      common /geo_nodes/ ds_vec1, ds_vec2, z_vec1, z_vec2
      common /geo_nodes2/ zn_vec1, zn_vec2, zs_vec1, zs_vec2
      common /bessel_fns/ besselK0, besselK1
      common /coefficients/ c0, c1, c2, c3, c4, c5
c
         pi = 4.0d0*datan(1.0d0)
c
c  Identity operator
         do i = 1, 2*nbk
            w(i) = u(i)
         end do
c
c  Calculate E o alpha
         do i = 1, nbk
            sigma1(i) = u(i)
            sigma2(i) = u(nbk+i)
         enddo
         call DIAG (k0, k, nd, nbk, sigma1, sigma2, dsdth,  
     1                   rkappa, Dalpha1, Dalpha2, Palpha2, zf, 
     2                   wsave)
c
c  Interpolate layer densities to nodes         
         call REAL_AT_NODES (k0, k, nd, nbk, vnodes, numpts, sigma1,   
     1                       Esig1_vec1, Esig1_vec2, zf, zf2, wsave)
         call REAL_AT_NODES (k0, k, nd, nbk, vnodes, numpts, alpha2,   
     1                       Esig2_vec1, Esig2_vec2, zf, zf2, wsave)
c
         call ALPERT_QUAD_GEN (k0, k, nd, nbk, alambda, z, z_vec1,    
     2            z_vec2, dsdth, ds_vec1, ds_vec2, zn, zn_vec1, zn_vec2, 
     3            zs, zs_vec1, zs_vec2, sigma1, Esig1_vec1,  
     4            Esig1_vec2, sigma2, Esig2_vec1, Esig2_vec2, c0,  
     5            c1, c2, c3, c4, c5, G11, G12, G21, G22, ds_k)
ccc        call prin2 (' G12 = *', G12, nbk)
ccc         err = 0.d0
ccc         do i = 1, nbk
ccc            err = max(err,dabs(G12(i)-G12_O(i)))
ccc         end do
ccc         call prin2 (' err g12 = *', err, 1)
ccc         stop
         do i = 1, nbk
            w(i) = Dalpha1(i) + G11(i) + G12(i)
            w(nbk+i) = Dalpha2(i) + G21(i) + G22(i)
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine BUILD_FOR_MAT (k0, k, nd, nbk, dsdth, pmat, zf, wsave)
c---------------
c Constructs Fourier Integration matrix
c
      implicit real*8 (a-h,o-z)
      dimension pmat(nbk,nbk), dsdth(nbk)
      complex*16 zf(nbk), wsave(*), eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
         n = nd/2
c
         call DCFFTI (nd, wsave)
c
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               do j = 1, i-1
                  zf(j) = 0.d0
               end do
               zf(i) = dsdth(ist+i)
               do j = i+1,nd
                  zf(j) = 0.d0
               end do
               call DCFFTF (nd, zf, wsave) 
               do j = 1, nd
                  zf(j) = zf(j)/nd
               end do
               zf(1) = 0.d0
               do kmode = 1, n-1
                  zf(kmode+1) = zf(kmode+1)/(eye*kmode)
                  zf(nd-kmode+1) = -zf(nd-kmode+1)/(eye*kmode)
               end do
               zf(n+1) = 0
               call DCFFTB (nd, zf, wsave)
               asum = dreal(zf(1))
               do j = 1, nd
                  pmat(ist+j,ist+i) = dreal(zf(j)) - asum
               end do
            end do
            ist = ist + nd
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine OPERATE_E (k0, k, nd, nbk, alpha1, alpha2, dsdth,  
     1                      rkappa, Ealpha1, Ealpha2, Palpha2, zf, 
     2                      wsave)
c---------------
c Calculates E o alpha, where
c     E = [2 I  4 kappa P
c            0        2 P ]
c and
c     P o alpha = \int^s alpha(s) ds
c
      implicit real*8 (a-h,o-z)
      dimension alpha1(nbk), alpha2(nbk), dsdth(nbk), rkappa(nbk)
      dimension Palpha2(nbk), Ealpha1(nbk), Ealpha2(nbk), ck(0:100)
      complex*16 zf(nbk), wsave(*), eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
         n = nd/2
c
         call DCFFTI (nd, wsave)
c
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               zf(i) = alpha2(ist+i)*dsdth(ist+i)
            end do
            call DCFFTF (nd, zf, wsave)
            do i = 1, nd
               zf(i) = zf(i)/nd
            end do
            f0 = dreal(zf(1))
ccc            call PRIN2 (' in OPERATE_E, f0 = *', f0, 1)
            zf(1) = 0.d0
            do kmode = 1, n-1
               zf(kmode+1) = zf(kmode+1)/(eye*kmode)
               zf(nd-kmode+1) = -zf(nd-kmode+1)/(eye*kmode)
            end do
            zf(n+1) = 0
            call DCFFTB (nd, zf, wsave)
            asum = dreal(zf(1))
ccc            ftot = f0*2.d0*pi
            do i = 1, nd
               th = 2.d0*pi*(i-1.d0)/nd
               if (kbod.eq.0) then
                  Palpha2(ist+i) = dreal(zf(i)) - asum 
                 else
ccc                  Palpha2(ist+i) = -dreal(zf(i)) + asum 
                  Palpha2(ist+i) = dreal(zf(i)) - asum 
               end if
ccc               Palpha2(ist+i) = Palpha2(ist+i) + ck(kbod)
            end do
            ist = ist+nd
         end do
c
c  Construct E o alpha
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               if (kbod.eq.0) then
                  Ealpha1(ist+i) = 2.d0*alpha1(ist+i) 
     1                              + 4.d0*rkappa(ist+i)*Palpha2(ist+i)
                 else
                  Ealpha1(ist+i) = 2.d0*alpha1(ist+i) 
ccc     1                              - 4.d0*rkappa(ist+i)*Palpha2(ist+i)
     1                              + 4.d0*rkappa(ist+i)*Palpha2(ist+i)
               end if
               Ealpha2(ist+i) = 2.d0*Palpha2(ist+i)
            end do
            ist = ist + nd
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine OPERATE_E_CHK (k0, k, nd, nbk, alpha1, alpha2, dsdth,  
     1                      rkappa, Ealpha1, Ealpha2, Palpha2, zf, 
     2                      wsave)
c---------------
c Calculates E o alpha, where
c     E = [2 I  4 kappa P
c            0        2 P ]
c and
c     P o alpha = \int^s alpha(s) ds
c
      implicit real*8 (a-h,o-z)
      dimension alpha1(nbk), alpha2(nbk), dsdth(nbk), rkappa(nbk)
      dimension Palpha2(nbk), Ealpha1(nbk), Ealpha2(nbk), ck(0:100)
      complex*16 zf(nbk), wsave(*), eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
         n = nd/2
c
         call DCFFTI (nd, wsave)
c
         ist = nd
         ck(0) = 0.d0
         do kbod = 1, k
            do i = 1, nd
               zf(i) = alpha1(ist+i)*dsdth(ist+i)
            end do
            call DCFFTF (nd, zf, wsave)
            ck(kbod) = 2.d0*pi*dreal(zf(1))/nd
            ist = ist + nd
         end do
c
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               zf(i) = alpha2(ist+i)*dsdth(ist+i)
            end do
            call DCFFTF (nd, zf, wsave)
            do i = 1, nd
               zf(i) = zf(i)/nd
            end do
            f0 = dreal(zf(1))
ccc            call PRIN2 (' in OPERATE_E, f0 = *', f0, 1)
            zf(1) = 0.d0
            do kmode = 1, n-1
               zf(kmode+1) = zf(kmode+1)/(eye*kmode)
               zf(nd-kmode+1) = -zf(nd-kmode+1)/(eye*kmode)
            end do
            zf(n+1) = 0
            call DCFFTB (nd, zf, wsave)
            asum = dreal(zf(1))
ccc            ftot = f0*2.d0*pi
            do i = 1, nd
               th = 2.d0*pi*(i-1.d0)/nd
               if (kbod.eq.0) then
                  Palpha2(ist+i) = dreal(zf(i)) - asum 
                 else
ccc                  Palpha2(ist+i) = -dreal(zf(i)) + asum 
                  Palpha2(ist+i) = dreal(zf(i)) - asum 
               end if
ccc               Palpha2(ist+i) = Palpha2(ist+i) + ck(kbod)
            end do
            ist = ist+nd
         end do
c
c  Construct E o alpha
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               if (kbod.eq.0) then
                  Ealpha1(ist+i) = 2.d0*alpha1(ist+i) 
     1                              + 4.d0*rkappa(ist+i)*Palpha2(ist+i)
                 else
                  Ealpha1(ist+i) = 2.d0*alpha1(ist+i) 
ccc     1                              - 4.d0*rkappa(ist+i)*Palpha2(ist+i)
     1                              + 4.d0*rkappa(ist+i)*Palpha2(ist+i)
               end if
               Ealpha2(ist+i) = 2.d0*Palpha2(ist+i)
            end do
            ist = ist + nd
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine DIAG (k0, k, nd, nbk, Ealpha1, Ealpha2, dsdth,  
     1                      rkappa, alpha1, alpha2, Palpha2, zf, 
     2                      wsave)
c---------------
c Calculates E o alpha, where
c     E = [2 I  4 kappa P
c            0        2 P ]
c and
c     P o alpha = \int^s alpha(s) ds
c
      implicit real*8 (a-h,o-z)
      dimension alpha1(nbk), alpha2(nbk), dsdth(nbk), rkappa(nbk)
      dimension Palpha2(nbk), Ealpha1(nbk), Ealpha2(nbk), ck(0:100)
      complex*16 zf(nbk), wsave(*), eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
         n = nd/2
c
         call DCFFTI (nd, wsave)
c
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               zf(i) = Ealpha2(ist+i)
            end do
            call DCFFTF (nd, zf, wsave)
            do i = 1, nd
               zf(i) = zf(i)/nd
            end do
            zf(1) = 0.d0
            do kmode = 1, n-1
               zf(kmode+1) = zf(kmode+1)*(eye*kmode)
               zf(nd-kmode+1) = -zf(nd-kmode+1)*(eye*kmode)
            end do
            zf(n+1) = 0
            call DCFFTB (nd, zf, wsave)
            do i = 1, nd
               if (kbod.eq.0) then
                  Palpha2(ist+i) = dreal(zf(i))/dsdth(ist+i) 
                 else
ccc                  Palpha2(ist+i) = -dreal(zf(i))/dsdth(ist+i) 
                  Palpha2(ist+i) = dreal(zf(i))/dsdth(ist+i) 
               end if
            end do
            ist = ist+nd
         end do
c
c  Construct E o alpha
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               if (kbod.eq.0) then
                  alpha1(ist+i) = 0.5d0*Ealpha1(ist+i) 
     1                              - rkappa(ist+i)*Ealpha2(ist+i)
                 else
                  alpha1(ist+i) = 0.5d0*Ealpha1(ist+i) 
     1                              - rkappa(ist+i)*Ealpha2(ist+i)
ccc     1                              + rkappa(ist+i)*Ealpha2(ist+i)
               end if
               alpha2(ist+i) = 0.5d0*Palpha2(ist+i)
            end do
            ist = ist + nd
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine EVAL_G1 (k0, k, nd, nbk, alambda, z, zn, dsdth,
     1                    sigma, z_tar, G1, zf, wsave)
c---------------
c  Evaluate \int sigma G1(y,z_tar) ds
c  it's assumed z_tar is not on the boundary
c
      implicit real*8 (a-h,o-z)
      dimension dsdth(nbk), rkappa(nbk), sigma(nbk), bk(0:1)
      complex*16 z(nbk), zn(nbk), z_tar, zdis, eye, zdir, 
     1           zf(nbk), wsave(*), zn_tar, zs_tar, zs_src
c
      pi=4.0d0*datan(1.0d0)
      eye = dcmplx(0.d0,1.d0)
c
c  mesh spacing
         dth = 2.d0*pi/nd
c
c Integrate
         G1 = 0.0d0
         rmin = 1.d10
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               zdis = z_tar - z(ist+i)
               r = cdabs(zdis)
               rmin = min(r,rmin)
               vn = dreal(zdis*dconjg(zn(ist+i)))
               az = alambda*r
               if (az.lt.100.d0) then
                  call RKBESL (alambda*r, 0.d0, 2, 1, bk, ncalc)
                 else
                  call RKBESL (alambda*r, 0.d0, 2, 2, bk, ncalc)
                  do ik = 0, 1
                     bk(ik) = dexp(-az)*bk(ik)
                  end do
               end if
                  if (ncalc.ne.2) then 
                     write (6,*) 'ERROR IN EVAL_G1! ', ncalc
                     stop
                  end if
               val1 = dlog(r)/(2.d0*pi) + (vn**2)*bk(0)/(pi*r**2)
               val2 = (2.d0*vn**2 - r**2)*(alambda*r*bk(1) - 1.d0)
     1                     /(pi*(alambda**2)*r**4)
               val = val1 + val2
               zf(i) = val*sigma(ist+i)*dsdth(ist+i)
            end do
            call DCFFTF (nd, zf, wsave)
            G1 = G1 + 2.d0*pi*dreal(zf(1))/nd
            ist = ist + nd
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine KERNEL_IE (alambda, z_tar, z_src, zn_tar, zn_src, 
     1                      zs_tar, zs_src, c0, c1, c2, c3, c4, c5, 
     2                      akern)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension akern(2,2)
      complex*16 z_tar, z_src, zn_tar, zn_src, zs_tar, zs_src
      complex*16 zdis, eye
      data euler /0.57721566490153286d0/
      common /cutoff/ rcut
c
      pi = 4.0d0*datan(1.0d0)
      eye = dcmplx(0.d0,1.d0)
c
            zdis = z_tar - z_src
            r = cdabs(zdis)
            rnu = dreal(zdis*dconjg(zn_src))
            rnux = dreal(zdis*dconjg(zn_tar))
            anunux = dreal(zn_src*dconjg(zn_tar))
            anutaux = dreal(zn_src*dconjg(zs_tar))
            rtaux = dreal(zdis*dconjg(zs_tar))
            rtau = dreal(zdis*dconjg(zs_src))
            az = alambda*r
            if (az.gt.rcut) then
               cc0 = c0
               cc1 = c1
               cc2 = c2
               cc3 = c3
               cc4 = c4
               cc5 = c5
             else
               alnz=dlog(2.d0/az)-euler
               az2=az**2
               az4=az2**2
               cc0=-0.5d0 + (alnz+0.75d0)*az2/8.d0 
     1             + (alnz+17.d0/12.d0)*az4/96.d0
               cc1=-(alnz-0.25d0)*az2/8.d0-(alnz+13.d0/12.d0)*az4/32.d0
               cc2=-1.d0+az2/8.d0-(alnz+11.d0/12.d0)*az4/48.d0
               cc3=-az2/8.d0+(alnz+7.d0/12.d0)*az4/16.d0
               cc4=-4.d0+az2/4.d0-az4/48.d0
               cc5=-(alnz+0.25d0)*az2/4.d0-(alnz+7.d0/6.d0)*az4/24.d0
            end if
            g11 = 2.d0*rnu*anunux*cc0 + 0.5d0*rnux*cc5
     1                - rnux*rnu**2*cc2/r**2
            g21 = 2.d0*rnu*anutaux*cc0 + 0.5d0*rtaux*cc5
     1                - rtaux*rnu**2*cc2/r**2
            g11 = g11/(pi*r**2)
            g21 = g21/(pi*r**2)
            g12 = -anunux*cc1 + 3.d0*anunux*rnu**2*cc2/r**2
     1               +rnux*rnu*cc3/r**2 - rnux*rnu**3*cc4/r**4
            g22 = -anutaux*cc1 + 3.d0*anutaux*rnu**2*cc2/r**2
     1               +rtaux*rnu*cc3/r**2 - rtaux*rnu**3*cc4/r**4
            g12 = g12/(pi*r**2)
            g22 = g22/(pi*r**2)
            akern(1,1) = g11
            akern(1,2) = g12
            akern(2,1) = g21
            akern(2,2) = g22
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine GRAD_KERNEL (alambda, z_tar, z_src, zn_src, grad_kern)
c---------------
c
c grad_G1 = (grad_kern(1,1) , grad_kern(1,2))
c grad_G2 = (grad_kern(2,1) , grad_kern(2,2))
      implicit real*8 (a-h,o-z)
      dimension grad_kern(2,2), bk(0:1)
      complex*16 z_tar, z_src, zn_src
      complex*16 zdis, eye
      data euler /0.57721566490153286d0/
      common /cutoff/ rcut
c
      pi = 4.0d0*datan(1.0d0)
      eye = dcmplx(0.d0,1.d0)
c
            zdis = z_tar - z_src
            r = cdabs(zdis)
c
            rx = dreal(zdis)
            ry = dimag(zdis)
            anux = dreal(zn_src)
            anuy = dimag(zn_src)
            rnu = dreal(zdis*dconjg(zn_src))
c
            az = alambda*r
            if (az.gt.rcut) then
               if (az.lt.100.d0) then
                  call RKBESL (alambda*r, 0.d0, 2, 1, bk, ncalc)
                 else
                  call RKBESL (alambda*r, 0.d0, 2, 2, bk, ncalc)
                  do ik = 0, 1
                     bk(ik) = dexp(-az)*bk(ik)
                  end do
               end if
               if (ncalc.ne.2) then 
                  write (6,*) 'ERROR IN GRAD_KERNEL! ', ncalc
                  stop
               end if
               cc0 = bk(0) + 2.d0*bk(1)/az - 2.d0/az**2
               cc1 = 3.d0*cc0 + az*bk(1) + 0.5d0
               cc2 = 4.d0*cc0 + az*bk(1)
               cc3 = 12.d0*cc0 + 5.d0*az*bk(1) + az**2*bk(0)
     1                   + 1.d0
               cc4 = 24.d0*cc0 + 8.d0*az*bk(1) + az**2*bk(0)
               cc5 = 2.d0*cc0 + az*bk(1)
             else
               alnz=dlog(2.d0/az)-euler
               az2=az**2
               az4=az2**2
               cc0=-0.5d0 + (alnz+0.75d0)*az2/8.d0 
     1             + (alnz+17.d0/12.d0)*az4/96.d0
               cc1=-(alnz-0.25d0)*az2/8.d0-(alnz+13.d0/12.d0)*az4/32.d0
               cc2=-1.d0+az2/8.d0-(alnz+11.d0/12.d0)*az4/48.d0
               cc3=-az2/8.d0+(alnz+7.d0/12.d0)*az4/16.d0
               cc4=-4.d0+az2/4.d0-az4/48.d0
               cc5=-(alnz+0.25d0)*az2/4.d0-(alnz+7.d0/6.d0)*az4/24.d0
            end if
            g11 = 2.d0*rnu*anux*cc0 + 0.5d0*rx*cc5
     1                - rx*rnu**2*cc2/r**2
            g21 = 2.d0*rnu*anuy*cc0 + 0.5d0*ry*cc5
     1                - ry*rnu**2*cc2/r**2
            g11 = g11/(pi*r**2)
            g21 = g21/(pi*r**2)
            g12 = -anux*cc1 + 3.d0*anux*rnu**2*cc2/r**2
     1               +rx*rnu*cc3/r**2 - rx*rnu**3*cc4/r**4
            g22 = -anuy*cc1 + 3.d0*anuy*rnu**2*cc2/r**2
     1               +ry*rnu*cc3/r**2 - ry*rnu**3*cc4/r**4
            g12 = g12/(pi*r**2)
            g22 = g22/(pi*r**2)
            grad_kern(1,1) = g11
            grad_kern(1,2) = g12
            grad_kern(2,1) = g21
            grad_kern(2,2) = g22
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine KERNEL_PSI (alambda, z_tar, z_src, zn_src, G)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension G(2), bk(0:1)
      complex*16 z_tar, z_src, zn_src
      complex*16 zdis, eye
      data euler /0.57721566490153286d0/
      common /cutoff/ rcut
c
         pi = 4.0d0*datan(1.0d0)
         eye = dcmplx(0.d0,1.d0)
c
         zdis = z_tar - z_src
         r = cdabs(zdis)
         rnu = dreal(zdis*dconjg(zn_src))
         az = alambda*r
         if (az.gt.rcut) then
               if (az.lt.100.d0) then
                  call RKBESL (alambda*r, 0.d0, 2, 1, bk, ncalc)
                 else
                  call RKBESL (alambda*r, 0.d0, 2, 2, bk, ncalc)
                  do ik = 0, 1
                     bk(ik) = dexp(-az)*bk(ik)
                  end do
               end if
               if (ncalc.ne.2) then 
                  write (6,*) 'ERROR IN KERNEL_PSI! ', ncalc
                  stop
               end if
               cc0 = bk(0) + 2.d0*bk(1)/az - 2.d0/az**2
               cc1 = 3.d0*cc0 + az*bk(1) + 0.5d0
               cc2 = 4.d0*cc0 + az*bk(1)
               cc3 = 12.d0*cc0 + 5.d0*az*bk(1) + az**2*bk(0)
     1                   + 1.d0
               cc4 = 24.d0*cc0 + 8.d0*az*bk(1) + az**2*bk(0)
               cc5 = 2.d0*cc0 + az*bk(1)
          else
            alnz=dlog(2.d0/az)-euler
            az2=az**2
            az4=az2**2
            cc0=-0.5d0 + (alnz+0.75d0)*az2/8.d0 
     1             + (alnz+17.d0/12.d0)*az4/96.d0
            cc1=-(alnz-0.25d0)*az2/8.d0-(alnz+13.d0/12.d0)*az4/32.d0
            cc2=-1.d0+az2/8.d0-(alnz+11.d0/12.d0)*az4/48.d0
            cc3=-az2/8.d0+(alnz+7.d0/12.d0)*az4/16.d0
            cc4=-4.d0+az2/4.d0-az4/48.d0
            cc5=-(alnz+0.25d0)*az2/4.d0-(alnz+7.d0/6.d0)*az4/24.d0
         end if
         g1 = -(0.5d0-(rnu/r)**2)*cc0/pi
         g2 = -rnu*cc1 + rnu**3*cc2/r**2
         g2 = g2/(pi*r**2)
c
         G(1) = g1
         G(2) = g2
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine ALPERT_MAT (k0, k, nd, nbk, alambda, z, z_vec1,  
     1               z_vec2, dsdth, ds_vec1, ds_vec2, zn, zn_vec1, 
     2               zn_vec2, zs, zs_vec1, zs_vec2, f_vec1, f_vec2,
     3               c0, c1, c2, c3, c4, c5, G11, G12, G21, 
     6               G22, ds_k)
c---------------
      implicit real*8 (a-h,o-z)
      complex*16 z(nbk), z_vec1(nbk,15), z_vec2(nbk,15)
      dimension dsdth(nbk), ds_vec1(nbk,15), ds_vec2(nbk,15) 
      complex*16 zn(nbk), zn_vec1(nbk,15), zn_vec2(nbk,15)
      complex*16 zs(nbk), zs_vec1(nbk,15), zs_vec2(nbk,15)
      dimension f_vec1(nd,nd,15), f_vec2(nd,nd,15), ds_k(k0:k)
      dimension c0(nbk,nbk), c1(nbk,nbk), c2(nbk,nbk), c3(nbk,nbk),
     1          c4(nbk,nbk), c5(nbk,nbk), G11(nbk,nbk), G12(nbk,nbk), 
     1          G21(nbk,nbk), G22(nbk,nbk)
      complex*16 zdis, eye
      complex*16 z_tar, z_src, zn_tar, zn_src, zs_tar, zs_src
      dimension u(16), v(16), akern(2,2)
      common /alpert/ u, v, numpts, ishift, iorder, iquad
c
         pi = 4.0d0*datan(1.0d0)
         eye = dcmplx(0.d0,1.d0)
c
c  mesh spacing
         dth = 2.d0*pi/nd
c
c zero everything
         do i = 1, nbk
            do j = 1, nbk
               G11(i,j) = 0.d0
               G12(i,j) = 0.d0
               G21(i,j) = 0.d0
               G22(i,j) = 0.d0
            end do
         end do
c
c  Trapezoid rule
         ist = 0
         rmin = 1.d10
         rmax = 0.d0
         do kbod = k0, k
            do i = 1, nd
               jst = 0
               do lbod = k0, k
                  ds_fact = ds_k(kbod)/ds_k(lbod)
                  do j = 1, nd
                     if (ist+i.ne.jst+j) then
                        z_tar = z(ist+i)
                        z_src = z(jst+j)
                        r = cdabs(z_tar - z_src)
                        rmax = max(r,rmax)
                        rmin = min(r,rmin)
                        zn_tar = zn(ist+i)
                        zn_src = zn(jst+j)
                        zs_tar = zs(ist+i)
                        zs_src = zs(jst+j)
                        cc0 = c0(ist+i,jst+j)
                        cc1 = c1(ist+i,jst+j)
                        cc2 = c2(ist+i,jst+j)
                        cc3 = c3(ist+i,jst+j)
                        cc4 = c4(ist+i,jst+j)
                        cc5 = c5(ist+i,jst+j)
                        call KERNEL_IE (alambda, z_tar, z_src, zn_tar,  
     1                           zn_src, zs_tar, zs_src, cc0, cc1, 
     2                           cc2, cc3, cc4, cc5, akern)
                        val = akern(1,1)*dsdth(jst+j)*ds_fact
                        G11(ist+i,jst+j) = val*dth
                        val = akern(1,2)*dsdth(jst+j)
                        G12(ist+i,jst+j) = val*dth
                        val = akern(2,1)*dsdth(jst+j)
                        G21(ist+i,jst+j) = val*dth
                        val = akern(2,2)*dsdth(jst+j)
                        G22(ist+i,jst+j) = val*dth
                     end if
                  end do
                  jst = jst + nd
               end do
            end do
            ist = ist + nd
         end do
c
c     this part of the routine will be edventually be done inside of the
c     FMM.  However, this way I don't have to change the method of doing
c     the quadrature correction where the nodes near the singularity are
c     subtracted
c
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               th = (i-1.d0)*dth
               do j = 1, nd
                  do inode = 1, numpts
c
c  1st term in Alpert's paper
                     z_src = z_vec1(ist+i,inode)
                     z_tar = z(ist+i)
                     zn_src = zn_vec1(ist+i,inode)
                     zn_tar = zn(ist+i)
                     zs_src = zs_vec1(ist+i,inode)
                     zs_tar = zs(ist+i)
                     r = cdabs(z_tar - z_src)
                        rmax = max(r,rmax)
                        rmin = min(r,rmin)
                     az = alambda*r
                     call GET_BESSEL_PNT (az, cc0, cc1, cc2, cc3, cc4, 
     1                                 cc5)
                     call KERNEL_IE (alambda, z_tar, z_src, zn_tar,  
     1                            zn_src, zs_tar, zs_src, cc0, cc1, 
     2                            cc2, cc3, cc4, cc5, akern)
                     ds_interp = ds_vec1(ist+i,inode)
                     val = akern(1,1)*ds_interp
                     G11(ist+i,ist+j) = G11(ist+i,ist+j) 
     1                  + u(inode)*f_vec1(i,j,inode)*val*dth
                     val = akern(1,2)*ds_interp
                     G12(ist+i,ist+j) = G12(ist+i,ist+j) 
     1                  + u(inode)*f_vec1(i,j,inode)*val*dth
                     val = akern(2,1)*ds_interp
                     G21(ist+i,ist+j) = G21(ist+i,ist+j) 
     1                  + u(inode)*f_vec1(i,j,inode)*val*dth
                     val = akern(2,2)*ds_interp
                     G22(ist+i,ist+j) = G22(ist+i,ist+j) 
     1                  + u(inode)*f_vec1(i,j,inode)*val*dth
               
c
c  3rd term in Alpert's paper
                     z_src = z_vec2(ist+i,inode)
                     z_tar = z(ist+i)
                     zn_src = zn_vec2(ist+i,inode)
                     zn_tar = zn(ist+i)
                     zs_src = zs_vec2(ist+i,inode)
                     zs_tar = zs(ist+i)
                     r = cdabs(z_tar - z_src)
                         rmax = max(r,rmax)
                         rmin = min(r,rmin)
                     az = alambda*r
                     call GET_BESSEL_PNT (az, cc0, cc1, cc2, cc3, cc4, 
     1                                 cc5)
                     call KERNEL_IE (alambda, z_tar, z_src, zn_tar,  
     1                            zn_src, zs_tar, zs_src, cc0, cc1, 
     2                            cc2, cc3, cc4, cc5, akern)
                     ds_interp = ds_vec2(ist+i,inode)
                     val = akern(1,1)*ds_interp
                     G11(ist+i,ist+j) = G11(ist+i,ist+j) 
     1                  + u(inode)*f_vec2(i,j,inode)*val*dth
                     val = akern(1,2)*ds_interp
                     G12(ist+i,ist+j) = G12(ist+i,ist+j) 
     1                  + u(inode)*f_vec2(i,j,inode)*val*dth
                     val = akern(2,1)*ds_interp
                     G21(ist+i,ist+j) = G21(ist+i,ist+j) 
     1                  + u(inode)*f_vec2(i,j,inode)*val*dth
                     val = akern(2,2)*ds_interp
                     G22(ist+i,ist+j) = G22(ist+i,ist+j) 
     1                  + u(inode)*f_vec2(i,j,inode)*val*dth
                  end do
               end do
c
c  2nd term in Alpert's paper - Subtract off the
c     terms that should not be included in the 'trapezoidal rule'
c     part of the quadrature
               mid_sum_end = nd-2*ishift+1
               do inode = mid_sum_end,nd-1
                  index1=mod(i+ishift+inode,nd)
                  if (index1 .eq. 0) index1 = nd
                  if (i .ne. index1) then
                     z_src = z(ist+index1)
                     z_tar = z(ist+i)
                     r = cdabs(z_tar-z_src)
                        rmax = max(r,rmax)
                        rmin = min(r,rmin)
                     zn_src = zn(ist+index1)
                     zn_tar = zn(ist+i)
                     zs_src = zs(ist+index1)
                     zs_tar = zs(ist+i)
                     cc0 = c0(ist+index1,ist+i)
                     cc1 = c1(ist+index1,ist+i)
                     cc2 = c2(ist+index1,ist+i)
                     cc3 = c3(ist+index1,ist+i)
                     cc4 = c4(ist+index1,ist+i)
                     cc5 = c5(ist+index1,ist+i)
                     call KERNEL_IE (alambda, z_tar, z_src, zn_tar,  
     1                           zn_src, zs_tar, zs_src, cc0, cc1, 
     2                           cc2, cc3, cc4, cc5, akern)
                     val11 = akern(1,1)*dsdth(ist+index1)
                     val12 = akern(1,2)*dsdth(ist+index1)
                     val21 = akern(2,1)*dsdth(ist+index1)
                     val22 = akern(2,2)*dsdth(ist+index1)
                    else
                     val11 = 0.d0
                     val12 = 0.d0
                     val21 = 0.d0
                     val22 = 0.d0
                  endif
                  G11(ist+i,ist+index1) = G11(ist+i,ist+index1) 
     1                              - val11*dth
                  G12(ist+i,ist+index1) = G12(ist+i,ist+index1) 
     1                               - val12*dth
                  G21(ist+i,ist+index1) = G21(ist+i,ist+index1) 
     1                               - val21*dth
                  G22(ist+i,ist+index1) = G22(ist+i,ist+index1) 
     1                               - val22*dth
               end do
            end do
            ist = ist + nd
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine ALPERT_QUAD_GEN (k0, k, nd, nbk, alambda, z, z_vec1,  
     1               z_vec2, dsdth, ds_vec1, ds_vec2, zn, zn_vec1, 
     2               zn_vec2, zs, zs_vec1, zs_vec2, Ealpha1,   
     4               Ealph1_vec1, Ealph1_vec2, Ealpha2, Ealph2_vec1, 
     5               Ealph2_vec2, c0, c1, c2, c3, c4, c5, G11, G12, G21, 
     6               G22, ds_k)
c---------------
      implicit real*8 (a-h,o-z)
      complex*16 z(nbk), z_vec1(nbk,15), z_vec2(nbk,15)
      dimension dsdth(nbk), ds_vec1(nbk,15), ds_vec2(nbk,15), ds_k(k0:k)
      complex*16 zn(nbk), zn_vec1(nbk,15), zn_vec2(nbk,15)
      complex*16 zs(nbk), zs_vec1(nbk,15), zs_vec2(nbk,15)
      dimension Ealpha1(nbk), Ealph1_vec1(nbk,15), Ealph1_vec2(nbk,15),  
     4          Ealpha2(nbk), Ealph2_vec1(nbk,15), Ealph2_vec2(nbk,15)
      dimension c0(nbk,nbk), c1(nbk,nbk), c2(nbk,nbk), c3(nbk,nbk),
     1          c4(nbk,nbk), c5(nbk,nbk)
      complex*16 zdis, eye
      complex*16 z_tar, z_src, zn_tar, zn_src, zs_tar, zs_src
      dimension u(16), v(16), G11(nbk), G12(nbk), G21(nbk), G22(nbk),
     1          akern(2,2)
      common /alpert/ u, v, numpts, ishift, iorder, iquad
c
         pi = 4.0d0*datan(1.0d0)
         eye = dcmplx(0.d0,1.d0)
c
c  mesh spacing
         dth = 2.d0*pi/nd
c
c  Trapezoid rule
         ist = 0
         rmin = 1.d10
         rmax = 0.d0
         do kbod = k0, k
            do i = 1, nd
               jst = 0
               G11(ist+i) = 0.d0
               G12(ist+i) = 0.d0
               G21(ist+i) = 0.d0
               G22(ist+i) = 0.d0
               do lbod = k0, k
                  ds_fact = ds_k(kbod)/ds_k(lbod)
                  do j = 1, nd
                     if (ist+i.ne.jst+j) then
                        z_tar = z(ist+i)
                        z_src = z(jst+j)
                        r = cdabs(z_tar - z_src)
                        rmax = max(r,rmax)
                        rmin = min(r,rmin)
                        zn_tar = zn(ist+i)
                        zn_src = zn(jst+j)
                        zs_tar = zs(ist+i)
                        zs_src = zs(jst+j)
                        cc0 = c0(ist+i,jst+j)
                        cc1 = c1(ist+i,jst+j)
                        cc2 = c2(ist+i,jst+j)
                        cc3 = c3(ist+i,jst+j)
                        cc4 = c4(ist+i,jst+j)
                        cc5 = c5(ist+i,jst+j)
                        call KERNEL_IE (alambda, z_tar, z_src, zn_tar,  
     1                           zn_src, zs_tar, zs_src, cc0, cc1, 
     2                           cc2, cc3, cc4, cc5, akern)
                        val = akern(1,1)*dsdth(jst+j)*ds_fact
                        G11(ist+i) = G11(ist+i) 
     1                               + Ealpha1(jst+j)*val*dth
                        val = akern(1,2)*dsdth(jst+j)
                        G12(ist+i) = G12(ist+i) 
     1                               + Ealpha2(jst+j)*val*dth
                        val = akern(2,1)*dsdth(jst+j)
                        G21(ist+i) = G21(ist+i) 
     1                               + Ealpha1(jst+j)*val*dth
                        val = akern(2,2)*dsdth(jst+j)
                        G22(ist+i) = G22(ist+i) 
     1                               + Ealpha2(jst+j)*val*dth
                     end if
                  end do
                  jst = jst + nd
               end do
            end do
            ist = ist + nd
         end do
c
c     this part of the routine will be edventually be done inside of the
c     FMM.  However, this way I don't have to change the method of doing
c     the quadrature correction where the nodes near the singularity are
c     subtracted
c
         ist = 0
         do kbod = k0, k
            do i = 1, nd
               th = (i-1.d0)*dth
               do inode = 1, numpts
c
c  1st term in Alpert's paper
                  z_src = z_vec1(ist+i,inode)
                  z_tar = z(ist+i)
                  zn_src = zn_vec1(ist+i,inode)
                  zn_tar = zn(ist+i)
                  zs_src = zs_vec1(ist+i,inode)
                  zs_tar = zs(ist+i)
                  r = cdabs(z_tar - z_src)
                        rmax = max(r,rmax)
                        rmin = min(r,rmin)
                  az = alambda*r
                  call GET_BESSEL_PNT (az, cc0, cc1, cc2, cc3, cc4, 
     1                                 cc5)
                  call KERNEL_IE (alambda, z_tar, z_src, zn_tar,  
     1                            zn_src, zs_tar, zs_src, cc0, cc1, 
     2                            cc2, cc3, cc4, cc5, akern)
                  ds_interp = ds_vec1(ist+i,inode)
                  Ealph1_int = Ealph1_vec1(ist+i,inode)
                  Ealph2_int = Ealph2_vec1(ist+i,inode)
                  val = akern(1,1)*ds_interp
                  G11(ist+i) = G11(ist+i) 
     1                            + u(inode)*Ealph1_int*val*dth
                  val = akern(1,2)*ds_interp
                  G12(ist+i) = G12(ist+i) 
     1                            + u(inode)*Ealph2_int*val*dth
                  val = akern(2,1)*ds_interp
                  G21(ist+i) = G21(ist+i) 
     1                            + u(inode)*Ealph1_int*val*dth
                  val = akern(2,2)*ds_interp
                  G22(ist+i) = G22(ist+i) 
     1                            + u(inode)*Ealph2_int*val*dth
c
c  3rd term in Alpert's paper
                  z_src = z_vec2(ist+i,inode)
                  z_tar = z(ist+i)
                  zn_src = zn_vec2(ist+i,inode)
                  zn_tar = zn(ist+i)
                  zs_src = zs_vec2(ist+i,inode)
                  zs_tar = zs(ist+i)
                  r = cdabs(z_tar - z_src)
                         rmax = max(r,rmax)
                        rmin = min(r,rmin)
                  az = alambda*r
                  call GET_BESSEL_PNT (az, cc0, cc1, cc2, cc3, cc4, 
     1                                 cc5)
                  call KERNEL_IE (alambda, z_tar, z_src, zn_tar,  
     1                            zn_src, zs_tar, zs_src, cc0, cc1, 
     2                            cc2, cc3, cc4, cc5, akern)
                  ds_interp = ds_vec2(ist+i,inode)
                  Ealph1_int = Ealph1_vec2(ist+i,inode)
                  Ealph2_int = Ealph2_vec2(ist+i,inode)
                  val = akern(1,1)*ds_interp
                  G11(ist+i) = G11(ist+i) 
     1                            + u(inode)*Ealph1_int*val*dth
                  val = akern(1,2)*ds_interp
                  G12(ist+i) = G12(ist+i) 
     1                            + u(inode)*Ealph2_int*val*dth
                  val = akern(2,1)*ds_interp
                  G21(ist+i) = G21(ist+i) 
     1                            + u(inode)*Ealph1_int*val*dth
                  val = akern(2,2)*ds_interp
                  G22(ist+i) = G22(ist+i) 
     1                            + u(inode)*Ealph2_int*val*dth
               end do
c
c  2nd term in Alpert's paper - Subtract off the
c     terms that should not be included in the 'trapezoidal rule'
c     part of the quadrature
               mid_sum_end = nd-2*ishift+1
               do inode = mid_sum_end,nd-1
                  index1=mod(i+ishift+inode,nd)
                  if (index1 .eq. 0) index1 = nd
                  if (i .ne. index1) then
                     z_src = z(ist+index1)
                     z_tar = z(ist+i)
                     r = cdabs(z_tar-z_src)
                        rmax = max(r,rmax)
                        rmin = min(r,rmin)
                     zn_src = zn(ist+index1)
                     zn_tar = zn(ist+i)
                     zs_src = zs(ist+index1)
                     zs_tar = zs(ist+i)
                     cc0 = c0(ist+index1,ist+i)
                     cc1 = c1(ist+index1,ist+i)
                     cc2 = c2(ist+index1,ist+i)
                     cc3 = c3(ist+index1,ist+i)
                     cc4 = c4(ist+index1,ist+i)
                     cc5 = c5(ist+index1,ist+i)
                     call KERNEL_IE (alambda, z_tar, z_src, zn_tar,  
     1                           zn_src, zs_tar, zs_src, cc0, cc1, 
     2                           cc2, cc3, cc4, cc5, akern)
                     val11 = akern(1,1)*dsdth(ist+index1)
                     val12 = akern(1,2)*dsdth(ist+index1)
                     val21 = akern(2,1)*dsdth(ist+index1)
                     val22 = akern(2,2)*dsdth(ist+index1)
                    else
                     val11 = 0.d0
                     val12 = 0.d0
                     val21 = 0.d0
                     val22 = 0.d0
                  endif
                  G11(ist+i) = G11(ist+i) 
     1                              - Ealpha1(ist+index1)*val11*dth
                  G12(ist+i) = G12(ist+i) 
     1                               - Ealpha2(ist+index1)*val12*dth
                  G21(ist+i) = G21(ist+i) 
     1                               - Ealpha1(ist+index1)*val21*dth
                  G22(ist+i) = G22(ist+i) 
     1                               - Ealpha2(ist+index1)*val22*dth
               end do
            end do
            ist = ist + nd
         end do
c
      return
      end
c
c*********************
c
      subroutine DUMP (nx,ny,ugrid,igrid,ireal,if)
      implicit real*8 (a-h,o-z)
      dimension ugrid(nx,ny), igrid(nx,ny)
c
         DO i = 1,NX
            do j = 1, NY
               if (ireal.eq.1) then 
                  write(if,'(e20.13,$)')(ugrid(I,J))
                  write (if,'(a)')  ''
                 else
                  write(if,'(i4,$)') (igrid(i,j))
                  write (if,'(a)')  ''
               end if
            end do
         ENDDO
c
      return
      end
