C
C
        SUBROUTINE RSPLOT(X,Y,N,NBOD,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION X(1),Y(1),V(4)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        X,Y (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      WRITE(IW,*) 'x = ['
      WRITE(IW,1400) (X(I),Y(I),I=1,N)
ccc      write (iw,1400) x(1),y(1)
      WRITE(IW,*) '];'
1300  FORMAT(D15.6)
1400  FORMAT(2D15.6)
      write(iw,*) ' xx = x(:,1);'
      write(iw,*) ' yy = x(:,2);'
      write(iw, *) ' plot(xx,yy)'
ccc      write(iw, *) ' plot(xx,yy,\'o\')'
      write(iw,*) 'axis equal'
ccc      write(iw, *) ' plot(x)'
      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
C
      entry rsplt1(x,n,iw)
      WRITE(IW,*) 'x = ['
      WRITE(IW,1300) (X(I),I=1,N)
      WRITE(IW,*) '];'
      RETURN
C
C
      ENTRY RSPINI(IW,XMIN,XMAX,YMIN,YMAX)
      II = 0
      V(1) = XMIN
      V(2) = XMAX
      V(3) = YMIN
      V(4) = YMAX
      write(IW,*) 'v = ['
      write(iw,1300) (v(i),i=1,4)
      write(iw,*) '];'
      write(iw,*) ' axis(v)'
ccc      write(iw,*) ' axis(\'square\')'
          RETURN
          END
C
C
        SUBROUTINE RSLOGPLOT(X,Y,N,NBOD,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION X(1),Y(1),V(4)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        X,Y (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      WRITE(IW,*) 'x = ['
      WRITE(IW,1400) (X(I),Y(I),I=1,N)
ccc      write (iw,1400) x(1),y(1)
      WRITE(IW,*) '];'
1300  FORMAT(D15.6)
1400  FORMAT(2D15.6)
      write(iw,*) ' xx = x(:,1);'
      write(iw,*) ' yy = x(:,2);'
      write(iw, *) ' semilogy(xx,yy)'
ccc      write(iw, *) ' plot(xx,yy,\'o\')'
ccc      write(iw,*) 'axis equal'
ccc      write(iw, *) ' plot(x)'
      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
      end
C
C
        SUBROUTINE RSCLINEPLOT(Z,C,time,N,NBOD,k,iframe,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
        complex*16 Z(*)
        DIMENSION V(4), C(*)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        X,Y (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      ymax = -1.d6
      ymin = 1.d6
      WRITE(IW,*) 'x = ['
      do i = 1, N
         WRITE(IW,1400) Z(I), C(I)
         ymin = min(ymin,dimag(z(i)))
         ymax = max(ymax,dimag(z(i)))
      end do
      write(IW,1400) (Z(1), C(1))
ccc      write (iw,1400) x(1),y(1)
      WRITE(IW,*) '];'
      if (nbod.eq.k) then
         write (iw,1200) time
      write (iw,*) "title([' t= ',num2str(time)]);"
      end if
1200  FORMAT ('time=',F6.2)
1300  FORMAT(D15.6)
1400  FORMAT(3D15.6)
      if (NBOD.eq.1) write (iw,*) 'figure(1);clf;'
      write(iw,*) ' xx = transpose(x(:,1));'
      write(iw,*) ' yy = transpose(x(:,2));'
      write(iw,*) ' c = transpose(x(:,3));'
      write (iw,*) ' z = zeros(size(xx));'
      write (iw,*) 'H=clinep(xx,yy,z,c,3)'

      write (iw,*) "axis ([-3 3, 0, 6.25]); "
      write (iw,*) "axis equal"
      write (iw,*) "colormap(new)"
      write (iw,*) "caxis([0 2])"
      write (iw,*) "colorbar"
      write (iw,*) "fname = sprintf('frame%.5d',", iframe,");"
      write (iw,*) "print('-dpng', '-r75', fname);"
ccc      write(iw, *) ' plot(xx,yy,\'o\')'
ccc      write(iw,*) 'axis equal'
ccc      write(iw, *) ' plot(x)'
      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
      end
C
C
        SUBROUTINE RSZ_MOVIE(Z,time,N,NBOD,k,iframe,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
        complex*16 Z(*)
        DIMENSION V(4)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        X,Y (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      WRITE(IW,*) 'x = ['
      WRITE(IW,1400) (Z(I), I=1,N)
      write(IW,1400) (Z(1))
ccc      write (iw,1400) x(1),y(1)
      WRITE(IW,*) '];'
      if (nbod.eq.k) then
         write (iw,1200) time
      write (iw,*) "title([' t= ',num2str(time)]);"
      end if
1200  FORMAT ('time=',F6.2)
1300  FORMAT(D15.6)
1400  FORMAT(2D15.6)
      if (NBOD.eq.1) write (iw,*) 'figure(1);clf;'
      write(iw,*) ' xx = transpose(x(:,1));'
      write(iw,*) ' yy = transpose(x(:,2));'
      write (iw,*) "H=plot(xx,yy,'LineWidth',2)"

      write (iw,*) "axis ([-3 3 0 6.25]); "
      write (iw,*) "axis equal"
      write (iw,*) "fname = sprintf('frame%.4d',", iframe,");"
      write (iw,*) "print('-dpng', '-r75', fname);"
ccc      write(iw, *) ' plot(xx,yy,\'o\')'
ccc      write(iw,*) 'axis equal'
ccc      write(iw, *) ' plot(x)'
      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
      end
CC
C
        SUBROUTINE RSCPLOT(Z,N,NBOD,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION V(4)
        complex*16 z(n)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        z = x + Iy (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      WRITE(IW,*) 'x = ['
      WRITE(IW,1400) (z(I),I=1,N)
ccc      write (iw,1400) x(1),y(1)
      WRITE(IW,*) '];'
1300  FORMAT(D15.6)
1400  FORMAT(2D15.6)
      write(iw,*) ' xx = x(:,1);'
      write(iw,*) ' yy = x(:,2);'
      write(iw, *) ' plot(xx,yy)'
ccc      write(iw, *) ' plot(xx,yy,\'o\')'
      write(iw,*) 'axis equal'
ccc      write(iw, *) ' plot(x)'
      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
      end
      