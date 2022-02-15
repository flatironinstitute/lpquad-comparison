      implicit real *8 (a-h,o-z)
      real *8 a(3,4)
      integer i1,i2,i3
      real *8 b(5,2)

      i1 = 7
      i2 = 3
      i3 = 1
      
      do i=1,3
        do j=1,4
          a(i,j) = 1
        enddo
      enddo

      do i=1,5
        do j=1,2
          b(i,j) = -1.1d0
        enddo
      enddo

      a(2,3) = 3
      b(3,1) = -5

      irl = 4+4+4+(12+10)*8
      

      open(unit=33,file='tunfw.dat',form='unformatted')
      write(33) i1
      write(33) i2
      write(33) i3
      write(33) a
      write(33) b
      close(33)



      stop
      end
