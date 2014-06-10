C     Reads grid file infile and writes it in gmsh format
C     infile is a single block grid file with three columns (x,y,z)
C     Author: Praveen C, http://math.tifrbng.res.in
C     Date:   16 April 2011
      program main
      parameter(imax=161, jmax=41)
      real*8 x(imax,jmax), y(imax,jmax)
      integer num(imax,jmax)
      integer count
      character*32 infile, outfile

      infile ='naca.struct'
      outfile='naca.msh'

      write(*,*)'Reading ',infile
      open(10, file=infile)
      read(10,*) ni, nj
      write(*,*) 'imax,jmax =', ni, nj
      read(10,*)((x(i,j),y(i,j),j=1,nj),i=1,ni)
      close(10)

c     Set node numbers
      count = 1
      do i=1,ni-1
         do j=1,nj
            num(i,j) = count
            count = count + 1
         enddo
      enddo

c     i index is cyclic
      do j=1,nj
         num(ni,j) = num(1,j)
      enddo

      write(*,*)'Writing ',outfile
      open(11, file=outfile)
      write(11,'("$MeshFormat")')
      write(11,'("2.1 0 8")')
      write(11,'("$EndMeshFormat")')
      write(11,'("$Nodes")')
      write(11,*) (ni-1)*nj
c     write coordinates
      do i=1,ni-1
         do j=1,nj
            write(11,*) num(i,j), x(i,j), y(i,j), 0.0
         enddo
      enddo
      write(11,'("$EndNodes")')
      write(11,'("$Elements")')
      write(11,*) (ni-1)*(nj-1) + 2*(ni-1)
c     write quads
      count = 1
      do i=1,ni-1
         do j=1,nj-1
            write(11,*) count, 3, 2, 2, 1, num(i,j), num(i+1,j),
     1                  num(i+1,j+1), num(i,j+1)
            count = count + 1
         enddo
      enddo
c     write boundary segments on airfoil
      do i=1,ni-1
         write(11,*) count, 1, 2, 0, 1, num(i,1), num(i+1,1)
         count = count + 1
      enddo
c     write boundary segments on outer boundary
      do i=1,ni-1
         write(11,*) count, 1, 2, 1, 1, num(i,nj), num(i+1,nj)
         count = count + 1
      enddo
      write(11,'("$EndElements")')
      close(11)

      stop
      end
