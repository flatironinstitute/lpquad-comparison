      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targo(:,:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      real *8 errs(6),ts(2)
      real *8, allocatable :: rfacs(:,:)
      complex *16 zk
      character *100 fname
      integer ipars(2)
      integer, allocatable :: row_ptr(:),col_ind(:)
      integer, allocatable :: iquad(:)
      integer, allocatable :: col_ptr(:),row_ind(:),iper(:)
      integer, allocatable :: irowind(:),icolind(:)


      real *8, allocatable :: slp_near(:),dlp_near(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      integer, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:),targs(:,:)
      real *8 xyz_out(3),xyz_in(3)
      real *8, allocatable :: sigma(:)
      real *8 dpars(2)
      character *200 fname1,fname2


      call prini(6,13)

      done = 1
      pi = atan(done)*4

c
c  levels of refinement
c
      iref = 2

c
c
c   average aspect ratio
c     iasp = 1, aavg = 1.6
c     iasp = 2, aavg = 2.3
c     iasp = 3, aavg = 4.6
c     iasp = 4, aavg = 9.3
c     iasp = 5, aavg = 14.0 

      iasp = 1


c       select geometry type
c       igeomtype = 1 => sphere
c       igeomtype = 2 => stellarator
c 

      igeomtype = 2

      if(iasp.eq.1) then
        ipars(1) = 5*(2**(iref))
        ipars(2) = 15*(2**(iref))

      else if(iasp.eq.2) then
        ipars(1) = 6*(2**(iref))
        ipars(2) = 12*(2**(iref))

      else if(iasp.eq.3) then
        ipars(1) = 9*(2**(iref))
        ipars(2) = 9*(2**(iref))

      else if(iasp.eq.4) then
        ipars(1) = 12*(2**(iref))
        ipars(2) = 6*(2**(iref))
      else
        ipars(1) = 15*(2**(iref))
        ipars(2) = 5*(2**(iref))
      endif

      

      igeomtype = 2

      npatches = 2*ipars(1)*ipars(2)


      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

      if(igeomtype.eq.1) then
        xyz_out(1) = 3.17d0
        xyz_out(2) = -0.03d0
        xyz_out(3) = 3.15d0

        xyz_in(1) = 0.17d0
        xyz_in(2) = 0.23d0
        xyz_in(3) = -0.11d0
      endif

      if(igeomtype.eq.2) then
        xyz_in(1) = -4.501d0
        xyz_in(2) = 1.7d-3
        xyz_in(3) = 0.00001d0

        xyz_out(1) = -3.5d0
        xyz_out(2) = 3.1d0
        xyz_out(3) = 20.1d0
      endif

      norder = 3 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(targs(12,npts))
      ifplot = 0



      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)



      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)

      ndtarg = 12
     
      do i=1,npts
        do j=1,12
          targs(j,i) = srcvals(j,i)
        enddo
      enddo

      allocate(ipatch_id(npts),uvs_targ(2,npts))
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, 
     1         ipatch_id,uvs_targ)

 
c
c    find near field
c
      iptype = 1
      call get_rfacs(norder,iptype,rfac,rfac0)
      do i=1,npatches 
        rad_near(i) = rads(i)*rfac
      enddo
      

      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)
      ntarg = npts
      allocate(row_ind(nnz),col_ptr(npatches+1),iper(nnz))
      call rsc_to_csc(npatches,ntarg,nnz,row_ptr,col_ind,col_ptr,
     1  row_ind,iper)

      nquad = iquad(nnz+1)-1
      allocate(slp_near(nquad),dlp_near(nquad))


      ndtarg = 12

      eps = 0.50001d-3

      ikerorder = -1
      zk = 0


      write(fname1,'(a,i1.1,a,i1.1,a,i1.1,a)')
     1   './data/stell_norder',norder,'_iref',iref,'_iasp',iasp,'.bin'
      write(fname2,'(a,i1.1,a,i1.1,a,i1.1,a)')
     1   './data/stell_quad_norder',norder,'_iref',iref,'_iasp',
     2   iasp,'.bin'

      open(unit=33,file=trim(fname1),form='unformatted')
      write(33) npatches
      write(33) npts
      write(33) nnz
      write(33) ndtarg
      write(33) npts
      write(33) nquad
      write(33) norders
      write(33) ixyzs
      write(33) iptype
      write(33) srcvals
      write(33) srccoefs
      write(33) wts
      write(33) targs
      write(33) ipatch_id
      write(33) uvs_targ
      write(33) row_ptr
      write(33) col_ind
      write(33) iquad
      write(33) col_ptr
      write(33) row_ind
      write(33) iper
      close(33)
      

      do i=1,nquad
        slp_near(i) = 0
        dlp_near(i) = 0
      enddo



      call cpu_time(t1)

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

      iquadtype = 1

      call getnearquad_lap_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,dpars,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,slp_near)

      
      dpars(1) = 0.0d0
      dpars(2) = 1.0d0
      call getnearquad_lap_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,dpars,iquadtype,
     1      nnz,row_ptr,col_ind,iquad,
     1      rfac0,nquad,dlp_near)
      
      call cpu_time(t2)
      tquadgen = t2-t1

      allocate(irowind(nquad),icolind(nquad))
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            irowind(jquadstart+l-1) = i
            icolind(jquadstart+l-1) = jstart+l-1
          enddo
        enddo
      enddo
      open(unit=34,file=trim(fname2),form='unformatted')
      write(34) nquad
      write(34) irowind
      write(34) icolind
      write(34) slp_near
      write(34) dlp_near
      close(34)
      
      stop
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri, 
     1      triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         iptr3,iptr4, norder,
     2         'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end

