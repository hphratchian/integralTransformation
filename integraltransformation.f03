INCLUDE 'integraltransformation_mod.f03'
      Program integraltransformation
!
!     This program reads AO integrals from a Gaussian matrix file and times
!     AO-to-MO integral transformations.
!
!     -H. P. Hratchian, 2022.
!
!
!     USE Connections
!
      use integraltransformation_mod
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64)::nCommands,iPrint=0,nAtoms,nBasis,nBasisUse,  &
        nElectrons,nElectronsAlpha,nElectronsBeta
      integer(kind=int64)::mu,nu,lambda,sigma,p,q,r,s,pq,rs,pqrs,iCount,i,j,a,b
      real(kind=real64)::time0,time1,deltaIJAB,numerator,E2AA,E2BB,E2AB,E2BA
      real(kind=real64),dimension(:),allocatable::moEnergiesAlpha,moEnergiesBeta
      real(kind=real64),dimension(:,:),allocatable::CAlpha,CBeta
      real(kind=real64),dimension(:,:,:,:),allocatable::aoInts,moInts,  &
        partialInts1,partialInts2
      type(MQC_Variable)::ERIs,mqcTmpArray
      character(len=512)::tmpString,matrixFilename
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      logical::fail=.false.,doN8=.false.,doSlowN5=.true.
!
!     Format Statements
!
 1000 Format(1x,'Enter Test Program integralTransformation.')
 1010 Format(3x,'Matrix File: ',A,/)
 1100 Format(1x,'nAtoms    =',I4,6x,'nBasis  =',I4,6x,'nBasisUse=',I4,/,  &
             1x,'nElectrons=',I4,6x,'nElAlpha=',I4,6x,'nElBeta  =',I4)
 1500 Format(/,1x,'Carrying out O(N^8) transformation.')
 1510 Format(/,1x,'Skipping O(N^8) transformation.')
 2000 Format(1x,'<',I3,',',I3,' || ',I3,',',I3,' > ... pq=',I3,'  rs=',I3,'  pqrs=',I3)
 3000 Format(/,1x,'E(2)-SS = ',f15.10,' a.u.',4x,'E(2)-OS = ',f15.10,' a.u.')
 5000 Format(1x,'Time (',A,'): ',f8.1,' s.')
 8999 Format(/,1x,'END OF PROGRAM integralTransformation.')
!
!
      write(IOut,1000)
!
!     Open the Gaussian matrix file and load the number of atomic centers.

      nCommands = command_argument_count()
      if(nCommands.eq.0)  &
        call mqc_error('No command line arguments provided. The input Gaussian matrix file name is required.')
      call get_command_argument(1,tmpString)
      call commandLineArgs(iPrint,matrixFilename,doN8,doSlowN5,fail)

      write(iOut,*)
      write(iOut,*)

      write(iOut,*)' iPrint = ',iPrint
      write(iOut,*)' matrixFilename = ',TRIM(matrixFilename)
      write(iOut,*)' doN8 = ',doN8
      write(iOut,*)' doSlowN5 = ',doSlowN5
      write(iOut,*)' fail = ',fail

      write(iOut,*)
      write(iOut,*)


      call GMatrixFile%load(matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      nAtoms = GMatrixFile%getVal('nAtoms')
      nBasis = Int(GMatrixFile%getVal('nbasis'))
      nBasisUse = Int(GMatrixFile%getVal('nbasisuse'))
      nElectrons = Int(GMatrixFile%getVal('nelectrons'))
      nElectronsAlpha = Int(GMatrixFile%getVal('nAlpha'))
      nElectronsBeta = Int(GMatrixFile%getVal('nBeta'))
      write(IOut,1100) nAtoms,nBasis,nBasisUse,nElectrons,  &
        nElectronsAlpha,nElectronsBeta
!
!     Load the orbital eigenvalues.
!
      call GMatrixFile%getArray('ALPHA ORBITAL ENERGIES',mqcVarOut=mqcTmpArray)
      call mqcTmpArray%print(header='Alpha MO Energies')
      moEnergiesAlpha = mqcTmpArray
      if(GMatrixFile%isUnrestricted()) then
        call mqc_error('UHF/UKS NYI.')
      else
        moEnergiesBeta = moEnergiesAlpha
      endIf
!
!     Load the MO coefficients.
!
      write(*,*)' Hrant -  isUnrestricted: ',GMatrixFile%isUnrestricted()
      call GMatrixFile%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=mqcTmpArray)
      CAlpha = mqcTmpArray
      if(GMatrixFile%isUnrestricted()) then
        call GMatrixFile%getArray('BETA  MO COEFFICIENTS',mqcVarOut=mqcTmpArray)
        CBeta = mqcTmpArray
      else
        CBeta = CAlpha
      endIf
!
!     Read in and report out the (AO) ERIs.
!
      call GMatrixFile%getArray('REGULAR 2E INTEGRALS',mqcVarOut=ERIs)
      if(iPrint.ge.2) call ERIs%print(IOut,' ERIs=')
!
!     Allocate space for AO and MO integrals. Then fill the intrinsic array of
!     AO integrals.
!
      allocate(aoInts(nBasis,nBasis,nBasis,nBasis),  &
        moInts(nBasisUse,nBasisUse,nBasisUse,nBasisUse))
      aoInts = ERIs
      if(iPrint.ge.2) call mqc_print_rank4Tensor_array_real(iOut,  &
        aoInts,header='Intrinsic AO Integrals')
!
!     Do N^8 AO --> MO transformation.
!
      if(doN8) then
        write(iOut,1500)
        call cpu_time(time0)
        moInts = float(0)
        call cpu_time(time1)
        write(iOut,5000) 'Init moInts',time1-time0
        flush(iOut)
!
        call cpu_time(time0)


        do p = 1,nBasisUse
          do q = 1,nBasisUse
            do r = 1,nBasisUse
              do s = 1,nBasisUse
                do mu = 1,nBasis
                  do nu = 1,nBasis
                    do lambda = 1,nBasis
                      do sigma = 1,nBasis
                        moInts(p,q,r,s) = moInts(p,q,r,s)  &
                          + CAlpha(mu,p)*CAlpha(nu,q)*CAlpha(lambda,r)*CAlpha(sigma,s)*aoInts(mu,nu,lambda,sigma)
                      endDo
                    endDo
                  endDo
                endDo
              endDo
            endDo
          endDo
        endDo



        call cpu_time(time1)
        write(iOut,5000) 'O(N^8) Transformation',time1-time0
      else
        write(iOut,1510)
      endIf
      flush(iOut)
!
!     Do quarter transformations now.
!
      moInts = float(0)
      Allocate(partialInts1(nBasisUse,nBasis,nBasis,nBasis))
      partialInts1 = float(0)
!
      if(doSlowN5) then
        call cpu_time(time0)
        do p = 1,nBasisUse
          do mu = 1,nBasis
            do nu = 1,nBasis
              do lambda = 1,nBasis
                do sigma = 1,nBasis
                  partialInts1(p,nu,lambda,sigma) = partialInts1(p,nu,lambda,sigma)  &
                    + CAlpha(mu,p)*aoInts(mu,nu,lambda,sigma)
                endDo
              endDo
            endDo
          endDo
        endDo
        call cpu_time(time1)
        write(iOut,5000) 'Quarter Transformation 1a',time1-time0
        flush(iOut)
        partialInts1 = float(0)
      endIf
!
      call cpu_time(time0)
      do nu = 1,nBasis
        do lambda = 1,nBasis
          do sigma = 1,nBasis
            do p = 1,nBasisUse
              do mu = 1,nBasis
                partialInts1(p,nu,lambda,sigma) = partialInts1(p,nu,lambda,sigma)  &
                  + CAlpha(mu,p)*aoInts(mu,nu,lambda,sigma)
              endDo
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)
      write(iOut,5000) 'Quarter Transformation 1b',time1-time0
      flush(iOut)
!
      Allocate(partialInts2(nBasisUse,nBasisUse,nBasis,nBasis))
      partialInts2 = float(0)
!
      if(doSlowN5) then
        call cpu_time(time0)
        do p = 1,nBasisUse
          do q = 1,nBasisUse
            do nu = 1,nBasis
              do lambda = 1,nBasis
                do sigma = 1,nBasis
                  partialInts2(p,q,lambda,sigma) = partialInts2(p,q,lambda,sigma)  &
                    + CAlpha(nu,q)*partialInts1(p,nu,lambda,sigma)
                endDo
              endDo
            endDo
          endDo
        endDo
        call cpu_time(time1)
        write(iOut,5000) 'Quarter Transformation 2a',time1-time0
        flush(iOut)
        partialInts2 = float(0)
      endIf
!
      call cpu_time(time0)
      do lambda = 1,nBasis
        do sigma = 1,nBasis
          do q = 1,nBasisUse
            do p = 1,nBasisUse
              do nu = 1,nBasis
                partialInts2(p,q,lambda,sigma) = partialInts2(p,q,lambda,sigma)  &
                  + CAlpha(nu,q)*partialInts1(p,nu,lambda,sigma)
              endDo
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)
      write(iOut,5000) 'Quarter Transformation 2b',time1-time0
      flush(iOut)
!
      DeAllocate(partialInts1)
      Allocate(partialInts1(nBasisUse,nBasisUse,nBasisUse,nBasis))
      partialInts1 = float(0)
!
      if(doSlowN5) then
        call cpu_time(time0)
        do p = 1,nBasisUse
          do q = 1,nBasisUse
            do r = 1,nBasisUse
              do lambda = 1,nBasis
                do sigma = 1,nBasis
                  partialInts1(p,q,r,sigma) = partialInts1(p,q,r,sigma)  &
                    + CAlpha(lambda,r)*partialInts2(p,q,lambda,sigma)
                endDo
              endDo
            endDo
          endDo
        endDo
        call cpu_time(time1)
        write(iOut,5000) 'Quarter Transformation 3a',time1-time0
        flush(iOut)
        partialInts1 = float(0)
      endIf
!
      call cpu_time(time0)
      do p = 1,nBasisUse
        do q = 1,nBasisUse
              do sigma = 1,nBasis
          do r = 1,nBasisUse
            do lambda = 1,nBasis
                partialInts1(p,q,r,sigma) = partialInts1(p,q,r,sigma)  &
                  + CAlpha(lambda,r)*partialInts2(p,q,lambda,sigma)
              endDo
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)
      write(iOut,5000) 'Quarter Transformation 3b',time1-time0
      flush(iOut)
!
      DeAllocate(partialInts2)
      moInts = float(0)
!
      call cpu_time(time0)
      do p = 1,nBasisUse
        do q = 1,nBasisUse
          do r = 1,nBasisUse
            do s = 1,nBasisUse
              do sigma = 1,nBasis
                moInts(p,q,r,s) = moInts(p,q,r,s)  &
                  + CAlpha(sigma,s)*partialInts1(p,q,r,sigma)
              endDo
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)
      write(iOut,5000) 'Quarter Transformation 4',time1-time0
      flush(iOut)
      DeAllocate(partialInts1)
      if(iPrint.ge.1) call mqc_print_rank4Tensor_array_real(iOut,  &
        moInts,header='Transformed MO Integrals 2')
      iCount = 0
      do p = 1,nBasisUse
        do q = 1,p
          pq = (p*(p-1))/2 + q
          do r = 1,nBasisUse
            do s = 1,r
              rs = (r*(r-1))/2 + s
              if(pq.ge.rs) then
                pqrs = (pq*(pq-1))/2 + rs
                iCount = iCount + 1
              endIf
            endDo
          endDo
        endDo
      endDo
      call GMatrixFile%getArray('AA MO 2E INTEGRALS',mqcVarOut=mqcTmpArray)
      if(iPrint.ge.1) call mqcTmpArray%print(header='MO Ints from Matrix File')
!
!     Evaluate the E(2) AA, BB, AB, and BA contribution.
!
      write(*,*)
      write(*,*)' Same Spin E2...'
      E2AA = float(0)
      call cpu_time(time0)
      do i = 1,nElectronsAlpha
        do j = 1,nElectronsAlpha
          do a = nElectronsAlpha+1,nBasisUse
            do b = nElectronsAlpha+1,nBasisUse
              deltaIJAB = moEnergiesAlpha(i) + moEnergiesAlpha(j)  &
                - moEnergiesAlpha(a) - moEnergiesAlpha(b)
              numerator = moInts(i,a,j,b) - moInts(i,b,j,a)
              numerator = numerator*numerator
              if(iPrint.ge.1) write(*,*)' num, denom = ',numerator,deltaIJAB
              E2AA = E2AA + numerator/(float(4)*deltaIJAB)
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)
      write(iOut,5000) 'Same Spin E2',time1-time0
      E2BB = E2AA
      write(*,*)
      write(*,*)' Opposite Spin E2...'
      E2AB = float(0)
      call cpu_time(time0)
      do i = 1,nElectronsAlpha
        do j = 1,nElectronsBeta
          do a = nElectronsAlpha+1,nBasisUse
            do b = nElectronsBeta+1,nBasisUse
              deltaIJAB = moEnergiesAlpha(i) + moEnergiesBeta(j)  &
                - moEnergiesAlpha(a) - moEnergiesBeta(b)
              numerator = moInts(i,a,j,b)*moInts(i,a,j,b)
              if(iPrint.ge.1) write(*,*)' num, denom = ',numerator,deltaIJAB
              E2AB = E2AB + numerator/deltaIJAB
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)
      write(iOut,5000) 'Opposite Spin E2',time1-time0
      E2BA = E2AB




      write(iOut,3000) E2AA,E2AB
      write(*,*)' Total E2 = ',E2AA+E2BB+E2AB






!
  999 Continue
      write(iOut,8999)
      end program integraltransformation
