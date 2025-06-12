!!! Program for simulating three-state model for a single filament, with polymerization, hydrolysis & severing!!!
!!! This code has been developed by Soumyadipta Ray !!!

module params
    implicit none
    double precision, parameter:: pi = 4.0d0*datan(1.0d0)
    double precision, parameter:: Tfinal = 20000d0, Tint = 1.0d0  !! Tfinal is the time upto which, simulation would run. Tint is time interval of data print
    double precision,parameter::h0=0.15d0,h1=1.0d0,gamma_T=0.00d0,gamma_D=0.00d0,gamma_DC=0.0d0,s0_boundary=0.00001d0,&
                                &c=0.6d0,r0=11.6d0,r_slow=0.0019d0,r_fast=0.19d0,k_ne=0.0d0,s0_bulk=0.0d0,&
                                &c_actin=0.8d0,nh=3.5d0
    integer, parameter:: R=100 
    !! R is the spatial/geometric range, defining the range from +1 (ADP-Pi), within which if a -1 (ADP-cof) lies, &!
    !& the +1(ADP-Pi) would be 0 (ADP) !                        
    ! r_slow= standard rate of ADP-Pi to ADP conversion; r_fast= faster rate of ADP-Pi to ADP conversion &!
    !& if a ADP-Pi is within a range R from a cofilin decorated ADP !
    ! gamma_D, gamma_DC= depolymerization rate base on tip subunit; c= cofilin
    ! concn.; r0= polymerization rate; c_actin= actin concentarion, nh= hill
    ! coefficient; s0_boundary= severing rate when it happen at the junction decorated(ADP+cofilin) and&!
    !& undecorated(ADP-Pi) pair; !                                 
end module params

module variables
    implicit none
    double precision:: Tprint !! Data printing time
    double precision:: T !! Real time (stochastic variable)
    integer, allocatable:: l(:), l_dummy(:) !! l is our filament; l_dummy is a dummy filament. 
end module variables

!!!!!########################################################################################################################!!!!

!! Main program starts here !!

program single_filament_dyn
    use params
    use variables
    implicit none    
    integer:: i, j, reaction, l0, funit, ensmbl, idum, NT, ND, NDC, m1, p1, NG, o, s3, s4
    integer:: BndryT_DC_cnt, BndryD_DC_cnt, severing_sites_cnt !, LT_not_at_T_DC_bndry_count
    integer:: LT_outside_geometric_rnge_cnt, LT_within_geometric_rnge_cnt, min_pos, max_pos
    double precision:: A(7), Atot, sum_A, tao, r1, r2, t1, t2, ran2, k1
    integer, allocatable:: LT(:), LD(:), LDC(:), LBndryT_DC(:), LBndryD_DC(:)
    !integer, allocatable:: LT_at_T_DC_bndry(:), LT_not_at_T_DC_bndry(:), 
    integer, allocatable:: LT_within_geometric_rnge(:), LT_outside_geometric_rnge(:), severing_sites_pos(:)
    character(len=100):: filename, filepath

    call cpu_time(t1)

    ensemble_loop: do ensmbl=2,4 
        filepath = 'avgd_data/with_geometric_factor/ho0.15_rslow0.0019_rfast0.19_cActin0.8_cof0.6_varyingR/'
write (filename,"('length_vs_time_actin0.8_h00.15_rslow0.0019_rfast0.19_R100_c0.6_gt0_gd0_2e4ite_cnf_',I0,'.txt')") ensmbl
        open(funit, file=trim(adjustl(filepath))//filename, status='unknown', position='append', action='write')

        NG = 0     !Number of ADP/ADP-cof subunits in th pool tracker
        T = 0.0d0  ! Initial real time
        Tprint = T ! Printing time is from the starting

        CALL SYSTEM_CLOCK(COUNT=idum)
        !r1 = ran2(idum) !! uniform random number between 0 & 1
        !r1 = ran2(idum)
        !l0 = nint(10 + (30-10)*r1) !! initial filament length
        !write(*,*) 10 + (30-10)*r1, r1, idum
        l0=40 ! Initial filament size
        if (allocated(l)) deallocate(l)
        if (allocated(l_dummy)) deallocate(l_dummy)
        allocate(l(l0))  !! Allocation of the filament array of size l0
        allocate(l_dummy(l0))
        l = 1 !! Assigning all array elements as +1 (ADP-Pi) of the initial filament
              !! Here, we denote ADP-Pi as +1, ADP as 0, and ADP-cof as -1 !!
        
        time_loop: do 

            if (allocated(LT)) deallocate(LT)
            allocate(LT(size(l)))
            if (allocated(LDC)) deallocate(LDC)
            allocate(LDC(size(l)))
            if (allocated(LD)) deallocate(LD)
            allocate(LD(size(l)))
            if (allocated(LBndryT_DC)) deallocate(LBndryT_DC)
            allocate(LBndryT_DC(size(l)))            
            if (allocated(LBndryD_DC)) deallocate(LBndryD_DC)
            allocate(LBndryD_DC(size(l)))            
            if (allocated(LT_within_geometric_rnge)) deallocate(LT_within_geometric_rnge)
            allocate(LT_within_geometric_rnge(size(l)))    
            if (allocated(LT_outside_geometric_rnge)) deallocate(LT_outside_geometric_rnge)
            allocate(LT_outside_geometric_rnge(size(l)))                        


            LT = 0                                ! ATP array: collects the position of ADP-Pi(+1) in the filament
            LDC = 0                               ! ADP-cof array: collects the position of ADP-cof(-1) in the filament 
            LD = 0                                ! ADP array: collects the position of ADP(0) in the filament 
            LBndryT_DC = 0                        ! Array of the ADP-Pi & ADP-cof boundary positions
            LBndryD_DC = 0                        ! Array of the ADP-ADPcof boundary positions 
            LT_within_geometric_rnge = 0
            LT_outside_geometric_rnge = 0
            NT = 0                                ! Tracker for counting the number of ADP-Pi monomers in the lattice (array) of the filamnet 
            NDC = 0                               ! Tracker for counting the number of ADP-cof monomers in the lattice (array) of the filamnet
            ND = 0                                ! Tracker for counting the number of ADP monomers in the lattice (array) of the filamnet
            BndryT_DC_cnt = 0                     ! count boundary nos of +1 and -1 and vice versa
            BndryD_DC_cnt = 0                     ! count boundary nos of 0 and -1 and vice versa
            LT_within_geometric_rnge_cnt = 0
            LT_outside_geometric_rnge_cnt = 0
            severing_sites_cnt = 0

            
            !!!#### total number & positional information of +1(D-Pi), 0(D), and -1(DC) #####!!!
            if( size(l) > 0 ) then                        ! Tracking the Number of ADP-Pi, only when filament length is >0
                do i = 1,size(l)
                    if (l(i)==1)then
                        NT = NT + 1
                        LT(NT)=i
                    else if (l(i) == -1) then             !Tracking the Number of ADP-cof, only when filament length is >0
                        NDC = NDC + 1
                        LDC(NDC)=i
                    else if (l(i) == 0) then             !Tracking the Number of ADP, only when filament length is >0
                        ND = ND + 1
                        LD(ND)=i    
                    end if
                end do
            end if
            !!!#### total number & positional information of +1(D-Pi), 0(D), and -1(DC) #####!!!

            
            !!!#### total number & positional info of +1/-1, and 0/-1 boundaries, and pos info of those +1 that are at +1/-1 boundary #####!!!
            if( size(l) > 1 ) then   
                do i=1,size(l)-1
                    if (abs(l(i)-l(i+1)).eq.2) then
                         BndryT_DC_cnt = BndryT_DC_cnt + 1
                         LBndryT_DC(BndryT_DC_cnt) = i               ! Keep track of  boundary of +1 and -1 and vice versa
                    else if ((l(i)+l(i+1)).eq.-1) then
                        BndryD_DC_cnt = BndryD_DC_cnt + 1
                        LBndryD_DC(BndryD_DC_cnt) = i               ! Keep track of  boundary of 0 and -1 and vice versa                            
                    end if
                   end do
            end if
            !!!#### total number & positional info of +1/-1, and 0/-1 boundaries, and pos info of those +1 that are at +1/-1 boundary #####!!!

            

            !!!### Total number & positional info of those +1(D-Pi) which lie within the geometric range(R) of a -1(DC) ###!!!
            if (NT.gt.0) then
                loop_over_LT: do i=1,NT
                                if (LT(i).lt.size(l)) then         !! strat to search in the LT(i)+R domain
                                    min_pos = LT(i) + 1
                                    max_pos = LT(i) + R
                                    if (max_pos.gt.size(l)) max_pos = size(l)

                                    if (any((LDC .ge. min_pos).and.(LDC .le. max_pos))) then
                                        LT_within_geometric_rnge_cnt = LT_within_geometric_rnge_cnt + 1
                                        LT_within_geometric_rnge(LT_within_geometric_rnge_cnt) = LT(i)  ! pos. inf. of the +1 those are within the geometric rnge
                                        cycle loop_over_LT
                                    end if    
                                end if

                                if (LT(i).gt.1) then         !! start to search in the LT(i)-R domain
                                    min_pos = LT(i) - R
                                    max_pos = LT(i) - 1
                                    if (min_pos.lt.1) min_pos = 1

                                    if (any((LDC .ge. min_pos).and.(LDC .le. max_pos))) then
                                        LT_within_geometric_rnge_cnt = LT_within_geometric_rnge_cnt + 1
                                        LT_within_geometric_rnge(LT_within_geometric_rnge_cnt) = LT(i)  !! pos. inf. of the +1 those are within the geometric rnge
                                        cycle loop_over_LT
                                    end if    
                                end if

                                LT_outside_geometric_rnge_cnt = LT_outside_geometric_rnge_cnt + 1
                                LT_outside_geometric_rnge(LT_outside_geometric_rnge_cnt) = LT(i)   !! pos. inf. of the +1 those are outside the geometric rnge
                            end do loop_over_LT                
            end if        
            !!!### Total number & positional info of those +1(D-Pi) which lie within the geometric range(R) of a -1(DC) ###!!!


            
            !!!#### total number & positional info of severing sites, i.e. +1/-1 and 0/-1 boundaries #####!!!
            if(size(l) > 1) then
                severing_sites_cnt = BndryD_DC_cnt + BndryT_DC_cnt      ! total number of severing sites
                if (allocated(severing_sites_pos)) deallocate(severing_sites_pos)
                allocate(severing_sites_pos(severing_sites_cnt))        ! array that would contain the positions of severing-sites        
                severing_sites_pos(1:BndryT_DC_cnt) = LBndryT_DC(1:BndryT_DC_cnt)
                severing_sites_pos(BndryT_DC_cnt+1:severing_sites_cnt) = LBndryD_DC(1:BndryD_DC_cnt)
            end if
            !!!#### total number & positional info of severing sites, i.e. +1/-1 and 0/-1 boundaries #####!!!

            
            !!!#### Calculation of the propensities ####!!!
            A(1) = r0*c_actin                            ! propensity for assembly; !! for finite pool, it would be proportional to (N-size(l)-NG).
            !A(3) = r_slow*LT_not_at_T_DC_bndry_count     ! hydrolysis at a slower rate
            A(3) = r_slow*LT_outside_geometric_rnge_cnt     ! hydrolysis at a slower rate
            !A(4) = r_fast*BndryT_DC_cnt                  ! hydrolysis at a faster rate
            A(4) = r_fast*LT_within_geometric_rnge_cnt    ! hydrolysis at a faster rate
            A(5) = ((h0*c**nh)/(h1 + c**nh))*(ND)        ! cofilin addition with ADP, the rate being a hill function
            A(6) = s0_boundary*severing_sites_cnt        ! propensity for severing at D-Pi/DC or D/DC boundaries             
            A(7) = k_ne*NG                               ! propensity for reverse hydrolysis in pool

            if (size(l) > 0 .and. l(size(l)) == 1 ) then          !Conditioning the rates with the length values
                A(2) = gamma_T    
            else if ( size(l) > 0 .and. l(size(l)) == 0) then
                A(2) = gamma_D
            else if ( size(l) > 0 .and. l(size(l)) == -1) then
                A(2) = gamma_DC    
            end if

            Atot = sum(A) !! Total propensity
            !!!#### Calculation of the propensities ####!!!


            CALL SYSTEM_CLOCK(COUNT=idum)
            r1 = ran2(idum) 
            r2 = ran2(idum) 
            tao = (1.0d0/Atot)*dlog(1/r1)
            T = T + tao
            if(T.gt.Tfinal) then
                write(*,*) T, tao, Tprint
                exit time_loop
            end if 
            
            r2 = r2*Atot
            sum_A = 0
            reaction_choice: do reaction=1, size(A)
                sum_A = sum_A + A(reaction)
                if(sum_A.ge.r2) exit reaction_choice
            end do reaction_choice

            
            select case (reaction)

                case(1)                                           ! Assembly at the left end
                    if (allocated(l_dummy)) deallocate(l_dummy)
                    allocate(l_dummy(size(l)))
                    l_dummy = l
                    deallocate(l)
                    allocate(l(size(l_dummy)+1))
                    l(2:size(l_dummy)) = l_dummy
                    l(1) = 1                      
                    !deallocate(l_dummy)
                
                case(2)                                             ! Decay at right end
                    if (allocated(l_dummy)) deallocate(l_dummy)  
                    allocate(l_dummy(size(l)))
                    l_dummy = l
                    deallocate(l)
                    allocate(l(size(l_dummy)-1))
                    l(:) = l_dummy(1:size(l_dummy)-1)
                    if ((l_dummy(size(l_dummy)) == -1).or.(l_dummy(size(l_dummy)) == 0)) NG = NG + 1  ! counting number of D or DC monomer in the pull.
            
                case(3)
                    k1 = ((LT_outside_geometric_rnge_cnt - 1)*ran2(idum) + 1) !Calling a random number [1,LT_not_at_T_DC_bndry_count], based upon the total number of D-Pi that are not at D-Pi/DC boundary
                    m1 = nint(k1)
                    p1 = LT_outside_geometric_rnge(m1)                     !Hydrolysis (at r_slow rate)
                    l(p1) = 0                                              !Inserting a 0(D monomer) in place of that +1(D-Pi monomer) 

                case(4)
                    k1 = ((LT_within_geometric_rnge_cnt - 1)*ran2(idum) + 1) !Calling a random number [1,LT_not_at_T_DC_bndry_count], based upon the total number of D-Pi that are not at D-Pi/DC boundary
                    m1 = nint(k1)
                    p1 = LT_within_geometric_rnge(m1)          !Hydrolysis (at r_fast rate)
                    l(p1) = 0                                  !Inserting a 0(D monomer) in place of that +1(D-Pi monomer) 

                case(5)                                        ! Cofilin addition: (0)D to (-1)DC conversion
                    k1 = ((ND - 1)*ran2(idum) + 1)
                    m1 = nint(k1)
                    p1 = LD(m1) 
                    l(p1) = -1                                 !Inserting a -1(DC monomer) in place of 0(D monomer) 

                case(6)                                                !severing at +1/-1 or 0/-1 boundary
                    o = nint(ran2(idum)*(severing_sites_cnt - 1) + 1)        ! Choose random no. betn [1,severing_sites_cnt]
                    s3 = severing_sites_pos(o)                                  ! s3 is the site where filament breaks
                    s4 = 0                                             ! initialisation; counts -1 right to the severing site
                                 
                    if (allocated(l_dummy)) deallocate(l_dummy)
                    allocate(l_dummy(size(l)))
                    l_dummy = l
                    deallocate(l)
                    allocate(l(s3))
                    l(:) = l_dummy(1:s3)

                    do i=s3+1,size(l_dummy)    
                    if ((l_dummy(i) == -1).or.(l_dummy(i) == 0)) s4 = s4 +1  ! counting number of D or DC monomers, that fall into the pull due to severing.
                    end do
                   
                    NG =  NG +  s4                                    ! counting number of D or DC monomer in the pull.  

                case(7)                                                !Nucleotide exchange in pool 
                    NG = NG - 1                                        ! reverse hydrolysis
                    
                case default
                    write(*,*) 'error in choosing reaction; reaction=', reaction
                    error stop 1
            end select        

           
            if(T.ge.Tprint) then
                write(funit,*) Tprint, size(l), NT, ND, NDC, severing_sites_cnt
                Tprint = Tprint + Tint
            end if    

        end do time_loop

    end do ensemble_loop    

    close(funit)
    call cpu_time(t2)
    write(*,*) t2-t1

end program single_filament_dyn

!!!!!########################################## FUNCTION AND SUBROUTINES USED #############################################!!!!!

!!!!/////// Uniform Random number generators////////////////////////////////////

FUNCTION ran2(idum)
    ! USE numz
     IMPLICIT NONE
     DOUBLE PRECISION:: ran2
     !INTEGER,INTENT(inout),OPTIONAL::idum
     INTEGER,INTENT(inout)::idum
     !INTEGER :: idum
     INTEGER,PARAMETER::IM1=2147483563,IM2=2147483399,IMM1=IM1-1
     INTEGER,PARAMETER::IA1=40014,IA2=40692,IQ1=53668
     INTEGER,PARAMETER::IQ2=52774,IR1=12211,IR2=3791   
     INTEGER,PARAMETER::NTAB=32,NDIV=1+IMM1/NTAB
     DOUBLE PRECISION,PARAMETER::AM=1.0d0/IM1,EPS=1.2e-7,RNMX=1.0d0-EPS
     INTEGER::idum2,j,k,iv(NTAB),iy
     SAVE iv,iy,idum2
     DATA idum2/123456789/, iv/NTAB*0/, iy/0/
     IF (idum<0) THEN
        idum=MAX(-idum,1)
        idum2=idum
         DO j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            IF (idum<0) idum=idum+IM1
            IF (j.LE.NTAB) iv(j)=idum
         ENDDO
         iy=iv(1)
      ENDIF
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      IF (idum<0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      IF (idum2<0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      IF(iy.LT.1)iy=iy+IMM1
      ran2=MIN(AM*iy,RNMX)
      RETURN
    END FUNCTION ran2
   

