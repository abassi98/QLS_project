   program main
     use debugging
     use bacteria
      implicit none
      integer*16 :: N_S, N_R, N_T, count,status, ii,jj
      double precision :: dt, threshold, theta, T, square_diff,tt
      double precision, dimension(:), allocatable ::m,c,k_n,k_t,rho,s,xi,K,tau,PHI_big,q
      double precision, dimension(:), allocatable ::r, temp_m,temp_c
      double precision, dimension(:,:),allocatable :: temp_phi
      double precision, dimension(:,:), allocatable::gamma_T,gamma,eta_T,eta,phi_T,phi
      logical :: debug, printer

      count=0
!!!!!!!! Parameter reading !!!!!!!!
      !Read debugging variables
      open(10,file="temp/debugging.txt",status='unknown')
      read(10,*) debug, printer
      close(10)

      !Read number of hours of simulation
      open(10,file="temp/parameters.txt",status='unknown')
      read(10,*) N_S,N_R,N_T,dt,threshold,theta
      close(10)


      !Allocate statical variables
      allocate(k_n(N_R),stat=status) !nutritional capacity
      call check_allocation(status,debug,printer)
      allocate(k_t(N_S),stat=status) !traslational capacity
      call check_allocation(status,debug,printer)
      allocate(gamma(N_S,N_R),stat=status) !ratio between them
      call check_allocation(status,debug,printer)
      allocate(gamma_T(N_R,N_S),stat=status) !ratio between them
      call check_allocation(status,debug,printer)
      allocate(rho(N_S),stat=status) !conversion factor
      call check_allocation(status,debug,printer)
      allocate(eta(N_S,N_R),stat=status) !eta parameter
      call check_allocation(status,debug,printer)
      allocate(eta_T(N_R,N_S),stat=status) !eta parameter
      call check_allocation(status,debug,printer)
      allocate(s(N_R), stat=status) !constant resource injection rate
      call check_allocation(status,debug,printer)
      allocate(xi(N_R), stat=status) !maximum catalytic rate
      call check_allocation(status,debug,printer)
      allocate(K(N_R),stat=status) !half saturation constant
      call check_allocation(status,debug,printer)
      allocate(tau(N_S),stat=status) !characteristic timescale of adaptive process
      call check_allocation(status,debug,printer)
      allocate(PHI_big(N_S),stat=status) !total  proteome fraction for groeth and met
      call check_allocation(status,debug,printer)
      allocate(q(N_S),stat=status) !death rate
      call check_allocation(status,debug,printer)

      !Allocate dynamical  variables
      allocate(m(N_S),stat=status) !species
      call check_allocation(status,debug,printer)
      allocate(temp_m(N_S),stat=status) !species temp vector
      call check_allocation(status,debug,printer)
      allocate(c(N_R),stat=status) !resources
      call check_allocation(status,debug,printer)
      allocate(r(N_R),stat=status) !output of monod
      call check_allocation(status,debug,printer)
      allocate(temp_c(N_R),stat=status) !species temp vector
      call check_allocation(status,debug,printer)
      allocate(phi(N_S,N_R),stat=status) !met startegies
      call check_allocation(status,debug,printer)
      allocate(phi_T(N_R,N_S),stat=status) !met startegies
      call check_allocation(status,debug,printer)
      allocate(temp_phi(N_S,N_R),stat=status) !met startegies
      call check_allocation(status,debug,printer)
     
      
      !Read stat variables
      open(10,file="temp/stat_var.txt",status='unknown')
      read(10,*) k_n,k_t,gamma_T,rho,eta_T,xi,K,PHI_big,q
      close(10)

      !Read s
      open(10,file="temp/s.txt", status='unknown')
      read(10,*) s
      close(10)

      !Read tau
      open(10,file="temp/tau.txt", status='unknown')
      read(10,*) tau
      close(10)

      !Read dynamical variables at time 0
      open(10,file="temp/dyn_var.txt",status='unknown')
      read(10,*) m,c,phi_T
      close(10)

      phi=transpose(phi_T)
      call tens_product(eta,1.0/rho,N_S,k_n,N_R)
      call tens_product(gamma,1.0/k_t,N_S,k_n,N_R)

      
!!!!!!!! Time evolution: constant proteome fractions !!!!!!!!    
    
      !Open files to save time evolution
      open(11,file="data/species.txt",status='unknown')
      open(12,file="data/resources.txt",status='unknown')
      open(13,file="data/proteome.txt", status='unknown')
      open(14,file="data/square_diff.txt",status='unknown')
      open(15,file="data/time.txt",status='unknown')
      open(16,file="data/m_dot.txt",status='unknown')
      
      do ii=0,N_T-1
         tt=ii*dt
         
         if(mod(ii,int(N_T/10000))==0) then
            write (11,*) m
            write(12,*) c
            write(13,*) transpose(phi)
            write(14,*) square_diff
            write(15,*) tt
            !write(16,*) temp_m
         endif
         
         !Time evolution
         call monod(r,c,K, N_R)
         call check_proteome(phi,PHI_big,r,gamma,N_S,N_R,threshold,square_diff,count)
         call Euler(temp_m,temp_c,temp_phi,m,c,phi,N_S,N_R,tt,dt,eta,gamma,r,s,q,xi,K,tau,count)
      enddo

      !close(11,12,13,14,15,16)
      deallocate(m,c,phi,k_n,k_t,rho,s,xi,K,tau,PHI_big,q,r,temp_m,temp_c,temp_phi)
         
      print *, "Total number of errors:", count
    end program main
    
    
    
    
    

      
