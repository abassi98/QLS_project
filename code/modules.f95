
module debugging
         contains
         !Every subroutine has an error counter 
        
         !initialize a debug variable. If it is true, then the debu gqill be performed
         subroutine init_debug(debug,printer)
            implicit none
            character :: yes, no, var 
            logical :: debug,printer !our debugging variable and printer
            printer=.false. !Default no warmings messages
            yes="y"
            print *, "Do you want to start the debugging mode?[y]/[n]"
            read *, var
            if (var.eq.yes) then
               debug=.true.
               print *, "=================="
               print *, "Debugging mode: ON"
               print *, "=================="
               print *, "Do you want  to print checkpoints?[y]/[n]"
               read *, var
               if(var.eq.yes) printer=.true.
            else
               debug=.false.
               print *, "==================="
               print *, "Debugging mode: OFF"
               print *, "==================="
            endif
            
          end subroutine init_debug
          

      !check if the dimensiions are positive
         subroutine check_dimensions(nn,mm,debug,printer,count)
            implicit none
            integer*16 :: nn,mm !dimensions of the matrix
            logical :: debug,printer
            integer*16 :: count !error counter
            if(debug.eqv..true.) then
               if (nn.ge.1 .and. mm.ge.1) then
                  if(printer.eqv..true.) then
                  print *,  "=================================="
                  print *,  "Check dimensions: OK"
                  print *,  "=================================="
                  endif
               else
                  if(printer.eqv..true.) then
                  print *, "====================================="
                  print *, "Check dimensions of: ERROR"
                  print *, "====================================="
                  endif
                  nn=1 !set default  values for dimensions if it encounters some error
                  mm=1 !try to interpret dimensions anyway
                  count=count+1
                  stop
               endif
            endif
          end subroutine check_dimensions
          
      
      !check if the input is integer and convert it to integer
      subroutine check_rtint(nn,mm,debug,printer,count)
      implicit none
      integer*16 :: count,mm
      double precision :: nn
      logical :: debug,printer
      if(debug.eqv..true.) then
         if(floor(nn).eq.nn) then
            if(printer.eqv..true.) then
               print *, "Check integerness: SUCCEEDED"
            endif
         else
            if(printer.eqv..true.) then
               print *, "Check integerness: FAILED"
            endif
            count=count+1
         endif
      endif
      mm=int(floor(nn)) !in any case convert nn in integer
      end subroutine check_rtint
    
      
      !Check the allocation status
         subroutine check_allocation(stat,debug,printer)
            implicit none
            integer*16 :: stat     !it takes and integr value. 0 is the allocation succeeded
            logical :: debug,printer
            if(debug.eqv..true.) then
               if(stat.ne.0) then
                  if(printer.eqv..true.) then
                  print *, "Allocation:  FAILED"
                  print *, "Program aborted"
               endif
               stop !if allocation fails, the program must be stopped
               else
                  if(printer.eqv..true.) then
                     print *, "Allocation: SUCCEEDED"
                  endif
               endif
            endif
          end subroutine check_allocation
          

!check eigenvalues  computation status
      subroutine check_eigen(info,debug,printer,count)
      implicit none
      integer :: info
      logical :: debug,printer
      integer*16 :: count !error counter
      if(debug.eqv..true.) then
         if(info.eq.0) then
            if(printer.eqv..true.) then
               print *, "Eigenvalues computation: SUCCEEDED"
            endif
         else
            if(printer.eqv..true.) then
               print *, "Eigenvalues computation: FAILED"
               print *, "INFO=",info
            endif
            count=count+1
         endif
      endif
      end subroutine check_eigen
    
         !This subroutine checks if the computation time exceeds a certain threshold.
         subroutine check_time(start,threshold,debug,printer)
            implicit none
            logical :: debug,printer
            double precision :: start, finish,threshold
            if(debug.eqv..true.) then
               call cpu_time(finish)
               if(finish-start .gt. threshold) then
                  if(printer.eqv..true.) then
                  print *, "Computation time exceeded"
                  print *,"Ran in: ",finish-start
                  print *, "Program aborted"
                  endif
                  stop
               endif
            endif
          end subroutine check_time

         subroutine check_diff(mat1,mat2,nn,mm,threshold,debug,printer,count)
            implicit none
            integer*16 :: nn,mm
            double complex, dimension(nn,mm) :: mat1,mat2
            logical :: debug,printer
            double precision :: threshold,diff
            integer*16 :: ii,jj,temp,count
            temp=0
            if(debug.eqv..true.) then
               do jj=1,mm
                  do ii=1,nn
                     diff=sqrt(dble((mat1(ii,jj)-mat2(ii,jj))*conjg(mat1(ii,jj)-mat2(ii,jj))))
                     if(diff.gt.threshold) then
                        temp=temp+1
                     endif    
                  enddo
               enddo
               if(temp.eq.0) then
                  if(printer.eqv..true.) then
                     print *, "Check equality: OK"
                  endif         
               else
                  if(printer.eqv..true.) then
                     print *, "Check equality: ERROR"
                  endif
                  count=count+1
               endif 
            endif 
          end subroutine check_diff
         
         
        end module debugging
        



      module bacteria
      contains
        subroutine rdn_uniform_vector(vec,dmn,a,b)
          implicit none
          integer :: ii,dmn
          double precision, dimension(dmn):: vec
          double precision :: u,a,b
          do ii=1,dmn
             call random_number(u)
             vec(ii)=(b-a)*u+a
          enddo
        end subroutine rdn_uniform_vector

        subroutine rdn_uniform_matrix(mat,nn,mm,a,b)
          implicit none
          integer::ii,jj,nn,mm
          double precision, dimension(nn,mm)::mat
          double precision :: u,a,b
          do jj=1,mm
             do ii=1,nn
                call random_number(u)
                mat(ii,jj)=(b-a)*u+a
             enddo
          enddo
        end subroutine rdn_uniform_matrix


        subroutine monod(r,c,K, N_R)
          implicit none
          integer*16 :: N_R
          double precision, dimension(N_R) :: r,c,K
          r=c/(c+K)
        end subroutine monod

        !This subroutine computes the tensor product of two vectors
        subroutine tens_product(prod,stateA,dimA,stateB,dimB)
          implicit none
          integer*16 :: dimA,dimB,dim,ii,jj
          double precision, dimension(dimA) :: stateA
          double precision, dimension(dimB) :: stateB
          double precision, dimension(dimA,dimB) :: prod      
          do ii=1,dimA
             do jj=1,dimB
                prod(ii,jj)=stateA(ii)*stateB(jj)
             enddo
          enddo
        end subroutine tens_product
        
        subroutine m_dot(temp_m,m,r,eta,phi,q,N_S,N_R)
          implicit none
          integer*16 :: N_S,N_R,ii,jj
          double precision, dimension(N_S)::temp_m,m,q
          double precision, dimension(N_R)::r
          double precision, dimension(N_S,N_R):: eta, phi
          temp_m=m*(matmul(eta*phi,r)-q)
        end subroutine m_dot

        subroutine c_dot(temp_c,m,r,s,xi,phi,N_S,N_R)
          implicit none
          integer*16::N_S,N_R
          double precision, dimension(N_S)::m
          double precision, dimension(N_R)::temp_c,r,s,xi
          double precision, dimension(N_S,N_R)::phi
          temp_c=s-xi*r*matmul(m,phi)
        end subroutine c_dot

        !!! This subroutines takes the diagonal of a square matrix
        subroutine get_diag(diag,mat,dmn)
          implicit none
          integer*16 :: dmn, ii
          double precision, dimension(dmn) :: diag
          double precision, dimension(dmn,dmn):: mat
          do ii=1,dmn
             diag(ii)=mat(ii,ii)
          enddo
        end subroutine get_diag
        
        !!! This subroutine computes phi_dot, returns temp_phi as phi_dot
        subroutine phi_dot(temp_phi,phi,eta,gamma,r,tau,K,c,temp_c,N_S,N_R,count)
          implicit none
          integer*16 :: N_S,N_R,count
          double precision, dimension(N_S) :: tau
          double precision, dimension(N_R) :: r,K,c,temp_c
          double precision, dimension(N_S,N_R) :: temp_phi,phi,eta,gamma
          double precision, dimension(:,:),allocatable ::v, temp,ones_mat,a,b,c_new
          double precision, dimension(:), allocatable :: ones_vecS,ones_vecR,diag
          allocate(temp(N_S,N_R))
          allocate(v(N_S,N_R))
          allocate(ones_mat(N_S,N_R))
          allocate(a(N_S,N_R))
          allocate(b(N_S,N_R))
          allocate(c_new(N_S,N_R))
          allocate(ones_vecR(N_R))
          allocate(diag(N_S))
          allocate(ones_vecS(N_S))
          ones_vecR=1.0
          ones_mat=1.0
          ones_vecS=1.0
          call tens_product(temp,1.0/tau,N_S,r,N_R)
          temp=eta*temp
          call tens_product(a,ones_vecS,N_S,r,N_R)
          v=ones_mat+gamma*a
          call get_diag(diag,matmul(phi,transpose(v**2)),N_S)
          call tens_product(a,diag,N_S,ones_vecR,N_R)
          a=v/a
          b=temp*v
          call tens_product(c_new,ones_vecS,N_S,(K/(c+K)**2)*temp_c,N_R)
          c_new=gamma*c_new
          call get_diag(diag,matmul(phi,transpose(b+c_new)),N_S)
          call tens_product(b,diag,N_S,ones_vecR,N_R)
          temp_phi=phi*(temp-a*b)
          deallocate(temp,ones_mat,ones_vecR,ones_vecS,v,a,b,c_new,diag)
        end subroutine phi_dot

        subroutine Euler(temp_m,temp_c,temp_phi,m,c,phi,N_S,N_R,tt,dt,eta,gamma,r,s,q,xi,K,tau,count)
          implicit none
          integer*16 :: N_S, N_R,count
          double precision :: dt,tt
          double precision, dimension(N_S) :: temp_m, m, q,tau
          double precision, dimension(N_R) :: temp_c, c,r,s,xi,K
          double precision, dimension(N_S,N_R)  :: temp_phi, phi,eta,gamma
          call monod(r,c,K, N_R)
          call m_dot(temp_m,m,r,eta,phi,q,N_S,N_R)
          call c_dot(temp_c,m,r,s,xi,phi,N_S,N_R)
          call phi_dot(temp_phi,phi,eta,gamma,r,tau,K,c,temp_c,N_S,N_R,count)
          m=m+dt*temp_m
          c=c+dt*temp_c
          phi=phi+dt*temp_phi
        end subroutine Euler

          
        subroutine check_proteome(phi,PHI_big,r,gamma,N_S,N_R,threshold,square_diff,count)
          implicit none
          integer*16 :: N_S,N_R,count,ii,temp_count
          double precision :: threshold, square_diff
          double precision, dimension(N_S):: PHI_big
          double precision, dimension(N_R)::  r
          double precision, dimension(:),allocatable:: ones,square
          double precision, dimension(N_S,N_R):: phi,gamma
          allocate(ones(N_R))
          ones=1.0
          allocate(square(N_S))
          square=(matmul(phi,ones)+matmul(phi*gamma,r)-PHI_big)**2
          square_diff=0
          temp_count=0
          do ii=1,N_S
             if(square(ii).gt.threshold) then
                temp_count=temp_count+1
             endif
             square_diff=square_diff+square(ii)
          enddo
          square_diff=square_diff/N_S
          if(temp_count>0) then
             count=count+1
          endif
          deallocate(ones)
        end subroutine check_proteome
        
        subroutine mat_print(matrix, nn,mm)
          implicit none
          integer*16 :: nn,mm,ii,jj
          double precision, dimension(nn,mm) :: matrix
          do ii=1,nn
             print *, (matrix(ii,jj), jj=1,mm)
          enddo
        end subroutine mat_print
        


      end module bacteria
      
  
        
        
