      program debug_program
      use debugging
      implicit none
      logical :: debug,printer
      call init_debug(debug,printer) !Initialize debug, printer
      open(10,file="temp/debugging.txt",status='unknown')
      write(10,*) debug,printer
      close(10)
      end
