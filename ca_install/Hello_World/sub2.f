      subroutine hello
c

c     	-- free format adds a space
      write(*,*) 'Hello SubWorld 2a'

c       -- this works ok
      write(*,'(A)') 'Hello SubWorld 2b'

      return
      end
