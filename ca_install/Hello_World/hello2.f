      program hello2
c

c     	-- free format adds a space
      write(*,*) 'Hello World 2a'

c       -- this works ok
      write(*,'(A)') 'Hello World 2b'

      end
