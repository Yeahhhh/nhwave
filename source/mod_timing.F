! http://people.sc.fsu.edu/~jburkardt/f_src/timer/timer.html

module timing
! use omp_lib
        parameter (NTIMER = 8)
        real(kind(0.0d0)) :: tt(NTIMER)
end module


subroutine TimeNow (t)
        implicit none
        real(kind(0.0d0)) :: t
        integer tt, clock_rate, clock_max

        call system_clock ( tt, clock_rate, clock_max )
        t = real ( tt ) / real ( clock_rate )
        ! t = omp_get_wtime ( )
        ! call cpu_time (t)
end subroutine


subroutine ZeroTimers
        use timing
        integer i

        do i = 1, NTIMER
           tt(i) = 0d0
        enddo
end subroutine 


subroutine PrintTimers
        use timing
        integer i

        do i = 1, NTIMER
           write(6,*) "timer", i, tt(i)
        enddo
end subroutine


