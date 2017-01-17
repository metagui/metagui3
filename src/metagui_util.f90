program metagui_util
   character*20 task
   character*100 inputfile
   call getarg(1,task)
   call getarg(2,inputfile)
   select case (trim(task))
     case ("MICROSTATES")
       call make_microstates(inputfile)
     case ("WHAM")
       call wham_on_microstates(inputfile)
     case("CLUSTERS")
       call smart_cluster(inputfile)
     case("KINETIC")
       call kinetic_basins(inputfile)
     case default
       write (6,*) "CASE ",trim(task),"unknown"
       write (6,*) "STOP"
       STOP
   end select
end program metagui_util
