module UseSW_mod

    USE block_mod
    USE mesh_mod

    implicit none

    external :: slave_init, slave_prepare

    type var_prepare
        integer*8 :: vdelon
        integer*8 :: vdelat
        integer*8 :: vlelon
        integer*8 :: vlelat
        integer*8 :: varea_cell
        integer :: vsize
    end type var_prepare


    contains

    subroutine sw_init(block)
        
        type(block_type), intent(inout), target :: block

        type(var_prepare) :: vp

        call slave_init()

        vp%vdelon      = loc(block%mesh%de_lon)
        vp%vdelat      = loc(block%mesh%de_lat)
        vp%vlelon      = loc(block%mesh%le_lon)
        vp%vlelat      = loc(block%mesh%le_lat)
        vp%varea_cell  = loc(block%mesh%area_cell)

        vp%vsize = size(block%mesh%de_lon) * 8

        call athread_spawn(slave_prepare, vp)
        call athread_join()

    end subroutine


end module UseSW_mod