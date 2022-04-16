module layer
    use m_MCTWorld,only: MCTWorld_init => init
    use m_MCTWorld,only: MCTWorld_clean => clean
    use m_GlobalSegMap,only: GlobalSegMap
    use m_GlobalSegMap,only: GlobalSegMap_init => init
    use m_GlobalSegMap,only: GlobalSegMap_lsize => lsize
    use m_GlobalSegMap,only: GlobalSegMap_clean => clean
    use m_GlobalSegMap,only: GlobalSegMap_Ordpnts => OrderedPoints
    use m_Rearranger,only: Rearranger
    use m_Rearranger,only: Rearranger_init => init
    use m_Rearranger,only: Rearranger_clean => clean
    use m_Rearranger,only: Rearranger_rearrange => rearrange  
    use m_AttrVect,only    : AttrVect
    use m_AttrVect,only    : AttrVect_init => init
    use m_AttrVect,only    : AttrVect_clean => clean
    use m_AttrVect,only    : AttrVect_indxR => indexRA
    use m_AttrVect,only    : AttrVect_importRAttr => importRAttr
    use m_AttrVect,only    : AttrVect_exportRAttr => exportRAttr
    use m_Router,only: Router
    use m_Router,only: Router_init => init
    use m_Router,only: Router_clean => clean
    use m_Transfer,only : MCT_Send => send
    use m_Transfer,only : MCT_Recv => recv
  
    use coupler_config
    use endnn_config

    integer :: layer_id
    integer :: layer_comm
    integer :: lon_num, lat_num ! 本层，本进程处理的网格数
    integer :: upper_lon_num, upper_lon_lat !本进程处理的上层网格数
    integer :: global_lon, global_lat ! 本层的总网格数
    integer :: upper_global_lon, upper_global_lat ! 上层的总网格数
    integer :: lon_ibeg, lon_iend, lat_ibeg, lat_iend ! 本层，本进程的网格范围，有重叠
    integer :: upper_lon_ibeg, upper_lon_iend, upper_lat_ibeg, upper_lat_iend ! 本进程处理的，上一层的网格范围，无重叠

    integer :: x_sec, x_res, y_sec, y_res

    type(GlobalSegMap) :: Map_up, Map_down_overlap, Map_down_no_overlap
    type(Router) :: Router_up, Router_down
    type(AttrVect) :: AV_up, AV_down_overlap, AV_down_no_overlap
    type(Rearranger) :: Re

    ! real, allocatable, dimension(:,:) :: my_grid !本层网格的数据
    ! real, allocatable, dimension(:,:) :: upper_grid !卷积后要向下层网格传输的数据

    integer :: comm_cart
    integer :: dims(2), coords(2)
    logical :: periods(2)  
    integer :: rank, lrank, lsize ! lsize 本层的总进程数

    integer, allocatable, dimension(:) :: my_start, my_length, upper_start, upper_length

    contains
    
    ! layer_num : endnn的总层数, 直接从config中获取 

    ! world_comm : MPI_DA_GROUP，DA的总通信子
    ! my_comm : 该层内部通信使用的MPI通信子
    ! my_layer_id : 该层的编号(编号从0~layer_num-1，默认最底层为0)
    ! my_mct_id : 该层MCT通信部件在所有部件中的id
    ! my_x, my_y, my_lev : 本层的网格维度
    ! next_x, next_y, next_lev : 上一层的网格维度

    subroutine layer_init(world_comm, my_comm, my_layer_id, my_mct_id, &
        my_x, my_y, my_lev, next_x, next_y, next_lev, var_num)
        integer, intent(in) :: world_comm, my_comm
        integer, intent(in) :: my_layer_id, my_mct_id
        integer, intent(in) :: my_x, my_y
        integer, intent(in) :: next_x, next_y
        integer :: proc_num, ierrs
        integer :: i

        layer_id = my_layer_id
        layer_comm = my_comm

        !!! MPI_INIT
        !! layer_0 考虑直接调用da_parallel_init
        ! parallel init from my_x, my_y, next_x, next_y

        dims(1:2) = 0.0
        call MPI_Comm_rank(my_comm, rank, ierrs)
        call MPI_Comm_size(my_comm, lsize, ierrs)
        call MPI_Dims_create(lsize, 2, dims,ierrs)
        call MPI_Cart_create(MPI_DA_GROUP, 2, dims, periods, 0, comm_cart, ierrs)
        call MPI_Comm_rank(comm_cart, lrank, ierrs)
        call MPI_Cart_coords(comm_cart, lrank, 2, coords, ierrs)

        global_lon = my_x
        global_lat = my_y
        upper_global_lon = next_x
        upper_global_lat = next_y


        ! 进程划分方式每层网格都相同
        x_sec = upper_global_lon / dims(1)
        x_res = mod(upper_global_lon, dims(1))
        y_sec = upper_global_lat / dims(2)
        y_res = mod(upper_global_lat, dims(2))

        if(coords(1) < x_res) then
            upper_lon_ibeg = coords(1) * (x_sec + 1) + 1
            upper_lon_iend = upper_lon_ibeg + x_sec
        else
            upper_lon_ibeg = coords(1) * x_sec + x_res + 1
            upper_lon_iend = upper_lon_ibeg + x_sec - 1
        end if
    
        if(coords(2) < y_res) then
            upper_lat_ibeg = coords(2) * (y_sec + 1) + 1
            upper_lat_iend = upper_lat_ibeg + y_sec 
        else
            upper_lat_ibeg = coords(2) * y_sec + y_res + 1
            upper_lat_iend = upper_lat_ibeg + y_sec - 1
        end if

        upper_lon_num = upper_lon_iend - upper_lon_ibeg + 1
        upper_lat_num = upper_lat_iend - upper_lat_ibeg + 1
        
        ! [(upper_xst, upper_yst), (upper_xed, upper_yed)] ->
        !      [(stride*upper_xst, stride*upper_yst), (stride*upper_xed+kx, stride*upper_yed+ky)]
        lon_ibeg = stride_x * upper_lon_ibeg
        lon_iend = stride_x * upper_lon_iend + kernel_x - 1
        lat_ibeg = stride_y * upper_lat_ibeg
        lat_iend = stride_y * upper_lat_iend + kernel_y - 1
        
        lon_num = lon_iend - lon_ibeg + 1
        lat_num = lat_iend - lat_ibeg + 1

        allocate(my_grid (lon_ibeg:lon_iend,lat_ibeg:lat_iend))

        if (layer_id /= layer_num - 1) then ! 非最上层
            allocate(upper_grid (upper_lon_ibeg:upper_lon_iend,upper_lat_ibeg:upper_lat_iend))
        end if
        

        !!! MCT_INIT
        
        call MCTWorld_init(layer_num, MPI_COMM_WORLD, my_comm, dnn_layer_id(my_id)) 
        ! 第0层要对外通信，所以只能用MPI_COMM_WORLD，上层应该都可以用，再议

        allocate (my_start(lon_ibeg:lon_iend))
        allocate (my_length(lon_ibeg:lon_iend))

        do i = lon_ibeg, lon_iend
            
        end do

        ! Map_down_init from my_x,my_y
        ! Map_up_init from upper_x, upper_y
        ! Router_up/down_init


    end subroutine layer_init


    subroutine layer_work(op_type)
        integer, intent(in) :: op_type

        ! call region_init(range, my_grid, upper_grid)

    end subroutine layer_work

    subroutine layer_fit()
    end subroutine

    subroutine layer_inference()
    end subroutine
    
    subroutine layer_sample()
    end subroutine

    subroutine layer_clean()
        
        if (allocated(my_grid)) deallocate(my_grid)
        if (allocated(upper_grid)) deallocate(upper_grid)

    end subroutine layer_clean
    

end module layer