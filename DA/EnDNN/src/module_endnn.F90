module endnn
    use mpi
    use layer
    use endnn_config
    use coupler_config
    use letkf_mod, only:letkf_init, letkf_run
    !use toy_model
    use pro_info
    use string
    use splitchar
    use var_info_mod
    use obsmaker_mod
    use da_namelist_mod, only:obs_hrz_interval_lon, &
    obs_hrz_interval_lat, obs_vtc_interval, parse_namelist, &
    da_halo
    
    integer, dimension(:), allocatable      :: color
    integer :: local_color
    integer :: new_comm
    integer :: proc_num
    integer, dimension(:), allocatable :: layer_x, layer_y, layer_lev !每层的网格维度，垂直层数
    integer, dimension(:), allocatable :: layer_proc !每层的进程数量
    integer :: var_num

    real, dimension(:,:,:), allocatable :: my_grid, upper_grid
    !! 数组均使用指针的方式传递
    !! 变量在endnn层定义，在layer层进行内存的分配，内部不需要再做处理
    !! 变量的维度顺序 x,y,lev,var (ens? gmcore传的数据带集合，endnn对带集合数据的处理?)
    !! var层在最内层还是

    contains

    subroutine endnn_init(namelist_path, MPI_DA_GROUP)
        character(*), intent(in) :: namelist_path
        integer, intent(in) :: MPI_DA_GROUP
        integer :: ierrs
        integer :: i, j, myrank
        integer :: layer_num_x, layer_num_y
        integer :: layer_sum, left_proc, left_layer, temp

        call MPI_Comm_size(MPI_DA_GROUP, proc_num, ierrs)
        call MPI_Comm_rank(MPI_DA_GROUP, myrank, ierrs)

        !!!  从namelist中读取参数，按照wiki上的公式，推算stride以及网格层数
        call parse_namelist(namelist_path)
        !从namelist中获取：kernel, 底层网格维度，垂直方向层数，变量个数(b)
        ! warning: namelist中需要增加的dnn相关参数：kernel,stride,dnn_layer_id,PCA主分量数m
        ! 目前暂定变量名为kernel,stride,dnn_layer_id,pca_m
        ! endnn_switch应加入在总的namelist里
        ! 观测场方面的设计待考虑

        ! stride = floor((sqrt(pca_m) - 1)*kernel/(kernel-1)) 
        ! stride也直接从配置里读 

        if (endnn_switch == 1) then

            layer_num_x = floor(log(num_lat) / log(stride)) !取下整
            layer_num_y = floor(log(num_lon) / log(stride)) 
            
            if (layer_num_x < layer_num_y) then
                layer_num = layer_num_x
            else 
                layer_num = layer_num_y
            end if

            allocate(layer_x (0:layer_num))
            allocate(layer_y (0:layer_num))
            allocate(layer_lev (0:layer_num))
            allocate(layer_proc (0:layer_num))
            ! 实际网格层的下标为 0 ~ layer_num - 1 , layer_num项为方便最顶层的网格划分用

            !！ 根据kernel和stride，计算各层网络中的网格规模

            layer_x(0) = num_lon
            layer_y(0) = num_lat
            layer_lev(0) = num_lev

            do i = 1, layer_num
                layer_x(i) = nint((layer_x(i-1) - kernel) / stride) + 1 !默认是整除的，保险起见用四舍五入
                layer_y(i) = nint((layer_y(i-1) - kernel) / stride) + 1
                layer_lev(i) = 1
            end do

            !! 根据各层的网格规模，进行进程和通信子划分
            layer_sum = 0
            do i = 0, layer_num
                layer_sum = layer_sum + layer_x(i) * layer_y(i)
            end do

            left_layer = layer_num 
            left_proc = proc_num
            do i = 0, layer_num - 1
                layer_proc(i) = floor(proc_num * layer_x(i) * layer_y(i) * / layer_sum)  !lev/ens/var
                if (layer_proc(i) == 0) then
                    layer_proc(i) = 1
                end if
                left_layer = left_layer - 1
                left_proc = left_proc - layer_proc(i)
                if (left_proc < left_layer) then ! 补正，防止某层网格分到0个进程
                    layer_proc(i) = layer_proc(i) - (lefy_layer - left_proc)
                end if
            end do

        else ! 只有底层网络和letkf
            
            layer_num = 1
            var_num = 4

            allocate(layer_x (0:layer_num))
            allocate(layer_y (0:layer_num))
            allocate(layer_lev (0:layer_num))
            allocate(layer_proc (0:layer_num))

            layer_x(0) = num_lon
            layer_y(0) = num_lat
            layer_lev(0) = num_lev
            layer_proc(0) = proc_num

            layer_x(1) = 0
            layer_y(1) = 0
            layer_lev(1) = 0
            layer_proc(1) = 0

        end if
    !    layer_proc(0) = layer_proc(0) + left_proc !取下整，所以可能有没分完的
        temp = 0

        allocate(color (0:proc_num-1))
        do i = 0, layer_num - 1
            do j = 1, layer_proc(i)
                color(temp) = i
                temp = temp + 1
            end do
        end do
        
        local_color = color(myrank)

        ! dnn_layer_id, 每层作为MCT通信组件的Id
        allocate(dnn_layer_id(1:layer_num))
        do i = 1, layer_num
            dnn_layer_id(i) = DEFAULT_Id + 1 + i
        end do

        call MPI_COMM_SPLIT(MPI_DA_GROUP, local_color, 0, new_comm, ierr)

        call layer_init(MPI_DA_GROUP, new_comm, local_color, dnn_layer_id(local_color), &
        layer_x(local_color), layer_y(local_color), layer_lev(local_color), &
        layer_x(local_color+1), layer_y(local_color+1), layer_lev(local_color+1), var_num)

    end subroutine

    subroutine endnn_run()

        integer :: is_da_time
        integer :: da_time

        da_time = 0
        !time schedule
        do while (1)
            if (is_da_time) then 

                call da_run()
                da_time = da_time + 1
            end if
            call self_training(my_grid)
            if (da_time == cal_da_time) then
                exit
            end if
        end do
    end subroutine

    !! fit/inference/sample参数类型不同的问题，拆解bottom_up和top_down
    subroutine self_training()
        call bottom_up_fit()
        call top_down_fit()
    end subroutine

    subroutine da_run()
        call bottom_up_inf()
        call top_down_sample()
    end subroutine

    subroutine bottom_up(op_type)
        integer, intent(in) :: op_type ! fit or inference
        if (local_color == 0) then
            call layer_Work(op_type, my_grid, upper_grid)
            call layer_Send_Up(upper_grid)
        else if (local_color == layer_num - 1) then
            call layer_Rev_from_bottom(my_grid)
            call layer_Work(op_type, my_grid, upper_grid)
        else 
            call layer_Rev_from_bottom(my_grid)
            call layer_Work(op_type, my_grid, upper_grid)
            call layer_Send_Up(op_type, my_grid, upper_grid)
        end if
    end subroutine

    subroutine bottom_up_fit(fit_param)

        integer, intent(in) :: fit_param ! fit

        !! fit/inference/sample参数类型不同的问题，拆解bottom_up和top_down

        if (local_color == 0) then
            call layer_fit(op_type, my_grid, upper_grid)
            call layer_Send_Up(upper_grid)
        else if (local_color == layer_num - 1) then
            call layer_Rev_from_bottom(my_grid)
            call layer_fit(op_type, my_grid, upper_grid)
        else 
            call layer_Rev_from_bottom(my_grid)
            call layer_fit(op_type, my_grid, upper_grid)
            call layer_Send_Up(op_type, my_grid, upper_grid)
        end if
    end subroutine

    subroutine top_down_fit(fit_param)
        integer, intent(in) :: fit_param

        if (local_color == 0) then
            call layer_Recv_from_top(my_grid, upper_grid)
            call layer_fit()
        else if (local_color == layer_num - 1) then
            call layer_fit()
            call layer_Send_Down(my_grid)
        else 
            call layer_Recv_from_top(my_grid, upper_grid)
            call layer_fit()
            call layer_Send_Down(my_grid)
        end if
    end subroutine

end module