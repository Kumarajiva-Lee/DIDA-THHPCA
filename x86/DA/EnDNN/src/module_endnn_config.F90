module endnn_config

    implicit none

    integer :: endnn_switch ! 是否使用endnn，1为是，0为否，否时逻辑应和letkf一样
    
    integer :: kernel_x, kernel_y
    integer :: stride_x, stride_y
    integer :: stride, kernel

    integer :: layer_num

    integer, dimension(:), allocatable  :: dnn_layer_id

    integer :: OP_Fit, OP_Inference, OP_Sample


end module endnn_config