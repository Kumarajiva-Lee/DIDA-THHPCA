import numpy as np

def cal_rmse(val0, val1):
    rmse = np.sqrt( np.mean( np.square( val0 - val1 ) ) ) #算rmse
    return rmse

def cal_mean_err(val0, val1):
    mean_err = np.square( np.mean(val0) - np.mean(val1) ) #算mean_err
    return mean_err

def cal_var_rlt(val0, val1):
    var_rlt = np.mean( np.square( val0 - np.mean(val0) ) ) \
                   -2*np.mean( ( val0 - np.mean(val0) ) * ( val1 - np.mean(val1) )) \
                   +np.mean( np.square( val1 - np.mean(val1) ) ) #算var_rlt
    return var_rlt

def cal_theta(bg, an, tr, obs, idx, NumXGrid, NumYGrid, nobs):
    bg_1d = bg.reshape([NumXGrid*NumYGrid, 1], order='F')  #背景场写成一列的变量
    an_1d = an.reshape([NumXGrid*NumYGrid, 1], order='F') 
    tr_1d = tr.reshape([NumXGrid*NumYGrid, 1], order='F')
    
    H_bg_1d = np.zeros(nobs) #设定观测空间
    H_an_1d = np.zeros(nobs)
    H_tr_1d = np.zeros(nobs)

    for i in range(nobs):
        H_bg_1d[i] = bg_1d[int(idx[i])] #把观测空间的模式变量提取出来
        H_an_1d[i] = an_1d[int(idx[i])]
        H_tr_1d[i] = tr_1d[int(idx[i])]

    dis_ob = obs - H_bg_1d  #在观测空间算omb,oma,amb
    dis_oa = obs - H_an_1d
    dis_ab = dis_ob - dis_oa

    osd = np.sqrt(np.dot(dis_ob, dis_oa)/nobs)
    bsd = np.sqrt(np.dot(dis_ab, dis_ob)/nobs)
    asd = np.sqrt(np.dot(dis_ab, dis_oa)/nobs) #算膨胀三角形

    asd_osd = np.arcsin(asd/osd)*180/np.pi
    asd_bsd = np.arccos(asd/bsd)*180/np.pi

    err_o = obs - H_tr_1d      #在观测空间算：观测误差、背景误差、分析误差
    err_b = H_bg_1d - H_tr_1d
    err_a = H_an_1d - H_tr_1d

    ormse = np.sqrt(np.dot(err_o, err_o)/nobs)  #在观测空间算：观测误差、背景误差、分析误差的rmse
    brmse = np.sqrt(np.dot(err_b, err_b)/nobs)
    armse = np.sqrt(np.dot(err_a, err_a)/nobs)
    return osd, asd, bsd, asd_osd, asd_bsd, ormse, brmse, armse