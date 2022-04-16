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

def cal_theta(bg_1d, an_1d, tr_1d, ob_1d):

    dis_ob = ob_1d - bg_1d  #在观测空间算omb,oma,amb
    dis_oa = ob_1d - an_1d
    dis_ab = dis_ob - dis_oa

    nobs = len(ob_1d)

    osd = np.sqrt((np.dot(dis_ob, dis_oa)/nobs))
    bsd = np.sqrt((np.dot(dis_ab, dis_ob)/nobs))
    asd = np.sqrt((np.dot(dis_ab, dis_oa)/nobs)) #算膨胀三角形 绝对值修改记得问一下

    asd_osd = np.arcsin(asd/osd)*180/np.pi
    asd_bsd = np.arccos(asd/bsd)*180/np.pi

    err_o = ob_1d - tr_1d      #在观测空间算：观测误差、背景误差、分析误差
    err_b = bg_1d - tr_1d
    err_a = an_1d - tr_1d

    ormse = np.sqrt(abs(np.dot(err_o, err_o)/nobs))  #在观测空间算：观测误差、背景误差、分析误差的rmse
    brmse = np.sqrt(abs(np.dot(err_b, err_b)/nobs))
    armse = np.sqrt(abs(np.dot(err_a, err_a)/nobs))
    return osd, asd, bsd, asd_osd, asd_bsd, ormse, brmse, armse
    #return ormse, brmse, armse
