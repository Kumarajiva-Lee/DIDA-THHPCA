import numpy as np
# from mpi4py import MPI
import time
import matplotlib.pyplot as plt
plt.switch_backend('agg')

import config as cf
from sklearn.externals import joblib
# import joblib
import gc
from skimage import data, img_as_float
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.metrics import peak_signal_noise_ratio
from skimage.util import random_noise

import scipy as sp
import scipy.ndimage


cf.init()
def get_path(tar,model = '',lev = 32,div = 1,bias = "",model_save_path = '/share1/liminyan/model'):
    if tar == 'phs':lev = 0
    if model == '0':path = model_save_path + '/'+ model +bias+tar+'_'+str(div) +'_'+str(lev)+'_0'
    elif model == '1' :path = model_save_path + '/'+ model +bias+tar+'_'+str(div) +'_'+str(lev)+'_1'
    else:path = model_save_path + '/'+ model +bias+tar+'_'+str(div) +'_'+str(lev)+''
    return path

def get_tar_lev(model = 'line',div = 1,model_save_path = '.'):
    layer_num = cf.get_layer_num()
    tar_l = cf.get_tar_l()
    ord_ = str(cf.get_ord_())

    l = []
    for layer_rank in range(97):
        if layer_rank == 0:
            tar = tar_l[0]
            lev = 0
        else:
            tar = tar_l[int((layer_rank - 1) / 32 + 1) ]
            lev = (layer_rank - 1) % 32

        path = get_path(tar,ord_,lev,div,model,model_save_path)
        l.append(path)

    return l

def get_tar_with_lev(num):

    if num == 0:
        return 'phs',0
    if num <=32:
        return 'pt', num - 1
    if num <=64:
        return 'u_Agrid', num - 1 - 32
    if num <=96:
        return 'v_Agrid', num - 1 - 32 - 32

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()


tarl = ['64','65','66','67']
x = 0
data_path = '/home/xuewei/ddn/share/svr_data/svr_data_100km_12month_ACGrid/'
saved_model_path = data_path +'model/'
l = get_tar_lev(model = 'svr',div = 2,model_save_path = saved_model_path)
test = 'mpas.360x180.forylm_01month.l32.h0_000001.nc'

import netCDF4 as nc

def drwa_model():
    dataset = nc.Dataset(data_path+test)
    for num in range(97):
        print(num)
        tar,lev = get_tar_with_lev(num)
        if num == 0:
            data = dataset.variables[tar][0]
        else:
            data = dataset.variables[tar][0][lev]

        data = np.flipud(data)
        plt.imshow(data)
        plt.colorbar(fraction=0.015, pad=0.01)
        map_key = str(num)
        if num < 10:
            map_key = '0'+map_key

        plt.savefig('model/'+map_key+'l_model.png')
        plt.close()

def draw_recv():
    recv = []
    result = []
    for x in range(4):
        recv.append(np.load('0_'+tarl[x]+'recive.npy'))
        result.append(np.load('0_'+tarl[x]+'result.npy'))
        shape = recv[x].shape
        recv[x] = recv[x].reshape(180,90,97)
        result[x] = result[x].reshape(180,90,97)

    for key in range(0,97):
        print(key)
        rec = []
        res = []
        for x in range(4):
            recv_phs = recv[x][:,:,key].T
            result_phs = result[x][:,:,key].T
            rec.append(recv_phs)
            res.append(result_phs)

        lr_rec0 = np.hstack((rec[0],rec[2]))
        lr_rec1 = np.hstack((rec[1],rec[3]))
        l = np.vstack((lr_rec0,lr_rec1))
        l = np.flipud(l)

        plt.imshow(l)
        plt.colorbar(fraction=0.015, pad=0.01)
        map_key = str(key)
        if key < 10:
            map_key = '0'+map_key

        plt.savefig('picture/'+map_key+'l_recv.png')
        plt.close()

# draw_recv()


def draw_sub(iter = 0,k = 0):
    recv = []
    result = []
    rec = []
    res = []
    key = k
    for x in range(4):
        recv.append(np.load(str(iter)+'_'+tarl[x]+'recive.npy'))
        result.append(np.load(str(iter)+'_'+tarl[x]+'result.npy'))
        shape = recv[x].shape
        recv[x] = recv[x].reshape(180,90,97)
        result[x] = result[x].reshape(180,90,97)
        recv_phs = recv[x][:,:,key].T
        result_phs = result[x][:,:,key].T
        rec.append(recv_phs)
        res.append(result_phs)

    lr_rec0 = np.hstack((rec[0],rec[2]))
    lr_rec1 = np.hstack((rec[1],rec[3]))
    l = np.vstack((lr_rec0,lr_rec1))
    l = np.flipud(l)

    plt.imshow(l)
    plt.colorbar(fraction=0.015, pad=0.01)
    plt.savefig('check/'+tarl[x]+'r_recv.png')
    plt.close()

    lr_res0 = np.hstack((res[0],res[2]))
    lr_res1 = np.hstack((res[1],res[3]))
    r = np.vstack((lr_res0,lr_res1))
    r = np.flipud(r)

    plt.imshow(r*0.3 + 0.7*l)
    plt.colorbar(fraction=0.015, pad=0.01)
    plt.savefig('check/'+tarl[x]+'r_result.png')
    plt.close()

    plt.imshow(r)
    plt.colorbar(fraction=0.015, pad=0.01)
    plt.savefig('check/'+tarl[x]+'r_org.png')
    plt.close()


    merge =  np.vstack((l,r*0.3 + 0.7*l))
    plt.imshow(merge)
    plt.colorbar(fraction=0.015, pad=0.01)
    plt.savefig('check/'+tarl[x]+'merge.png')
    plt.close()



def check_model():
    dataset = nc.Dataset(data_path+test)
    all_data = []
    for num in range(97):
        tar,lev = get_tar_with_lev(num)
        if num == 0:
            data = dataset.variables[tar][0]
        else:
            data = dataset.variables[tar][0][lev]
        # data = np.flipud(data)
        # data = data.T

        all_data.append(data.reshape(180*360))
    all_data = np.array(all_data)
    path = saved_model_path
    l_model = get_tar_lev('line',4,path)
    print(all_data.shape)
    lev = l_model[0]
    res = []
    for i in range(180 * 360):
        x = int(i / 360)
        y = int(i % 360)
        key = str(x)+'_'+str(y)
        file = lev + '/' + key +'.m'
        model = joblib.load(file)
        print(file)
        res.append(model.predict([all_data[:,i]]))
        # print(model)
        print(all_data[:,i])
        print(res[0])
        return

        if x % 10 == 0 and y == 0:
            print(x,'/ 180')
            del model
            gc.collect()

    res = np.array(res)
    res = res.reshape((180,360))
    res = np.flipud(res)

    np.save('np_res.npy', x)
    plt.imshow(res)
    plt.colorbar(fraction=0.015, pad=0.01)
    plt.savefig('check/model_res.png')
    plt.close()

    # print(all_data.shape)

# check_model()

def check_dida():
    recv = []
    result = []
    for x in range(4):
        recv.append(np.load('0_'+tarl[x]+'recive.npy'))
        result.append(np.load('0_'+tarl[x]+'result.npy'))
        shape = recv[x].shape
        recv[x] = recv[x].reshape(180,90,97)
        result[x] = result[x].reshape(180,90,97)

    all_data = []
    for key in range(0,97):
        print(key)
        rec = []
        res = []
        for x in range(4):
            recv_phs = recv[x][:,:,key].T
            result_phs = result[x][:,:,key].T
            rec.append(recv_phs)
            res.append(result_phs)

        lr_rec0 = np.hstack((rec[0],rec[2]))
        lr_rec1 = np.hstack((rec[1],rec[3]))
        l = np.vstack((lr_rec0,lr_rec1))
        all_data.append(l)

    all_data = np.array(all_data)
    pr_res = []
    path = saved_model_path
    l_model = get_tar_lev('line',4,path)
    print(all_data.shape)
    lev = l_model[0]

    for i in range(180 * 360):

        x = int(i / 360)
        y = int(i % 360)
    # for x in range(180):
        # print(x,'/ 180')
        # for y in range(360):
        key = str(x)+'_'+str(y)
        file = lev + '/' + key +'.m'
        model = joblib.load(file)
        pr_res.append(model.predict([all_data[:,x,y]]))
    if x % 10 == 0:
        del model
        gc.collect()


    pr_res = np.array(pr_res)
    pr_res = pr_res.reshape((180,360))
    pr_res = np.flipud(pr_res)
    np.save('np_recv.npy', x)
    plt.imshow(pr_res)
    plt.colorbar(fraction=0.015, pad=0.01)
    plt.savefig('check/dida_res.png')
    plt.close()


def check_err_model():
    model = 'dline'
    import os
    path = saved_model_path
    l_model = get_tar_lev(model,4,path)
    print(len(l_model))

    for lev in range(2,8):
        l_model = get_tar_lev(model,lev*2,path)
        for i in l_model:
            if os.path.exists(i):
                file = os.listdir(i)
                if len(file) != 64800:
                    print('err',i)
            else:
                print('err',i)





def denoise(data):

    noisy = data
    sigma_est = np.mean(estimate_sigma(noisy))
    patch_kw = dict(patch_size=5,      # 5x5 patches
                    patch_distance=6,  # 13x13 search area
                    )
    denoise_fast = denoise_nl_means(noisy, h=0.8 * sigma_est, fast_mode=True,**patch_kw)
    return denoise_fast


def denoise_gass(data):

    
    noisy = data
    sigma = [5, 3]
    denoise_fast = sp.ndimage.filters.gaussian_filter(noisy, sigma, mode='constant')

    return denoise_fast


def convolution2d(image, kernel):
    m, n = kernel.shape
    if (m == n):
        y, x = image.shape
        y = y - m + 1
        x = x - m + 1
        new_image = np.zeros((y,x))
        for i in range(y):
            for j in range(x):
                new_image[i][j] = np.sum(image[i:i+m, j:j+m]*kernel) 
    return new_image

key = 60
layer = 0

# for key in range(97):
#     time.sleep(6)
# for div in range(1,4):
recv = []
result = []

# data_sz = np.load('data_sz.npy')
# print(data_sz)

for x in range(4):
    # recv.append(np.load(str(layer)+'_'+tarl[x]+'send.npy'))
    # shape = recv[x].shape
    # recv[x] = recv[x].reshape(180,90,97,4)
    tmp = np.load(str(layer)+'_'+tarl[x]+'result.npy')[2]
    result.append(tmp)
    print('tmp',tmp.shape)

    # print(result[x].shape)
    # result[x] = result[x].reshape(190,90,97)
    # result[x] = result[x].transpose(2,0,1)

rec = []
res = []
for x in range(4):
    res.append(result[x][key])

lr_res0 = np.hstack((res[0],res[2]))
lr_res1 = np.hstack((res[1],res[3]))
r = np.vstack((lr_res0,lr_res1))
r = np.flipud(r)
print(r.shape)


# denoise_fast = sp.ndimage.filters.gaussian_filter(r, [1,1], mode='reflect')
# denoise_fast = denoise_gass(r)
# print(np.mean((r - denoise_fast)**2))
# r = denoise_fast - r
# r = np.c_[r,denoise_fast]

plt.imshow(r)


# plt.imshow(r,vmin=310, vmax=350)
plt.colorbar(fraction=0.015, pad=0.01)
plt.savefig('check/'+tarl[x]+'r_recv_check.png')
plt.close()

