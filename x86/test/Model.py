# from sklearn.svm import SVR,LinearSVR
# import torch 4
# import torch.nn as nn
import time
import numpy as np
import multiprocessing
import time
import mpi4py.MPI as MPI
import os

import joblib
# from sklearn import preprocessing

# from sklearn.linear_model import SGDRegressor
# from sklearn.pipeline import make_pipeline
# from sklearn.preprocessing import StandardScaler
# from sklearn.kernel_approximation import Nystroem

# import scipy as sp
# import scipy.ndimage

# from skimage.restoration import denoise_nl_means, estimate_sigma
# 
import torch
import torch.nn as nn
import torch.utils.data as Data
import copy

# from sklearn.neural_network import MLPRegressor
# from sklearn.datasets import make_regression
# from sklearn.model_selection import train_test_split
import config as cf
from pathlib import Path
class Model(object):
    """
    <<interface>> Model，接口类
    """

    def fit(self):
        """
        模型训练
        """
        pass

    def predict(self):
        """
        模型预测
        """
        pass

    def save(self, path):
        """
        模型保存
        """
        pass

    def load(self, path):
        """
        模型加载
        """
        pass

def get_tar_lev(model = 'svr',div = 1,model_save_path = '.'):
    layer_num = cf.get_layer_num()
    tar_l = cf.get_tar_l()
    ord_ = str(cf.get_ord_())

    l = []
    for layer_rank in range(layer_num):
        if layer_rank == 0:
            tar = tar_l[0]
            lev = 0
        else:
            tar = tar_l[int((layer_rank - 1) / 32 + 1) ]
            lev = (layer_rank - 1) % 32

        path = get_path(tar,ord_,lev,div,model,model_save_path)
        l.append(path)

    return l


def get_path(tar,model = '',lev = 32,div = 1,bias = "",model_save_path = '/share1/liminyan/model'):
    if tar == 'phs':lev = 0
    if model == '0':path = model_save_path + '/'+ model +bias+tar+'_'+str(div) +'_'+str(lev)+'_0'
    elif model == '1' :path = model_save_path + '/'+ model +bias+tar+'_'+str(div) +'_'+str(lev)+'_1'
    else:path = model_save_path + '/'+ model +bias+tar+'_'+str(div) +'_'+str(lev)+''
    return path

class ConvDeconvNeuralNetwork(Model):

    """
    用于集合后处理的卷积反卷积网络
    Attributes:
        kernel_model: torch框架下的自建模型
        learning_rate: float 学习率
        weight_decay: L2范数正则项在损失函数中的权重
        optimizer: Adam优化器
        batch_size: int 批样本数量
        epochs: int 所有样本上训练的次数
        patiences: int 早停策略下的容忍度 n/次epochs
    """

    def __init__(self, arch=None, learning_rate = 0.005, weight_decay=0.005,batch_size = 32, epochs = 30, patiences = 100):
        self.kernel_model = arch
        self.learning_rate = learning_rate
        self.weight_decay = weight_decay
        self.optimizer = torch.optim.Adam(self.kernel_model.parameters(), lr=self.learning_rate, weight_decay=self.weight_decay)
        self.batch_size = batch_size
        self.epochs = epochs
        self.patiences = patiences
        self.tar = ''
        self.lev = 0
        self.config = None
        self.div = 1


    def unnormalized(self,img,teg = 'Y'):

        if teg == 'X':
            img = img * self.std_X + self.mu_X
        elif teg == 'Y':
            img = img * self.std_Y + self.mu_Y
        return img

    def predict(self, data):
        """
        模型预测
        :param DataProcessor: 数据处理器实例，主要使用其中的 x_test
        :return: 返回预测值或实例本身
        """
        # 将模型切换为预测状态，将测试集数据转换为 torch 数据进行预测，预测结果再转换为 numpy.ndarray 数据
        self.kernel_model.eval()
        process_data = (data - self.mu_X)/self.std_X
        x_test = torch.from_numpy(process_data).type(torch.FloatTensor)
        y_pred = self.kernel_model(x_test)
        return self.unnormalized(y_pred.detach().numpy())

    def save(self, path):
        """
        模型保存
        :param path: 模型保存的路径，保存的文件为 .pt 格式，如 0302_PRAVG_CNN.pt
        :return: 返回实例本身
        """
        torch.save({
            'model_state_dict': self.kernel_model.state_dict(),
            # 'optimizer_state_dict': self.optimizer.state_dict(),
            'mu_X': self.mu_X,
            'std_X': self.std_X,
            'mu_Y': self.mu_Y,
            'std_Y': self.std_Y
        }, path)



        return self

    def get_item_mu(self,data,ibegin,iend,jbegin,jend,shape = 3,halo = 10):
    




        ibegin += halo -1
        iend += halo 
        ibegin -= halo
        iend += halo

        jbegin += halo -1
        jend += halo 
        jbegin -= halo
        jend += halo

        # print('org',data.shape)
        # print('log',ibegin,iend,jbegin,jend)

        real = data[:,5:-5,5:-5]
        left = real[:,:,-halo:]
        right = real[:,:,0:halo]

        real = np.append(left,real,axis = 2)
        real = np.append(real,right,axis = 2)
        up = np.repeat(real[:,[0]], repeats=halo, axis=1) 
        down = np.repeat(real[:,[-1]], repeats=halo, axis=1)

        # print('left',left.shape)
        # print('right',right.shape)
        # print('up',up.shape)
        # print('down',down.shape)
        real = np.append(up,real,axis=1)
        real = np.append(real,down,axis=1)
        # print('real',real.shape)
        mid = real[:,ibegin:iend,jbegin:jend]

        # print('mid',mid.shape)
        return mid
        # mid = mid.reshape(shape[0],1,shape[1],shape[2])

    def load(self, path,ibegin,iend,jbegin,jend,halo = 5):
        """
        模型加载
        :param path: 加载模型的路径，加载的文件为 .pt 格式，如 0302_PRAVG_CNN.pt
        :return: 返回实例本身
        """
        checkpoint = torch.load(path)
        self.kernel_model.load_state_dict(checkpoint['model_state_dict'])
        self.mu_X = np.load(self.model_save_path + '/mu_X.npy')
        self.std_X = np.load(self.model_save_path + '/std_X.npy')

        
        level = self.lev
        if self.tar =='pt':
            level += 1
        elif self.tar =='u_Agrid':
            level += 1 + 32
        elif self.tar =='v_Agrid':
            level += 1 + 32 + 32
    
        # self.mu_X = checkpoint['mu_X']
        # self.std_X = checkpoint['std_X']
        # self.mu_Y = checkpoint['mu_Y']
        # self.std_Y = checkpoint['std_Y']

        self.kernel_model.eval()

        comm = MPI.COMM_WORLD
        comm_rank = comm.Get_rank()
        comm_size = comm.Get_size()
        # comm_size = 120 
        # max_size = 180 * 360
        # per_size = int(max_size/comm_size)
        # begin = comm_rank*per_size   
        # end = (comm_rank + 1)*per_size 
        # i_begin = int(begin/360)
        # i_end = int(end/360)
        # j_begin = int(begin%360)
        # j_end = int(end%360)

        self.mu_X = self.get_item_mu(data = self.mu_X,ibegin=ibegin,iend = iend,jbegin = jbegin,jend = jend,halo = halo)
        self.std_X = self.get_item_mu(data = self.std_X,ibegin=ibegin,iend = iend,jbegin = jbegin,jend = jend,halo = halo)

        self.mu_Y = self.mu_X[level]
        self.std_Y =self.std_X[level]

        if comm_rank == 70:
            np.save('mux.npy',self.mu_X)
            np.save('std.npy',self.std_X)
       
        # self.mu_Y = self.get_item_mu(self.mu_Y,ibegin,iend,jbegin,jend,shape = 2)
        # self.std_Y = self.get_item_mu(self.std_Y,ibegin,iend,jbegin,jend,shape = 2)

        return self


    @staticmethod
    def weigth_init(m):
        """
        初始化模型参数
        :param m: torch中神经网络层
        """
        if isinstance(m, torch.nn.Conv2d):
            m.weight.data.normal_(0.1)
            m.bias.data.zero_()

class CDNNArch(torch.nn.Module):
    """
    torch架构下的卷积反卷积神经网络架构
    """
    def __init__(self):
        super(CDNNArch, self).__init__()

        self.input = torch.nn.Sequential(
            torch.nn.BatchNorm2d(97)
            )

        self.conv1 = torch.nn.Sequential(
            torch.nn.Conv2d(97, 8, 3),
            torch.nn.BatchNorm2d(8),
            # torch.nn.Tanh()
            torch.nn.Tanh()
            # torch.nn.ReLU()
            )

        self.conv2 = torch.nn.Sequential(
            torch.nn.Conv2d(8, 16, 3),
            # torch.nn.BatchNorm2d(16),
            # torch.nn.Tanh()
            torch.nn.Tanh()
            # torch.nn.ReLU()
            )

        self.conv3 = torch.nn.Sequential(
            torch.nn.Conv2d(16, 32, 3),
            # torch.nn.BatchNorm2d(32),
            # torch.nn.Tanh()
            torch.nn.Tanh()
            # torch.nn.ReLU()
            )

        # self.conv4 = torch.nn.Sequential(
        #     torch.nn.Conv2d(32, 64, 3),
        #     # torch.nn.BatchNorm2d(64),
        #     torch.nn.Tanh()
        #     )

        # self.conv5 = torch.nn.Sequential(
        #     torch.nn.Conv2d(64, 128, 3),
        #     # torch.nn.BatchNorm2d(64),
        #     torch.nn.Tanh()
        #     )

        # self.conv6 = torch.nn.Sequential(
        #     torch.nn.Conv2d(128, 256, 3),
        #     # torch.nn.BatchNorm2d(64),
        #     torch.nn.Tanh()
        #     )


        # self.conv6_merge = torch.nn.Sequential(
        #     torch.nn.Conv2d(256, 128, 3,padding=1),
        #     # torch.nn.BatchNorm2d(32),
        #     torch.nn.Tanh()
        #     )


        # self.conv5_merge = torch.nn.Sequential(
        #     torch.nn.Conv2d(128, 64, 3,padding=1),
        #     # torch.nn.BatchNorm2d(32),
        #     torch.nn.Tanh()
        #     )


        # self.conv4_merge = torch.nn.Sequential(
        #     torch.nn.Conv2d(64, 32, 3,padding=1),
        #     # torch.nn.BatchNorm2d(32),
        #     torch.nn.Tanh()
        #     )



        self.conv3_merge = torch.nn.Sequential(
            torch.nn.Conv2d(32, 16, 3,padding=1),
            # torch.nn.BatchNorm2d(16),
            torch.nn.Tanh()
            )

        self.conv2_merge = torch.nn.Sequential(
            torch.nn.Conv2d(16, 8, 3,padding=1),
            # torch.nn.BatchNorm2d(8),
            # torch.nn.Tanh()
            torch.nn.Tanh()

            # torch.nn.ReLU()

            )

        # self.deconv_6 = torch.nn.Sequential(
        #     torch.nn.ConvTranspose2d(256, 128, 3),
        #     # torch.nn.BatchNorm2d(16),
        #     # torch.nn.Tanh()
        #     torch.nn.Tanh()
        #     # torch.nn.ReLU()
        #     )

        # self.deconv_5 = torch.nn.Sequential(
        #     torch.nn.ConvTranspose2d(128, 64, 3),
        #     # torch.nn.BatchNorm2d(16),
        #     # torch.nn.Tanh()
        #     torch.nn.Tanh()
        #     # torch.nn.ReLU()
        #     )


        # self.deconv_4 = torch.nn.Sequential(
        #     torch.nn.ConvTranspose2d(64, 32, 3),
        #     # torch.nn.BatchNorm2d(16),
        #     # torch.nn.Tanh()
        #     torch.nn.Tanh()
        #     # torch.nn.ReLU()

        #     )


        self.deconv_3 = torch.nn.Sequential(
            torch.nn.ConvTranspose2d(32, 16, 3),
            # torch.nn.BatchNorm2d(16),
            # torch.nn.Tanh()
            torch.nn.Tanh()
            # torch.nn.ReLU()

            )

        self.deconv_2 = torch.nn.Sequential(
            torch.nn.ConvTranspose2d(16, 8, 3),
            # torch.nn.BatchNorm2d(8),
            # torch.nn.Tanh()
            torch.nn.Tanh()

            # torch.nn.ReLU()

            )

        self.deconv_1 = torch.nn.Sequential(
            torch.nn.ConvTranspose2d(8, 1, 3),
            # torch.nn.BatchNorm2d(1)
            # torch.nn.Tanh()
            )
            # torch.nn.ReLU()
            # torch.nn.ReLU()
            # )

    def forward(self, x):
        x = self.input(x)#10
        c1 = self.conv1(x)#9--------------------------|c1
        c2 = self.conv2(c1)#8-----------------|c2     |
        c3 = self.conv3(c2)#7                 |       |
        dco3 = self.deconv_3(c3)#8------------|dco3   |
        concat3 = torch.cat([dco3,c2],dim=1)#8|       |
        merge3 = self.conv3_merge(concat3)#8---------<|merge3

        dco2 = self.deconv_2(merge3)#9----------------|dco2
        concat2 = torch.cat([dco2,c1],dim=1)#9        |
        merge2 = self.conv2_merge(concat2)#9---------<|merge2

        dco1 = self.deconv_1(merge2)#10
        out = dco1.squeeze(1)
        return out

        
class MulityCNN(object):
    """docstring for MulityCNN"""
    def __init__(self, div_range,comm):
        super(MulityCNN, self).__init__()
        self.div_range = div_range
        self.comm = comm
        self.model = {}
        self.data_sz = None
        self.halo = 10

    def get_arg_from_rank(self,key,begin_div = 2):

        if key == 0:
            tar = 'phs'
            lev = 0
        elif key <= 32:
            tar = 'pt'
            lev = key - 1 
        elif key <= 64:
            tar = 'u_Agrid'
            lev = key - 1 - 32
        else:
            tar = 'v_Agrid'
            lev = key - 1 - 64

        return tar,lev

    def set_input_m(self,tar,model = '',lev = 32,div = 1,bias = "",model_save_path = '/share1/liminyan/model'):
        # comm = MPI.COMM_WORLD
        comm_rank = self.comm.Get_rank()
        comm_size = self.comm.Get_size()
        num = comm_rank
        num = 0

        last = '.npy'
        if bias == 'cnn':
            last = '.pt'
        if tar == 'phs':
            lev = 0
        if model == '0':
            path = model_save_path + '/'+ model +bias+tar+'_'+str(div) +'_'+str(lev)+'_[0]'+'_'+str(num)+last
        else:
            path = model_save_path + '/'+ model +bias+tar+'_'+str(div) +'_'+str(lev)+''+last


        return path

    def load(self,path):
        ord_ = 0
        s = time.time()
        num = self.comm.Get_rank()

        print(self.data_sz)
        for x in range(1,self.div_range+1):

            div = x

            div_eatch = int(32 / self.div_range )

            div = x*div_eatch

            tar,lev =  self.get_arg_from_rank(num,begin_div = 2)
            svr = ConvDeconvNeuralNetwork(CDNNArch(),epochs = 50)
            svr.tar = tar
            svr.lev = lev
            model_path = self.set_input_m(tar,
                model = str(ord_),lev = lev,div = div * 16,
                model_save_path = path,
                bias = 'cnn')

            ibegin = self.data_sz[-1-5]
            iend = self.data_sz[-1-4]
            jbegin = self.data_sz[-1-3]
            jend = self.data_sz[-1-2]

            svr.model_save_path = path
            svr.load(model_path,ibegin,iend,jbegin,jend,self.halo)
            key = str(div) +'_' + str(num)
            self.model[key] = svr

            
        e = time.time()
        if num == 0:
            print('load time',e-s)

    def convolution2d(self,image, kernel):
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

    def denoise(self,data):

        noisy = data
        sigma_est = np.mean(estimate_sigma(noisy))
        patch_kw = dict(patch_size=5,      # 5x5 patches
                        patch_distance=6,  # 13x13 search area
                        )
        denoise_fast = denoise_nl_means(noisy, h=0.8 * sigma_est, fast_mode=True,**patch_kw)
        return denoise_fast


    def predict(self,data):
        num = self.comm.Get_rank()
        size = self.comm.Get_size()

        ibegin = self.data_sz[-1-5]
        iend = self.data_sz[-1-4]
        jbegin = self.data_sz[-1-3]
        jend = self.data_sz[-1-2]

        result = np.empty([self.div_range,self.data_sz[-1-1],iend - ibegin+1, jend - jbegin+1])

        # kernel2 = np.ones([2,2]) / 4.0
        # kernel3 = np.ones([3,3]) / 9.0
        # kernel5 = np.ones([5,5]) / 25.0
        # kernel7 = np.ones([7,7]) / 49.0


        for x in range(1,self.div_range+1):
            div = x 

            div_eatch = int(32 / self.div_range )

            div = x*div_eatch

            key = str(div) +'_' + str(num)


            mid_res =  self.model[key].predict(data)

            if self.data_sz[-1-1] == 1:
                result[x - 1]= mid_res[0,self.halo:-self.halo,self.halo:-self.halo]
            else:
                result[x - 1]= mid_res[:,self.halo:-self.halo,self.halo:-self.halo]

            # mean_std = np.mean(self.model[key].std_X)
            # noisy = mid_res[0,halo:-halo,halo:-halo]
            # sigma = [mean_std*10, mean_std*3]
            # result[x - 1][0] = sp.ndimage.filters.gaussian_filter(mid_res[0], [1,1], mode='reflect')[halo:-halo,halo:-halo]

            # result[x - 1]  = self.denoise(noisy)
            # result[x - 1][0] = self.denoise(mid_res[0,halo : -halo,halo : -halo])
            # result[x -1][0] = self.convolution2d(mid_res[0,halo :-halo+1,halo :-halo+1],kernel2)
            # result[x - 1][0] = self.convolution2d(mid_res[0,halo-1:-halo+1,halo - 1:-halo+1],kernel3)
            # result[x - 1][0] = self.convolution2d(mid_res[0,halo-2:-halo+2,halo - 2:-halo+2],kernel5)
            # result[x - 1][0] = self.convolution2d(mid_res[0,halo-3:-halo+3,halo - 3:-halo+3],kernel7)
            # None

        return result



