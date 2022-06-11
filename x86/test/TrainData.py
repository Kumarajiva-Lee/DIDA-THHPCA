import numpy as np
import time

class TrainData(object):
    """docstring for TrainData"""
    def __init__(self, ensemble_nums=4, lats=180, lons=360):

        self.ensemble_nums = ensemble_nums
        self.lats = lats
        self.lons = lons
        self.train_months = None
        self.test_months = None
        self.x_train = None
        self.y_train = None
        self.x_test = None
        self.y_test = None
        self.data_label = {'x_train':0,'y_train':0,'x_test':0,'y_test':0,'train_months':0,'test_months':0}


    def show(self):
        """
        打印数据处理器信息，包括变量，输入输出纬度，数据加载状态
        """
        print('-------------------------------------------------')
        print('input shape: (%d, %d, %d)' % (self.ensemble_nums,self.lats,self.lons))
        print('output shape: (%d, %d)' % (self.lats,self.lons))
        print('data loading statu:')

        print('  loaded:')
        for key in self.data_label.keys():
            if self.data_label[key]:
                print('    ' + key + ': %s' % str(eval('self.'+key).shape))

        print('  unloaded:')
        for key in self.data_label.keys():
            if not self.data_label[key]:
                print('    ' + key + ': %s' % str(eval('self.'+key)))

        print('-------------------------------------------------')


    def loadTrainData(self,x_train=None,y_train=None):
        if type(x_train) is np.ndarray or type(x_train):
            self.x_train = x_train
            self.data_label['x_train'] = 1
        if type(y_train) is np.ndarray or type(y_train):
            self.y_train = y_train
            self.data_label['y_train'] = 1


    def loadData(self,from_path=False,x_train=None,y_train=None,x_test=None,y_test=None):
        s = time.time()
        if from_path:
            if x_train:
                self.x_train = np.load(x_train)
                self.data_label['x_train'] = 1
            if y_train:
                self.y_train = np.load(y_train)
                self.data_label['y_train'] = 1
            if x_test:
                self.x_test = np.load(x_test)
                self.data_label['x_test'] = 1
            if y_test:
                self.y_test = np.load(y_test)
                self.data_label['y_test'] = 1

            # return self

        else:
            if type(x_train) is np.ndarray :
                self.x_train = x_train
                self.data_label['x_train'] = 1
            if type(y_train) is np.ndarray :
                self.y_train = y_train
                self.data_label['y_train'] = 1
            if type(x_test) is np.ndarray or x_test.shape != None:
                self.x_test = x_test
                self.data_label['x_test'] = 1
            if type(y_test) is np.ndarray or y_test.shape != None:
                self.y_test = y_test
                self.data_label['y_test'] = 1
        e = time.time()
        print('loadData time',e - s)


        
