import numpy as np

class PostProcessor(object):
    """
    集合后处理处理器
    Attributes:
        data: DataProcessor, 涉及需使用的相关数据与数据操作
        model: 后处理算法类，为进行后处理所使用的算法
        results: numpy.ndarray 保存上一次预测的结果

    """
    def __init__(self, DataProcessor=None, Model=None):
        self.data = DataProcessor
        self.model = Model
        self.results = None

    def fit(self):
        self.model.fit(self.data)
        return self

    def partial_fit(self,bias):
        self.model.fit(self.data,partial_fit = bias)
        return self


    def predict(self):
        self.results = self.model.predict(self.data)
        return self.results


    def save(self, path):
        self.model.save(path)
        return self

    def load(self, path):
        self.model.load(path)
        return self

    def setData(self, newDataProcessor):
        self.data = newDataProcessor
        return self

    def setModel(self, newModel):
        self.model = newModel
        return self

    def getResults(self):
        return self.results


    def calRMSE(self):
        """
        计算均方根误差
        :return: 返回均方根误差
        """
        rmse = np.sqrt(np.mean((self.results-self.data.y_test)**2))
        return rmse







