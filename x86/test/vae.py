import torch as torch
import torch.nn as nn
import torch.utils.data

import numpy as np

mean = np.array([[95938.22234576],[1237.82497192],[957.16751403],[766.82583786],[641.08802845],[557.7870546],[500.78660246],[459.90315296],[430.17509624],[408.41704764],[391.42358218],[376.72257395],[363.25575287],[350.92643651],[340.15165774],[331.28988181],[324.29501513],[318.71699772],[313.99168821],[309.92327864],[306.43138158],[303.32326981],[300.47765614],[297.80472253],[295.28024783],[292.87532706],[290.53989666],[288.25810809],[286.04063418],[283.89449504],[281.77156853],[279.77418142],[277.01563663],[-13.00100971],[-3.63625814],[-0.2777924],[0.24179602],[-0.10180372],[0.31154],[2.04887822],[4.29837176],[6.28166199],[7.92324461],[9.3251559],[10.35242359],[10.74551872],[10.46715403],[9.75890739],[8.95465254],[8.19509963],[7.41520482],[6.63904785],[5.95421329],[5.38970177],[4.94399309],[4.61093324],[4.36176098],[4.14283342],[3.9336097],[3.75903115],[3.62946098],[3.52878517],[3.44534878],[3.37386557],[3.27971596],[-0.058163],[-0.04868889],[-0.01941438],[-0.0127029],[-0.00612373],[-0.00391849],[-0.00430699],[-0.00251479],[-0.00168673],[0.00082516],[-0.00190716],[-0.00346824],[-0.00721775],[-0.00699383],[-0.00745479],[-0.00367647],[-0.00039112],[-0.00461568],[-0.00791594],[-0.00869924],[-0.0088988],[-0.00948584],[-0.00899857],[-0.01027775],[-0.01128108],[-0.01374207],[-0.01841859],[-0.02324699],[-0.03211814],[-0.0418252],[-0.05279116],[-0.08090834]])

#mean = mean.reshape(97,1)

var = np.array([[99354104.80310565],[3397.75083504],[1743.48150421],[1022.71603525],[816.97346816],[779.6639258],[825.09996499],[862.5972393],[806.78515835],[652.00681061],[450.44505928],[256.22754482],[106.17649716],[26.14682646],[21.3759642],[67.23187],[119.56182543],[148.28246856],[151.41817451],[142.72535944],[132.0510347],[122.33262423],[114.05598232],[106.76461873],[100.3007061],[94.20679268],[87.87622511],[81.44554],[75.74285546],[71.37486742],[69.14780101],[69.54546221],[71.73932825],[516.99705041],[249.25721178],[152.75880768],[90.32309117],[71.2454967],[98.70235404],[121.51841153],[135.46149132],[152.98329962],[172.53215503],[189.94177073],[204.16625333],[218.13134354],[229.29508931],[229.01935439],[210.42516532],[179.17779284],[148.32002339],[124.83590774],[108.2761915],[95.48698158],[84.76517439],[75.29226165],[67.18936458],[61.10794967],[56.90455317],[53.99989178],[52.13929187],[51.20675668],[51.19147633],[51.83615518],[54.9917528],[128.33409431],[59.40633249],[35.30220863],[27.04200536],[22.40721403],[19.33071618],[17.45435352],[16.40935014],[15.96067447],[15.79778001],[16.06506595],[16.65673558],[17.20626117],[17.37561032],[17.02852118],[16.16846399],[14.92540678],[13.57498318],[12.409377],[11.51748427],[10.821789],[10.31406311],[9.95255533],[9.74080425],[9.64376702],[9.65178423],[9.7976026],[10.10133298],[10.5848731],[11.33865671],[12.35788995],[14.87182646]])

#var = var.reshape(97,1)

class VAE(nn.Module):
    def __init__(self):
        super(VAE, self).__init__()
        self.encoder = nn.Sequential(
            nn.Conv2d(97,128,kernel_size=3,stride=2,padding=1),
            nn.LeakyReLU(0.2,inplace=True),
            nn.Conv2d(128,64,kernel_size=3,stride=2,padding=1),
            nn.LeakyReLU(0.2,inplace=True),
            nn.Conv2d(64,32,kernel_size=3,stride=1,padding=1),
            nn.LeakyReLU(0.2,inplace=True),
        )
        self.encoder_rnn = nn.RNN(32*45*90, 1000, 1, batch_first=True)
        self.fc_mu = nn.Linear(1000, 500)
        self.fc_var = nn.Linear(1000, 500)
        
        #self.z2h = nn.Linear(500, 1000)
        self.norm = nn.BatchNorm1d(500)
        self.decoder_rnn = nn.RNN(500, 1000, 1, batch_first=True)
        self.decoder_fc = nn.Linear(1000, 32*45*90)
        self.decoder = nn.Sequential(
            nn.ConvTranspose2d(32, 16, 4, 2, 1),
            nn.ReLU(inplace=True),
            nn.ConvTranspose2d(16, 1, 4, 2, 1),
        )

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5*logvar)
        eps = torch.randn_like(std)
        return mu + eps*std

    def forward(self, x):
        batch_size = x.shape[0]
        x = self.encoder(x.view(-1,97,180,360))
        input = x.view(batch_size, 8, -1)
        _, x = self.encoder_rnn(input)
        
        x = x.view(batch_size, -1)

        mu = self.norm(self.fc_mu(x))
        logvar = self.fc_var(x)

        z = self.reparameterize(mu, logvar)
        #z = z.mean(dim=0)

        #hidden = self.z2h(z)
        #decoder_input = torch.zeros(mu.shape[0], 1, 1, requires_grad=True)
        output, _ = self.decoder_rnn(z.view(batch_size,1,500).repeat(1,8,1))
        output = self.decoder_fc(output.view(batch_size,8, -1))
        output = self.decoder(output.view(-1, 32, 45, 90))
        return output.view(batch_size,8,180,360), mu, logvar

class Model(object):
    def __init__(self, comm):
        super(Model, self).__init__()
        self.comm = comm
        self.model = VAE()

    def load(self,path):
        num = self.comm.Get_rank()
        self.model.load_state_dict(torch.load('/ddnstor/xuewei/zzh/VAE/4/models/VAE{}.pt'.format(num)))

    def predict(self,data):#data.shape=(1,97,180,360)
        num = self.comm.Get_rank()
        #scaler1 = pickle.load(open("/ddnstor/xuewei/zzh/test/scalers/scaler", 'rb'))
        #scaler2 = pickle.load(open("/ddnstor/xuewei/zzh/test/scalers/scaler{}".format(num), 'rb'))
        #scaler1 = joblib.load("/ddnstor/xuewei/zzh/conv/scalers/scaler")
        #scaler2 = joblib.load("/ddnstor/xuewei/zzh/conv/scalers/scaler{}".format(num))
        data = data.reshape(97,180*360)
        data = (data - mean)/np.sqrt(var)
        data = data.reshape(1,1,97,180,360)
        initial = data[0,0,num]
        #data = np.resize(scaler1.transform(data.T).T ,(1,1,97,180,360))
        truth = torch.from_numpy(np.repeat(data,8,axis=1)).float()

        result = np.empty([31,1,180, 360])
        self.model.eval()

        batch_size = truth.shape[0]
        x = self.model.encoder(truth.view(-1,97,180,360))
        input = x.view(batch_size, 8, -1)
        _, x = self.model.encoder_rnn(input)
        
        x = x.view(batch_size, -1)

        mu = self.model.norm(self.model.fc_mu(x))
        logvar = self.model.fc_var(x)

        z = self.model.reparameterize(mu, logvar)

        for i in range(31):
            z = 0.7*z+0.3*torch.randn(z.shape)
            output, _ = self.model.decoder_rnn(z.view(1,1,500).repeat(1,8,1))
            output = self.model.decoder_fc(output.view(8, -1))
            output = self.model.decoder(output.view(output.shape[0], 32, 45, 90))
            output = output.detach().numpy()[:,0]
            #output = np.resize(output,(8*180*360,1))
            output = output*np.sqrt(var[num]) + mean[num]
            output = np.resize(output, (8,180,360))
            diff = output[7] - initial
            #result[i] = 5*(output[7] - initial)+initial
            result[i] = output[8-1]
        result_mean = np.mean(result,axis=0,keepdims=True)
        result = 10 * (result -result_mean) + result_mean
        return result