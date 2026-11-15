from google.colab import drive
drive.mount('/content/drive')
%cd /content/drive/My\ Drive/Ding_Lab

import os
from torch import from_numpy
from torch.utils.data import Dataset, DataLoader
from natsort import natsorted
from torchvision.transforms import ToTensor
from torch.nn.utils.parametrizations import spectral_norm
import numpy as np
import pickle

#%matplotlib inline
import argparse
import glob
import random
import torch
import torch.nn as nn
import torch.nn.parallel
import torch.backends.cudnn as cudnn
import torch.optim as optim
import torchvision.transforms as transforms
import torchvision.utils as vutils
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML

# Set random seed for reproducibility
manualSeed = 999
# manualSeed = random.randint(1, 10000) # use if you want new results
print("Random Seed: ", manualSeed)
random.seed(manualSeed)
torch.manual_seed(manualSeed)

torch.autograd.set_detect_anomaly(True)

# #@title
# tot = len(glob.glob('ML_data/train/*.npz'))
# for i_file,file in enumerate(sorted(glob.glob('ML_data/train/*.npz'))):
#     print(i_file, '/', tot)
#     temp_dic = dict(np.load(file))

#     sorted_inds = np.empty(temp_dic['tile_ST'].shape)
#     temp_dic['ST_target'] = np.zeros(temp_dic['tile_ST'].shape)
#     for i in range(sorted_inds.shape[1]):
#         for j in range(sorted_inds.shape[2]):
#             sorted_inds[:,i,j] = np.argsort(-temp_dic['tile_ST'][:,i,j])
#             if temp_dic['tile_ST'][-1,i,j]==1: # tissue
#                 temp_dic['ST_target'][sorted_inds[:100,i,j].astype(int), i, j] = 1
#     file_base = file.split('.npz')[0].split('/')[-1]
#     with open('ML_data/train_pkl/' + file_base + '.pickle', 'wb') as handle:
#         pickle.dump(temp_dic, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Define your own class LoadFromFolder
class LoadFromFolder(Dataset):
    def __init__(self, main_dir):

        # Set the loading directory
        self.main_dir = main_dir

        # List all images in folder and count them
        all_imgs = os.listdir(main_dir)
        self.total_imgs = natsorted(all_imgs)

    def __len__(self):
        # Return the previously computed number of images
        return len(self.total_imgs)

    def __getitem__(self, idx):
        data_loc = os.path.join(self.main_dir, self.total_imgs[idx])

        with open(data_loc, 'rb') as handle:
            data = pickle.load(handle)

        downsampled_HE = np.load(os.path.join(self.main_dir, '../' + str(data['patient_id']) + '_' + str(data['sample_id']) + '_downsampled_tiff.npy'))
        mask = data['mask']

        coordinate = np.zeros((38,62))
        coordinate[int(data['coordinate'][0]), int(data['coordinate'][1])] = 1

        return data_loc, downsampled_HE, mask, from_numpy(data['HE_tile']).float(), from_numpy(coordinate).float(), from_numpy(data['tile_ST']).float(), from_numpy(data['ST_target']).float()

BATCH_SIZE = 2
trainset = LoadFromFolder(main_dir="./ML_data/train_pkl/")
trainloader = DataLoader(trainset, batch_size=BATCH_SIZE, shuffle=True, num_workers=2, pin_memory=True)
# valset = LoadFromFolder(main_dir="./data/Val/")
# valloader = DataLoader(valset, batch_size=BATCH_SIZE, shuffle=False, num_workers=2, pin_memory=True)

data = next(iter(trainloader))
print(data[0][0])
print('downsampled_HE -- ', data[1].shape)
print('mask -- ', data[2].shape)
print('HE_tile -- ', data[3].shape)
print('coordinate -- ', data[4].shape)
print('ST_tile -- ', data[5].shape)
print('ST_target -- ', data[6].shape)

# Size of z latent vector (i.e. size of generator input)
nz = 100

# Beta1 hyperparam for Adam optimizers
beta1 = 0.5

# Number of GPUs available. Use 0 for CPU mode.
ngpu = 1

# Decide which device we want to run on
device = torch.device("cuda:0" if (torch.cuda.is_available() and ngpu > 0) else "cpu")

# custom weights initialization called on netG and netD
def weights_init(m):
    classname = m.__class__.__name__
    if classname.find('Conv') != -1:
        nn.init.normal_(m.weight.data, 0.0, 0.02)
    elif classname.find('BatchNorm') != -1:
        nn.init.normal_(m.weight.data, 1.0, 0.02)
        nn.init.constant_(m.bias.data, 0)

class Dense_Block(nn.Module):
  def __init__(self, in_channels, num_layers, kernal_size, padding, dropout):
      super(Dense_Block, self).__init__()
      self.num_layers = num_layers
      self.leakyrelu = nn.LeakyReLU(inplace = True)
      self.drop_first = nn.Dropout2d(dropout)

      self.conv_first = spectral_norm(nn.Conv2d(in_channels = in_channels, out_channels = 32, kernel_size = kernal_size, padding = padding))
      self.bn_first = nn.BatchNorm2d(32)

      self.conv_layers = nn.ModuleList()
      self.bn_layers = nn.ModuleList()
      self.do_layers = nn.ModuleList()
      channels = 32+in_channels
      for i in range(num_layers):
        self.conv_layers.append(spectral_norm(nn.Conv2d(in_channels = channels, out_channels = 32, kernel_size = kernal_size, padding = padding)))
        self.bn_layers.append(nn.BatchNorm2d(32))
        self.do_layers.append(nn.Dropout2d(dropout))
        channels += 32

  def forward(self, x):
      intermediate = torch.cat((self.drop_first(self.leakyrelu(self.bn_first(self.conv_first(x)))), x), 1)

      for i in range(self.num_layers):
        intermediate = torch.cat([intermediate, self.do_layers[i](self.leakyrelu(self.bn_layers[i](self.conv_layers[i](intermediate))))], 1)

      return intermediate

def conv_output_size(input_size, kernal, padding, stride):
    return np.divide((input_size-kernal+(2*padding)),stride)+(1,1)

input_size = np.array((6,3))
kernal = np.array((6,3))
padding = np.array((0,0))
stride = np.array((1,1))
conv_output_size(input_size, kernal, padding, stride)

class HE_Encoder(nn.Module):
    def __init__(self, ngpu):
        super(HE_Encoder, self).__init__()
        self.ngpu = ngpu

        # Input input is BATCH X (3) x 752 X 432
        self.dense_block_t4 = Dense_Block(in_channels=3, num_layers=3, kernal_size=(3,3), padding=(1,1), dropout=0.2)
        # state size. BATCH x (131) x 752 x 432

        self.he_net_t4 = nn.Sequential(
            spectral_norm(nn.Conv2d(131,10, (4,4), padding=(1,1), stride=2)),
            nn.BatchNorm2d(10),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (10) x 376 X 216
        )

        self.dense_block_t3 = Dense_Block(in_channels=10, num_layers=3, kernal_size=(3,3), padding=(1,1), dropout=0.2)
        # state size. BATCH x (138) x 376 x 216

        self.dis_net_t3 = nn.Sequential(
            spectral_norm(nn.Conv2d(138,20,(4,4), padding=(1,1), stride=2)),
            nn.BatchNorm2d(20),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. (20) x 188 x 108
        )

        self.dense_block_t2 = Dense_Block(in_channels=20, num_layers=3, kernal_size=(3,3), padding=(1,1), dropout=0.2)
        # state size. BATCH x (148) x 188 x 108

        self.dis_net_t2 = nn.Sequential(
            spectral_norm(nn.Conv2d(148,30,(4,4), padding=(1,1), stride=2)),
            nn.BatchNorm2d(30),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (30) x 94 x 54
        )

        self.dense_block_t1 = Dense_Block(in_channels=30, num_layers=3, kernal_size=(3,3), padding=(1,1), dropout=0.2)
        # state size. BATCH x (158) x 94 x 54

        self.final = nn.Sequential(
            spectral_norm(nn.Conv2d(158,64,(4,4), padding=(2,2), stride=2)),
            nn.BatchNorm2d(64),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (64) x 48 x 28

            spectral_norm(nn.Conv2d(64,64,(4,6), padding=(1,0), stride=2)),
            nn.BatchNorm2d(64),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (64) x 24 x 12

            spectral_norm(nn.Conv2d(64,64,(4,4), padding=(1,1), stride=2)),
            nn.BatchNorm2d(64),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (64) x 12 x 6

            spectral_norm(nn.Conv2d(64,64,(4,4), padding=(1,1), stride=2)),
            nn.BatchNorm2d(64),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (64) x 6 x 3

            spectral_norm(nn.Conv2d(64,64,(6,3), padding=(0,0), stride=1)),
            nn.BatchNorm2d(64),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (64) x 1 x 1
        )

    def forward(self, HE_tile):
        intermediate = self.he_net_t4(self.dense_block_t4(torch.permute(HE_tile, (0,3,1,2))))
        intermediate = self.dis_net_t3(self.dense_block_t3(intermediate))
        intermediate = self.dis_net_t2(self.dense_block_t2(intermediate))
        intermediate = self.final(self.dense_block_t1(intermediate))
        return intermediate

HE_en = HE_Encoder(ngpu).to(device)

# Handle multi-gpu if desired
if (device.type == 'cuda') and (ngpu > 1):
    HE_en = nn.DataParallel(HE_en, list(range(ngpu)))

# Apply the weights_init function to randomly initialize all weights
#  to mean=0, stdev=0.02.
HE_en.apply(weights_init)

# Print the model
print(HE_en)

data = next(iter(trainloader))
print(data[0][0])
downsampled_HE = data[1]
mask = data[2]
HE_tile = data[3].to(device)
coordinate = data[4]
ST_tile = data[5]

HE_en(HE_tile).shape

class ST_Encoder(nn.Module):
    def __init__(self, ngpu):
        super(ST_Encoder, self).__init__()
        self.ngpu = ngpu

        self.st_net_conv = nn.Sequential(
            # input size BATCH x (17942) x 3 x 3
            spectral_norm(nn.Conv2d(17942,2000, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(2000),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (2000) x 3 x 3

            spectral_norm(nn.Conv2d(2000,1000, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(1000),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (1000) x 3 x 3

            spectral_norm(nn.Conv2d(1000,500, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(500),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (500) x 3 x 3

            spectral_norm(nn.Conv2d(500,250, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(250),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (250) x 3 x 3

            spectral_norm(nn.Conv2d(250,125, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(125),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (125) x 3 x 3

            spectral_norm(nn.Conv2d(125,64, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(64),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (64) x 3 x 3

            spectral_norm(nn.Conv2d(64,64, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(64),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (64) x 3 x 3

            spectral_norm(nn.Conv2d(64,64, (3,3), padding=(0,0), stride=1)),
            nn.BatchNorm2d(64),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (64) x 1 x 1
        )

    def forward(self, ST_tile):
        return self.st_net_conv(ST_tile)

ST_en = ST_Encoder(ngpu).to(device)

# Handle multi-gpu if desired
if (device.type == 'cuda') and (ngpu > 1):
    ST_en = nn.DataParallel(ST_en, list(range(ngpu)))

# Apply the weights_init function to randomly initialize all weights
#  to mean=0, stdev=0.02.
ST_en.apply(weights_init)

# Print the model
print(ST_en)

data = next(iter(trainloader))
print(data[0][0])
downsampled_HE = data[1]
mask = data[2]
HE_tile = data[3].to(device)
coordinate = data[4]
ST_tile = data[5]

ST_en(ST_tile).shape

input_size = np.array((1022,1150))
kernal = np.array((3,20))
padding = np.array((0,0))
stride = np.array((1,1))
conv_output_size(input_size, kernal, padding, stride)

class HE_Decoder(nn.Module):
    def __init__(self, ngpu):
        super(HE_Decoder, self).__init__()
        self.ngpu = ngpu

        # self.coordinate_net = nn.Sequential(
        #     # input is BATCH x (1) x 38 x 62
        #     nn.Conv2d(1,4, (4,4), padding=(2,2), stride=2, bias=False),
        #     nn.BatchNorm2d(4),
        #     nn.LeakyReLU(0.2, inplace=True),
        #     nn.Dropout2d(p=0.0),
        #     # state size Batch x 4 x 20 x 32

        #     nn.Conv2d(4,8,(4,4), padding=(1,1), stride=2, bias=False),
        #     nn.BatchNorm2d(8),
        #     nn.LeakyReLU(0.2, inplace=True),
        #     nn.Dropout2d(p=0.0),
        #     # state size. Batch x (8) x 10 x 16

        #     nn.Conv2d(8,16,(2,4), padding=(1,1), stride=2, bias=False),
        #     nn.BatchNorm2d(16),
        #     nn.LeakyReLU(0.2, inplace=True),
        #     nn.Dropout2d(p=0.0),
        #     # state size. BATCH x (16) x 6 x 8

        #     nn.Conv2d(16,32,(3,8), padding=(1,1), stride=1, bias=False),
        #     nn.BatchNorm2d(32),
        #     nn.LeakyReLU(0.2, inplace=True),
        #     nn.Dropout2d(p=0.0),
        #     # state size. BATCH x (32) x 6 x 3

        #     nn.Conv2d(32,1,(4,3), padding=(0,0), stride=2, bias=False),
        #     nn.BatchNorm2d(1),
        #     nn.LeakyReLU(0.2, inplace=True),
        #     nn.Dropout2d(p=0.0),
        #     # state size. BATCH x (1) x 2 x 1
        # )

        self.he_decoder = nn.Sequential(
            # input state BATCH x (2) x 8 x 8
            nn.ConvTranspose2d(2, 3, kernel_size=(6,3), stride=(2,1), padding=0),
            nn.BatchNorm2d(3),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (3) x 20 x 10

            nn.ConvTranspose2d(3, 3, kernel_size=(5,5), stride=2, padding=0),
            nn.BatchNorm2d(3),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (3) x 43 x 23

            nn.ConvTranspose2d(3, 3, kernel_size=(6,6), stride=2, padding=0),
            nn.BatchNorm2d(3),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (3) x 90 x 50

            nn.ConvTranspose2d(3, 3, kernel_size=(6,6), stride=2, padding=0),
            nn.BatchNorm2d(3),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (3) x 184 x 104

            nn.ConvTranspose2d(3,3, kernel_size=(6,6), stride=2, padding=0),
            nn.BatchNorm2d(3),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (3) x 372 x 212

            nn.ConvTranspose2d(3,3, kernel_size=(6,6), stride=2, padding=0),
            nn.BatchNorm2d(3),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. BATCH x (3) x 748 x 428

            nn.Conv2d(3,3,(3,3), padding=(2,2), stride=1, bias=False),
            nn.BatchNorm2d(3),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. BATCH x (3) x 750 x 430

            nn.Conv2d(3,3,(3,3), padding=(2,2), stride=1),
            nn.BatchNorm2d(3),
            nn.Sigmoid(),
            # state size. BATCH x (3) x 752 x 432
        )


    def forward(self, HE_latent, ST_latent):
        # coordinate = self.coordinate_net(torch.unsqueeze(coordinate, 1))
        intermediate = torch.cat((HE_latent.reshape(-1,1,8,8), ST_latent.reshape(-1,1,8,8)), dim=1)
        self.he_decoder(intermediate)
        return self.he_decoder(intermediate)

HE_De = HE_Decoder(ngpu).to(device)

# Handle multi-gpu if desired
if (device.type == 'cuda') and (ngpu > 1):
    HE_De = nn.DataParallel(HE_De, list(range(ngpu)))

# Apply the weights_init function to randomly initialize all weights
#  to mean=0, stdev=0.2.
HE_De.apply(weights_init)

# Print the model
print(HE_De)

data = next(iter(trainloader))
print(data[0][0])
downsampled_HE = data[1]
mask = data[2]
HE_tile = data[3].to(device)
coordinate = data[4].to(device)
ST_tile = data[5]

HE_De(HE_en(HE_tile), ST_en(ST_tile)).shape

ST_tile.shape

class ST_Decoder(nn.Module):
    def __init__(self, ngpu):
        super(ST_Decoder, self).__init__()
        self.ngpu = ngpu

        self.st_decoder = nn.Sequential(
            # input size BATCH x (2) x 8 x 8
            spectral_norm(nn.Conv2d(2,4, (3,3), padding=(0,0), stride=1)),
            nn.BatchNorm2d(4),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (4) x 6 x 6

            spectral_norm(nn.Conv2d(4,8, (3,3), padding=(0,0), stride=1)),
            nn.BatchNorm2d(8),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (8) x 4 x 4

            spectral_norm(nn.Conv2d(8,16, (2,2), padding=(0,0), stride=1)),
            nn.BatchNorm2d(16),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (16) x 3 x 3

            spectral_norm(nn.Conv2d(16,32, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(32),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (32) x 3 x 3

            spectral_norm(nn.Conv2d(32,64, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(64),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (64) x 3 x 3

            spectral_norm(nn.Conv2d(64,128, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(128),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (128) x 3 x 3

            spectral_norm(nn.Conv2d(128,256, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(256),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (256) x 3 x 3

            spectral_norm(nn.Conv2d(256,512, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(512),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (512) x 3 x 3

            spectral_norm(nn.Conv2d(512,1024, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(1024),
            nn.LeakyReLU(inplace=True),
            nn.Dropout2d(p=0.2),
            # state size. BATCH x (1024) x 3 x 3

            spectral_norm(nn.Conv2d(1024,2048, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(2048),
            nn.LeakyReLU(inplace=True),
            # state size. BATCH x (2048) x 3 x 3

            spectral_norm(nn.Conv2d(2048,17942, (3,3), padding=(1,1), stride=1)),
            nn.BatchNorm2d(17942),
            nn.Sigmoid(),
            # state size. BATCH x (17942) x 3 x 3
        )

    def forward(self, HE_latent, ST_latent):
        # coordinate = self.coordinate_net(torch.unsqueeze(coordinate, 1))
        intermediate = torch.cat((HE_latent.reshape(-1,1,8,8), ST_latent.reshape(-1,1,8,8)), dim=1)
        self.st_decoder(intermediate)
        return self.st_decoder(intermediate)

ST_De = ST_Decoder(ngpu).to(device)

# Handle multi-gpu if desired
if (device.type == 'cuda') and (ngpu > 1):
    ST_De = nn.DataParallel(ST_De, list(range(ngpu)))

# Apply the weights_init function to randomly initialize all weights
#  to mean=0, stdev=0.2.
ST_De.apply(weights_init)

# Print the model
print(ST_De)

data = next(iter(trainloader))
print(data[0][0])
downsampled_HE = data[1]
mask = data[2]
HE_tile = data[3].to(device)
coordinate = data[4].to(device)
ST_tile = data[5]

ST_De(HE_en(HE_tile), ST_en(ST_tile)).shape

torch.random()

def ST_tile_auto_loss(output, target):
    

# Initialize BCELoss function
HE_tile_auto_loss = nn.MSELoss()
cos_sim_loss = nn.CosineSimilarity(dim=1, eps=1e-6)

##good combos####
##############
optimizerE = optim.Adam(HE_en.parameters(), lr=0.0005, betas=(beta1, 0.999))
optimizerD = optim.Adam(net_Decoder.parameters(), lr=0.0005, betas=(beta1, 0.999))
##############


### SGD
# optimizerD = torch.optim.SGD(netD.parameters(), lr=0.001, momentum=0.9)
# optimizerG = torch.optim.SGD(netG.parameters(), lr=0.01, momentum=0.9)

# Setup Adam optimizers for both G and D
# optimizerD = optim.Adam(netD.parameters(), lr=0.005, betas=(beta1, 0.999))
# optimizerG = optim.Adam(netG.parameters(), lr=0.005, betas=(beta1, 0.999))

# Training Loop

num_epochs=1000
# Lists to keep track of progress
losses = []

print("Starting Training Loop...")
# For each epoch
for epoch in range(num_epochs):
    # For each batch in the dataloader
    for i, data in enumerate(trainloader, 0):

        HE_en.zero_grad()
        net_Decoder.zero_grad()

        downsampled_HE = data[1].float()
        mask = data[2].float()
        downsampled_HE[torch.where(mask!=1)]=0
        downsampled_HE = downsampled_HE.to(device) # size BATCH x 1013 x 1063 x 3
        HE_tile = data[3].float().to(device)
        coordinate = data[4].float().to(device)
        ST_tile = data[5].float().to(device)

        HE_latent = HE_en(HE_tile) # size BATCH x 1 x 2 x 1
        output = torch.permute(output, (0,2,3,1)) # size BATCH x 1013 x 1063 x 3
        loss = criterion(output, downsampled_HE)
        # Calculate gradients for D in backward pass
        loss.backward()
        optimizerE.step()
        optimizerD.step()
        losses.append(loss.item())

        # Output training stats
        if i % 2 == 0:
            print('[%d/%d][%d/%d]  Loss:%.2f' % (epoch, num_epochs, i, len(trainloader), loss.item()))

netG.eval()
for v, data in enumerate(valloader, 0):
  with torch.no_grad():

      # Real Batch
      val_vocals = torch.unsqueeze(data[1],1).to(device)
      val_chords = torch.unsqueeze(data[2],1).to(device)
      val_section = data[3].to(device)

      b_size = val_vocals.size(0)
      print('b_size:', b_size)
      labels_real = torch.full((b_size,), real_label, dtype=torch.float, device=device)

      # Fake Batch
      noise = torch.randn(b_size, nz, device=device)
      fake = netG(noise, val_chords, val_section)

q=6
print(data[0][q])
plt.figure(figsize=(8, 6), dpi=100)
plt.yticks(np.arange(len(cmajor_small_ref)), cmajor_small_ref)
plt.imshow(np.load(data[0][q])['chords_array'].T, interpolation='nearest', aspect='auto', origin='lower')

print(data[0][q])
plt.figure(figsize=(8, 6), dpi=100)
plt.yticks(np.arange(len(cmajor_small_ref)), cmajor_small_ref)
im = np.squeeze(fake[q].detach().cpu()).T
# im[im<0.5]=0
plt.imshow(im, interpolation='nearest', aspect='auto', origin='lower')

print(data[0][q])
plt.figure(figsize=(8, 6), dpi=100)
plt.yticks(np.arange(len(cmajor_small_ref)), cmajor_small_ref)
plt.imshow(np.load(data[0][q])['vocals_array'].T, interpolation='nearest', aspect='auto', origin='lower')

!mkdir /content/temp

for q in range(30):
  filename = data[0][q].split('./data/Train/')[1]
  print(filename)
  np.savez_compressed('/content/temp/temp.npz', vocals_array=torch.squeeze(fake.cpu().detach())[q])
  upload_file(filename, '/content/temp/temp.npz', '1FOvSNxUxNv04rr96bNLgnL-y-EkOwYRl')
  

torch.save({
            'epoch': epoch,
            'netG_state_dict': netG.state_dict(),
            'netD_state_dict': netD.state_dict(),
            'optimizerG_state_dict': optimizerG.state_dict(),
            'optimizerD_state_dict': optimizerD.state_dict(),
            'G_losses': G_losses,
            'D_losses': D_losses,
            'valD_losses': valD_losses,
            }, 'Models/dense_vanilla_midi_conditional_DCGAN_network_pytorch_with_section/model_0005_005_batch32_' + str(epoch).zfill(6) + '.pt')

plt.figure(figsize=(10,5))
plt.title("Discriminator Val Loss During Training")
plt.plot(valD_losses,label="D-val")
plt.xlabel("iterations")
plt.ylabel("Loss")
plt.legend()
plt.show()

plt.figure(figsize=(10,5))
plt.title("Generator and Discriminator Loss During Training")
plt.plot(G_losses,label="G")
plt.plot(D_losses,label="D")
plt.xlabel("iterations")
plt.ylabel("Loss")
# plt.ylim(0,100)
plt.legend()
plt.show()

#%%capture
fig = plt.figure(figsize=(6,6))
plt.axis("off")
ims = [[plt.imshow(i.T, interpolation='nearest', aspect='auto', origin='lower', animated=True)] for i in img_list]
ani = animation.ArtistAnimation(fig, ims, interval=1000, repeat_delay=1000, blit=True)
plt.yticks(np.arange(len(cmajor_small_ref)), cmajor_small_ref)

a = HTML(ani.to_jshtml())

html = a.data
with open('Maroon 5_Payphone_chorus_0.html', 'w') as f:
    f.write(html)

def plot_grad_flow(named_parameters):
    '''Plots the gradients flowing through different layers in the net during training.
    Can be used for checking for possible gradient vanishing / exploding problems.
    
    Usage: Plug this function in Trainer class after loss.backwards() as 
    "plot_grad_flow(self.model.named_parameters())" to visualize the gradient flow'''
    ave_grads = []
    max_grads= []
    layers = []
    for n, p in named_parameters:
        if(p.requires_grad) and ("bias" not in n):
            layers.append(n)
            try:
              ave_grads.append(p.grad.abs().mean())
              max_grads.append(p.grad.abs().max())
            except:
              pass
    plt.bar(np.arange(len(max_grads)), max_grads, alpha=0.1, lw=1, color="c")
    plt.bar(np.arange(len(max_grads)), ave_grads, alpha=0.1, lw=1, color="b")
    # plt.hlines(0, 0, len(ave_grads)+1, lw=2, color="k" )
    plt.xticks(range(0,len(ave_grads), 1), layers, rotation="vertical")
    # plt.xlim(left=0, right=len(ave_grads))
    # plt.ylim(bottom = -0.001, top=0.02) # zoom in on the lower gradient regions
    plt.xlabel("Layers")
    plt.ylabel("average gradient")
    plt.title("Gradient flow")

plot_grad_flow(netD.named_parameters())
