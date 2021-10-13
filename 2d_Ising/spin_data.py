import os
import sys
import random
import numpy as np
sys.path.append(".")
from montecarlo import *
from montecarlo import metropolis

class Config:
    DEBUG = True
    SEED = 2021

#シードの固定
def fix_seed(seed):
    random.seed(seed)
    np.random.seed(seed)

fix_seed(Config.SEED)



N_dis = 10**2
N_trj = 10**3 + N_dis
N_save = 10 #10回ごとに配位を保存する
prm_list = [
    [0.90, N_trj, N_dis, "conf/L32b090_", N_save],
    [0.85, N_trj, N_dis, "conf/L32b085_", N_save],
    [0.80, N_trj, N_dis, "conf/L32b080_", N_save],
    [0.75, N_trj, N_dis, "conf/L32b075_", N_save],
    [0.70, N_trj, N_dis, "conf/L32b070_", N_save],
    [0.65, N_trj, N_dis, "conf/L32b065_", N_save],
    [0.60, N_trj, N_dis, "conf/L32b060_", N_save],
    [0.55, N_trj, N_dis, "conf/L32b055_", N_save],
    [0.50, N_trj, N_dis, "conf/L32b050_", N_save],
    [0.45, N_trj, N_dis, "conf/L32b045_", N_save],
    [0.40, N_trj, N_dis, "conf/L32b040_", N_save],
    [0.35, N_trj, N_dis, "conf/L32b035_", N_save],
    [0.30, N_trj, N_dis, "conf/L32b030_", N_save],
    [0.25, N_trj, N_dis, "conf/L32b025_", N_save],
    [0.20, N_trj, N_dis, "conf/L32b020_", N_save],
    [0.15, N_trj, N_dis, "conf/L32b015_", N_save],
    [0.10, N_trj, N_dis, "conf/L32b010_", N_save],
    [0.05, N_trj, N_dis, "conf/L32b005_", N_save],
    [0.00, N_trj, N_dis, "conf/L32b000_", N_save],
]


if Config.DEBUG:
    N_dis = 10
    N_trj = 10**2 + N_dis
    N_save = 10
    N_x = 4
    N_y = 4

os.makedirs('conf', exist_ok=True)
nprm = len(prm_list)
betas = []
mags = []
mags_er = []
sc = np.ones((Nx, Ny))



for ibeta in range(nprm):
    beta = prm_list[ibeta][0]
    Nsweep = prm_list[ibeta][1]
    Ndiscard = prm_list[ibeta][2]
    file_name = prm_list[ibeta][3]
    save_every = prm_list[ibeta][4]
    conf_cnt = 0

    print(f'beta={beta} {Nsweep}')

    for isweep in range(Nsweep):
        sc = metropolis(sc, beta=beta) #各逆温度における配位が返ってくる

        if (isweep%save_every == 0) & (isweep>Ndiscard):
            scn = np.array(sc)
            np.savetxt(f"{file_name}{conf_cnt}", scn, fmt="%d")
            conf_cnt += 1
    
    print(f"generate spin Data {ibeta}")


"""


nconf = 100 #各温度に対する配位の生成数
betacr = 0.440686
data = []
labels = []
betas = []
nprm = len(prm_list)


for ibeta in  range(nprm):
    beta = prm_list[ibeta][0]
    file_name = prm_list[ibeta][3]

    for itrj in range(nconf):
        npsc = np.load(f"{file_name}{itrj}.npy") #各温度に対して100回分配位データを読み込む
        data.append(npsc)

        if beta > betacr:
            label = [0,1]
            labels.append(label)
        else:
            label = [1,0]
            labels.append(label)
        betas.append(beta)

data = np.array(data)
labels = np.array(labels)
        

"""

