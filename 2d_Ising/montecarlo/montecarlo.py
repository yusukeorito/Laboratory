# coding: UTF-8
"""
２次元イジングモデルのメトロポリス法によるモンテカルロシミュレーション
"""
from random import random, randrange
import random
import numpy as np
import matplotlib.pyplot as plt
from numpy.core.numerictypes import nbytes

J = 1 #強磁性相互作用
beta = 0.5
Nx = 32
Ny = 32
steps = 10**3

def shuffle_list(lst):
    lst2 = lst[:]
    random.shuffle(lst2)
    return lst2



#隣接サイトのスピンの合計を計算する関数
def neighbor_spin_sum(S,x,y):
    x_right = x + 1
    x_left = x - 1
    y_up = y + 1
    y_down = y - 1 

    #周期的境界条件
    if x_right >= Nx:
        x_right -= Nx
    if x_left < 0:
        x_left += Nx 
    if y_up >= Ny:
        y_up -= Ny 
    if y_down < 0:
        y_down += Ny 

    neighbor_spin_sum = S[x_right][y] + S[x_left][y] + S[x][y_up] + S[x][y_down]
    return neighbor_spin_sum


#スピン配位からエネルギーを計算する関数
def cal_E(S, Nx=Nx, Ny=Ny, J=J):
    E = 0
    for x in range(Nx):
        for y in range(Ny):
            E += - neighbor_spin_sum(S, x, y)/2
    E += np.sum(S)
    return E

#メトロポリス法
#ランダムに全てのサイトをフリップするアルゴリズム
def metropolis(S, beta=beta):
    xs = list(range(Nx))
    xs = shuffle_list(xs)
    ys = list(range(Ny))
    ys = shuffle_list(ys)
    
    for x in xs:
        for y in ys:
            sum = neighbor_spin_sum(S,x,y)
            transproba = np.exp(-2*beta*S[x][y]*sum)

            if np.random.random() <= transproba:
                S[x][y] = - S[x][y]
    return S

#モンテカルロ法
def mcmc(S, steps=1000, beta=beta,interval=10, burn_in=100):
    ms = []
    energies = []
    square_energies = []
    for step in range(steps):
        S = metropolis(S, beta=beta) #１montecalro_stepごとの配位が返ってくる
        if (step%interval == 0) & (step>=burn_in):

            m = np.sum(S) / (Nx*Ny)
            energy = cal_E(S, Nx=Nx,Ny=Ny)
            square_energy = energy * energy

            ms.append(m)
            energies.append(energy)
            square_energies.append(square_energy)

    return ms, energies, square_energies


S = np.ones((Nx,Ny)).tolist()
ms, energies, square_energies = mcmc(S, steps=steps, beta=beta)

time = np.arange(1, len(ms)+1)
fig, ax = plt.subplots(figsize=(6,4))

##熱力学量の温度変化をシミュレーションする関数
def beta_mcmc(S, betas, total_steps=10**6, burn_in=10**5, interval=20):
    #各逆温度に対して熱力学量の平均をプロットする
    beta_ms = []
    beta_energies = []
    beta_capacities = []

    for beta in betas:
        ms, energies, square_energies = mcmc(S, total_steps, beta, interval=interval, burn_in=burn_in)
        b_mag = np.mean(ms)
        b_energy = np.mean(energies)
        b_square_energy = np.mean(square_energies)
        b_capacity = (beta**2)*(b_square_energy - b_energy**2)
        
        beta_ms.append(b_mag)
        beta_energies.append(b_energy)
        beta_capacities.append(b_capacity)

    return beta_ms, beta_energies, beta_capacities


if __name__ == '__main__':

    
    
    ax.plot(time, ms)
    ax.set_ylim([-1.1, 1.1])
    ax.set_xlabel('MCS')
    ax.set_ylabel('Magnetization',fontsize=15)
    ax.set_title('Magnetization',fontsize=15)
    #ax.legend('beta = 0.1')
    #ax.grid()
    plt.tick_params(axis='y', width=1, length=4, pad=5, color='k', direction='in', labelsize=15, labelcolor='k')
    plt.tick_params(axis='x', width=1, length=4, pad=5, color='k', direction='in', labelsize=15, labelcolor='k')
    plt.show()
    
    

    Nx = 32
    Ny = 32
    
    S = np.ones((Nx, Ny)).tolist()
    betas = np.linspace(0.05, 1, 20).tolist()

    beta_ms, beta_energies, beta_capacities = beta_mcmc(S, betas, total_steps=10**3)

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12,4))
    betas = np.array(betas)
    beta_ms = np.array(beta_ms)
    beta_capacities = np.array(beta_capacities)
    ax1.plot(betas, beta_ms)
    ax1.set_xlabel(r'$\beta$',)
    ax1.set_ylabel('magnetization')
    ax1.grid()

    ax2.plot(betas, beta_capacities,c='red')
    ax2.set_xlabel(r'$\beta$',)
    ax2.set_ylabel('specific heat')
    ax2.grid()
    #plt.tick_params(axis='y', width=1, length=4, pad=5, color='k', direction='in', labelsize=15, labelcolor='k')
    #plt.tick_params(axis='x', width=1, length=4, pad=5, color='k', direction='in', labelsize=15, labelcolor='k')
    
    plt.show()
    


