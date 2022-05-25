//　4体相互作用を持つ模型のマルコフ連鎖モンテカルロ法 //
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
const int L = 80;
const int N = L*L;
const int MCS = N;
const int sample_num = 100000;
const int max_step = sample_num*10 + MCS;
const double J = 1.0;
const double h = 0.20;

//周期的境界条件//
void make_pair(int pair[][4]){
    for(int i = 0; i<N; i++){
        if((i+1)%L == 0){
            pair[i][0] = i- L+1;
        }else{
            pair[i][0] = i+1;
        }
        if(L*L -L <= i){
            pair[i][1] = i-(L*L -L);
        }else{
            pair[i][1] = i+L;
        }
        if(i%L == 0){
            pair[i][2] = i+(L-1);
        }else{
            pair[i][2] = i-1;
        }
        if(i < L){
            pair[i][3] = i+(L*L-L);
        }else{
            pair[i][3] = i-L;
        }
    }
}

// 初期状態all up //
void make_initial_up(int S[]){
    for (int i=0; i<N; i++){
        S[i] = 1;
    }
}

//　初期状態random //
void make_initial_random(int S[], mt19937 & engine){
    uniform_real_distribution<double> dist(0.1);
    for (int i=0; i<N; i++){
        if (dist(engine)<0.50){
            S[i] = -1;
        }else{
            S[i] = 1;
        }
    }
}

double Energy(int S[], int pair[][4]){
    double E1 = 0.0;
    double E2 = 0.0;
    for (int i =0; i<N; i++){
        for (int j=0; j<4; j++){
            E1 += -J*S[i]*S[pair[i][j]];
        }
    }
    
    for(int i=0; i<N; i++){
        for(int j=0; j<3; j+=2){
            for(int k=1; k<4; k+=2){
                E2 += -h*S[i]*S[pair[i][j]]*S[pair[pair[i][j]][k]]*S[pair[i][k]];
            }
        }
    }
    

   E1 /=2;
   E2 /=4;
   double E = E1 + E2;
    /* E /= (double)N; */
    return E;
}

double Magnetization(int S[]){
    double M = 0.0;
    for(int i=0; i<N; i++){
        M += S[i];
    }
    M/= (double)N; //単位磁化
    return M;
}

void MCMC(int S[], int pair[][4], double& E, double&M, double beta, mt19937 &engine){
    uniform_int_distribution<int> dist1(0, N-1);
    uniform_real_distribution<double> dist2(0, 1);

    int chosen_site = dist1(engine);
    int h1 = 0;
    double h2 = 0.0;

    for(int j=0; j<4; j++){
        h1 += J*S[pair[chosen_site][j]];//再近接相互作用の変化
    }

    for(int j=0; j<3; j+=2){
        for (int k=1; k<4; k+=2){
            h2 += h*S[pair[chosen_site][j]]*S[pair[pair[chosen_site][j]][k]]*S[pair[chosen_site][k]];
            //4体相互作用による変化//
        }
    }

    double delE = 2.0*(S[chosen_site]*h1 + S[chosen_site]*h2);
    double delM = -2.0*S[chosen_site];
    

    if (dist2(engine) < exp((-1*beta*delE))){
        S[chosen_site] *= -1;
        E += delE;
        M = M + delM/(double)N;
    }   
}

int main(){
    double E;
    double M;
    

    int seed =2028;
    int S[N];
    int pair[N][4];
    random_device seed_gen;
    mt19937 engine(seed);

    double T= 2.49;
    double beta = 1.0/T; 
    ofstream file1;
    file1.open("./output/L80/seed2028/Energy_Magnetization_L" + to_string(L) + "_T" + to_string(T).substr(0,4) +  "_" + to_string(seed) +".dat", ios::app);

    make_initial_up(S);
    make_pair(pair);

    E = Energy(S, pair);
    M = Magnetization(S);
    //file1 << E << " " <<  abs(M) << endl; 
    for (int step=0; step<100000*MCS; step++){
        MCMC(S, pair, E, M, beta, engine);
        if (step%MCS==0){
            file1 << E << " "  << abs(M) << endl;
        }
    }
    file1.close();  
}