//２次元Ising模型 wolffアルゴリズム実装 //
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
const int L = 40;
const int N = L*L;
const int MCS = N;
const int sample_num = 100000;
const int max_step = sample_num*10 + MCS;

const double J = 1.0;

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

void make_initial_up(int S[]){
    for (int i=0; i<N; i++){
        S[i] = 1;
    }
}


double Energy(int S[], int pair[][4]){
    double E = 0.0;
    for (int i =0; i<N; i++){
        for (int j=0; j<4; j++){
            E += -1*S[i]*S[pair[i][j]];
        }
    }
    E /= 2.0;
    /* E /= (double)N; */
    return E;
}


//単位磁化を計算する//
double Magnetization(int S[]){
    double M = 0.0;
    for(int i=0; i<N; i++){
        M += S[i];
    }
    M/= (double)N; //単位磁化
    return M;
}


int make_cluster(int S[N], double J, double T, int& n_cluster, int(&i_cluster)[N][1], mt19937 &engine){
    int in_out[N];
    int pair[N][4];
    make_pair(pair);


    //1ならばクラスターに属さず、0ならばクラスターに属する//
    for(int i=0; i<N; i++){
        in_out[i] =1; //初めは全てクラスターに属さない
    }

    uniform_int_distribution<int> dist1(0, N-1);
    uniform_real_distribution<double> dist2(0, 1);
    int chosen_site = dist1(engine); //ランダムに点を１つ選ぶ
    int spin_cluster = S[chosen_site];//選んだ点のスピンの値をspin_clusterに格納spin_cluster = +-1

    in_out[chosen_site] = 0; //選んだ点をクラスターに追加//
    i_cluster[0][0] = chosen_site;
    n_cluster = 1;

    double prob = 1.0 - exp(-2.0*J/T);
    int k=0;
    while(k<n_cluster){
        chosen_site = i_cluster[k][0];
        for(int j=0; j<4; j++){
            if(S[pair[chosen_site][j]] == spin_cluster){
                if(in_out[pair[chosen_site][j]]==1){//まだクラスターに追加されていないものが対象
                    if(dist2(engine) < prob){
                        i_cluster[n_cluster][0] = pair[chosen_site][j];
                        n_cluster = n_cluster+1;
                        in_out[pair[chosen_site][j]] = 0;
                    }
                }
            }
        }
        k = k+1;
    }
    return spin_cluster;
}

void MCMC(int S[N], int pair[][4], int n_cluster, int spin_cluster, int i_cluster[N][1], double beta, mt19937 &engine){
    uniform_real_distribution<double> dist2(0, 1);

    for (int k=0; k<n_cluster; k++){
        int i = i_cluster[k][0]; //クラスターに入っているスピンのインデックスを取り出す
        S[i] *= -1; //スピンをフリップ
    }  
}



int main(){
    int S[N];
    int pair[N][4];
    int seed = 2;


    double E;
    double M;
    double T=2.26;
    double beta = 1.0/T;

    random_device seed_gen;
    mt19937 engine(seed);

    ofstream file1;
    file1.open("../output/2d_Ising/wolff_energy_mag_"+ to_string(seed)+"_" + to_string(T).substr(0,4) + ".dat", ios::app);

    make_initial_up(S);
    make_pair(pair);

    E = Energy(S, pair);
    M = Magnetization(S);

    file1 << E << " " << abs(M) << endl;

    for(int step=0; step<1000*MCS; step++){
        int n_cluster; //クラスターに属するスピンの個数
        int i_cluster[N][1]; //クラスターに属するスピンのインデックス
        int spin_cluster = make_cluster(S, J, T, n_cluster, i_cluster, engine);
        //make_cluster１回ごとにi_clusterとn_clusterの値が書き換わる//

        MCMC(S, pair, n_cluster, spin_cluster, i_cluster, beta, engine);
        E = Energy(S, pair);
        M = Magnetization(S);

        if (step%1 == 0){
            file1 << E << " " << abs(M) << endl;
        }
    }
    file1.close();
    return 0;
}