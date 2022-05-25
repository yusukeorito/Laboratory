//制限付きの自己学習モンテカルロ
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>

using namespace std;
const int L = 240;
const int N = L*L;
const int MCS = N;
const int sample_num = 100000;
const int max_step = sample_num*10 + MCS;


const double J = 1.0;
const double h = 0.20;
const double J_eff = 1.1064;

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

//4体相互作用のエネルギーを計算//
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

//有効模型によるエネルギーの計算//
double Energy_eff(int S[], int pair[][4]){
    double E = 0.0;
    for (int i =0; i<N; i++){
        for (int j=0; j<4; j++){
            E += -J_eff*S[i]*S[pair[i][j]];
        }
    }
    E /= 2.0;
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

int get_euclid_distance(int i, int j){
    int x_i = i%L;
    int y_i = i/L;
    int x_j = j%L;
    int y_j = j/L;
    int dx = abs(x_i - x_j);
    int dy = abs(y_i - y_j);
    if (dx > L/2) dx = L - dx ;
    if (dy > L/2) dy = L - dy;
    int distance = dx + dy;
    return distance;
}

int make_cluster(int S[N], double J_eff, double T, double& weight,  int& n_cluster, int(&i_cluster)[N][1], mt19937 &engine){
    int in_out[N];
    int pair[N][4];
    make_pair(pair);
    
    //1ならばクラスターに属さず、0ならばクラスターに属する//
    for(int i=0; i<N; i++){
        in_out[i] =0; //初めは全てクラスターに属さない
    }

    uniform_int_distribution<int> dist1(0, N-1);
    uniform_real_distribution<double> dist2(0, 1);
    int chosen_site = dist1(engine); //ランダムに点を１つ選ぶ
    int spin_cluster = S[chosen_site];//選んだ点のスピンの値をspin_clusterに格納spin_cluster = +-1
    int distance;

    in_out[chosen_site] = 1; //選んだ点をクラスターに追加//
    i_cluster[0][0] = chosen_site; //ランダムに選んだクラスターの最初のサイト
    n_cluster = 1;

    double prob = 1.0 - exp(-2.0*J_eff/T);
    int k=0;
    while(k<n_cluster){
        chosen_site = i_cluster[k][0];
        for(int j=0; j<4; j++){
            if(S[pair[chosen_site][j]] == spin_cluster){//同符号のスピン
                if(in_out[pair[chosen_site][j]]==0){//まだクラスターに追加されていないものが対象
                    distance = get_euclid_distance(i_cluster[0][0],pair[chosen_site][j]);//初めのサイトとのユークリッド距離を計算
                    /*
                    if(distance<r){
                        if(dist2(engine) < prob){
                            i_cluster[n_cluster][0] = pair[chosen_site][j];
                            n_cluster = n_cluster+1;
                            in_out[pair[chosen_site][j]] = 1;
                        }
                    }else{
                        weight *= exp((-2/T)*J_eff*S[chosen_site]*S[pair[chosen_site][j]]);
                    }*/
                    if(dist2(engine) < prob){
                        i_cluster[n_cluster][0] = pair[chosen_site][j];
                        n_cluster = n_cluster+1;
                        in_out[pair[chosen_site][j]] = 1;
                    }
                }
            }
        }

        k = k+1;
    }
    return spin_cluster;
}




void restricted_SLMC(int S[N],  int pair[][4], int n_cluster, int spin_cluster, int i_cluster[N][1], double beta, double weight, mt19937 &engine){
    uniform_real_distribution<double> dist2(0, 1);
    int S_[N];
    for (int i=0; i<N; i++){
        S_[i] = S[i];
    }

    double E_org = Energy(S, pair); //flip前のoriginal模型のエネルギ-
    double E_eff = Energy_eff(S, pair);  //flip前の有効模型のエネルギ-
    double del_alpha = E_org - E_eff;
 

    for (int k=0; k<n_cluster; k++){
        int i = i_cluster[k][0]; //クラスターに入っているスピンのインデックスを取り出す
        S[i] *= -1; //スピンをフリップ   
    }  

    double E_org_ = Energy(S_, pair); //flip後のoriginal模型のエネルギ-
    double E_eff_ = Energy_eff(S_, pair); //flip後の有効模型のエネルギ-
    double del_beta = E_org_ - E_eff_;
    // double accept_ratio = exp(-1*beta*(del_beta-del_alpha));
    

    if (dist2(engine) < exp(-1*beta*(del_beta-del_alpha))*weight){
        for (int k=0; k<n_cluster; k++){
            int i = i_cluster[k][0]; //クラスターに入っているスピンのインデックスを取り出す
            S[i] *= -1; //スピンをフリップ 
        } 
    }   
}

int main(){
    int S[N];
    int pair[N][4];
   
    double E;
    double M;
    double T=2.49;
    double beta = 1.0/T;
    for (int seed=2022; seed<2027; seed++){
        random_device seed_gen;
        mt19937 engine(seed);
        ofstream file1;

        file1.open("./output/L240/seed"+to_string(seed)+"_restSLMC/energy_mag_"+ to_string(seed)+"_" + to_string(T).substr(0,4) + ".dat", ios::app);

        make_initial_random(S, engine);
        make_pair(pair);

        E = Energy(S, pair);
        M = Magnetization(S);

        file1 << E << " " << abs(M) << endl;

        for(int step=0; step<50000; step++){
            int n_cluster; //クラスターに属するスピンの個数
            int i_cluster[N][1]; //クラスターに属するスピンのインデックス

            double weight =1.0;//重みの初期化
            int spin_cluster  = make_cluster(S, J_eff, T, weight, n_cluster, i_cluster, engine);
            //make_cluster１回ごとにi_clusterとn_clusterの値が書き換わる//
            
            restricted_SLMC(S, pair, n_cluster, spin_cluster, i_cluster, beta, weight, engine);
            E = Energy(S, pair);
            M = Magnetization(S);

            if (step%1 == 0){
                cout << E << " " << abs(M) << endl;
            }
        }
        file1.close();
    }
    return 0;
}