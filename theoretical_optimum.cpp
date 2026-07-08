#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

// ビットの立っている数をカウント（ハミング重み）
int popcount(int v) {
    int count = 0;
    while (v) { count += v & 1; v >>= 1; }
    return count;
}

// グレイ符号の生成
int gray_code(int m, int M) {
    if (M == 2) return m;
    int half = M / 2;
    if (m < half) return gray_code(m, half);
    else return half + gray_code(M - 1 - m, half);
}

// 組み合わせ nCr の計算
long long nCr(int n, int r) {
    if (r > n) return 0;
    if (r == 0 || r == n) return 1;
    if (r > n / 2) r = n - r;
    long long res = 1;
    for (int i = 1; i <= r; ++i) {
        res = res * (n - i + 1) / i;
    }
    return res;
}

// 事前計算した差分ベクトルと重み
struct DiffVector {
    vector<double> v;
    double kappa;
};

// 差分ベクトルの再帰的な全列挙（任意のKに対応）
void build_vectors(int dim, int K, int M, const vector<int>& N_vals, const vector<int>& D_vals, vector<int>& current_v, vector<DiffVector>& valid_vectors) {
    if (dim == K) {
        bool all_zero = true;
        for (int x : current_v) if (x != 0) all_zero = false;
        if (all_zero) return; // 原点(0,0,..,0)は除外
        
        double H = 0;
        for (int k = 0; k < K; ++k) {
            double term = D_vals[current_v[k] + M - 1];
            for (int j = 0; j < K; ++j) {
                if (j != k) term *= N_vals[current_v[j] + M - 1];
            }
            H += term;
        }
        double kappa = H / (K * log2(M));
        if (kappa > 0) {
            DiffVector dv;
            dv.v.resize(K);
            for(int k=0; k<K; ++k) dv.v[k] = current_v[k];
            dv.kappa = kappa;
            valid_vectors.push_back(dv);
        }
        return;
    }
    // 各次元の差分 -(M-1) から M-1 まで
    for (int val = -(M - 1); val <= M - 1; ++val) {
        current_v[dim] = val;
        build_vectors(dim + 1, K, M, N_vals, D_vals, current_v, valid_vectors);
    }
}

// 任意のK次元に対するギブンス回転行列の積を生成
vector<vector<double>> get_rotation_matrix(int K, double theta) {
    // 単位行列で初期化
    vector<vector<double>> R(K, vector<double>(K, 0.0));
    for(int i=0; i<K; ++i) R[i][i] = 1.0;
    
    double c = cos(theta);
    double s = sin(theta);
    
    // すべての軸ペア(i, j)に対するGivens回転を適用
    for(int i=0; i<K-1; ++i) {
        for(int j=i+1; j<K; ++j) {
            vector<vector<double>> G(K, vector<double>(K, 0.0));
            for(int k=0; k<K; ++k) G[k][k] = 1.0;
            G[i][i] = c;
            G[i][j] = -s;
            G[j][i] = s;
            G[j][j] = c;
            
            // R = G * R
            vector<vector<double>> next_R(K, vector<double>(K, 0.0));
            for(int row=0; row<K; ++row) {
                for(int col=0; col<K; ++col) {
                    for(int k=0; k<K; ++k) {
                        next_R[row][col] += G[row][k] * R[k][col];
                    }
                }
            }
            R = next_R;
        }
    }
    return R;
}

// 一般化されたシミュレーション実行関数
void run_simulation(int K, int M, string name) {
    cout << "Starting simulation for K=" << K << ", " << name << " (M=" << M << ")..." << endl;
    
    // 式(70)のBER変換用のスケーリング係数を計算 (Eb/N0 = 30dB)
    double snr_db = 30.0;
    double snr_linear = pow(10.0, snr_db / 10.0);
    double coeff_base = (M * M - 1.0) / (12.0 * M * log2(M) * snr_linear);
    double scaling_factor = pow(coeff_base, K) * nCr(2 * K, K);
    
    // 1次元ごとの出現数(N)とハミング距離和(D)を事前計算
    vector<int> N_vals(2 * M - 1, 0);
    vector<int> D_vals(2 * M - 1, 0);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            int diff = j - i;
            N_vals[diff + M - 1]++;
            D_vals[diff + M - 1] += popcount(gray_code(i, M) ^ gray_code(j, M));
        }
    }
    
    // 有効な差分ベクトルと重みkappaを全列挙
    vector<DiffVector> valid_vectors;
    vector<int> current_v(K, 0);
    build_vectors(0, K, M, N_vals, D_vals, current_v, valid_vectors);
    
    string filepath = "C:\\Users\\Ide Nanako\\Desktop\\result3\\K" + to_string(K) + "_" + name + ".csv";
    ofstream out(filepath);
    out << "Angle(deg),BER_Bound_Eq70\n";
    
    double min_BER = 1e30;
    double opt_angle = 0;
    
    for (double deg = 0; deg <= 45.0001; deg += 0.01) {
        double th = deg * M_PI / 180.0;
        vector<vector<double>> R = get_rotation_matrix(K, th);
        
        double fR = 0;
        for (const auto& dv : valid_vectors) {
            double P = 1.0;
            // r = R * v を計算し、成分の2乗の積をとる
            for (int k = 0; k < K; ++k) {
                double rk = 0;
                for (int j = 0; j < K; ++j) {
                    rk += R[k][j] * dv.v[j];
                }
                P *= (rk * rk);
            }
            
            // ゼロ除算の回避（重なりペナルティ）
            if (P < 1e-15) {
                fR += 1e30; 
            } else {
                fR += dv.kappa / P;
            }
        }
        
        double BER_bound = (fR >= 1e30) ? 1e30 : fR * scaling_factor;
        
        out << fixed << setprecision(2) << deg << "," << scientific << setprecision(8) << BER_bound << "\n";
            
        if (BER_bound < min_BER) {
            min_BER = BER_bound;
            opt_angle = deg;
        }
    }
    
    out.close();
    cout << "Finished K=" << K << " " << name << ". Optimal Angle: " << fixed << setprecision(2) << opt_angle << " deg" << endl;
}

int main() {
    // 任意のK次元を指定して実行
    int K = 3; 
    
    run_simulation(K, 2, "BPSK");
    run_simulation(K, 2, "4QAM");
    run_simulation(K, 4, "16QAM");
    run_simulation(K, 8, "64QAM");
    run_simulation(K, 16, "256QAM");
    
    return 0;
}
