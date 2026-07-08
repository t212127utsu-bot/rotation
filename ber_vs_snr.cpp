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
        if (all_zero) return; // 原点は除外
        
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
    for (int val = -(M - 1); val <= M - 1; ++val) {
        current_v[dim] = val;
        build_vectors(dim + 1, K, M, N_vals, D_vals, current_v, valid_vectors);
    }
}

// 任意のK次元に対するギブンス回転行列の積を生成
vector<vector<double>> get_rotation_matrix(int K, double theta) {
    vector<vector<double>> R(K, vector<double>(K, 0.0));
    for(int i=0; i<K; ++i) R[i][i] = 1.0;
    
    double c = cos(theta);
    double s = sin(theta);
    
    for(int i=0; i<K-1; ++i) {
        for(int j=i+1; j<K; ++j) {
            vector<vector<double>> G(K, vector<double>(K, 0.0));
            for(int k=0; k<K; ++k) G[k][k] = 1.0;
            G[i][i] = c;
            G[i][j] = -s;
            G[j][i] = s;
            G[j][j] = c;
            
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

// SNRに対するBERを計算する関数
void run_snr_simulation(int K, int M, string name, double opt_angle_deg) {
    cout << "Starting SNR sweep for K=" << K << ", " << name << " at optimal angle " << opt_angle_deg << " deg..." << endl;
    
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
    
    vector<DiffVector> valid_vectors;
    vector<int> current_v(K, 0);
    build_vectors(0, K, M, N_vals, D_vals, current_v, valid_vectors);
    
    string filepath = "C:\\Users\\Ide Nanako\\Desktop\\result3\\ber_eb\\SNR_vs_BER_K" + to_string(K) + "_" + name + ".csv";
    ofstream out(filepath);
    out << "SNR(dB),BER\n";
    
    // 回転行列の準備
    double th_ssd = opt_angle_deg * M_PI / 180.0;
    vector<vector<double>> R_ssd = get_rotation_matrix(K, th_ssd);
    
    // SNRを0dBから40dBまで1dB刻みでループ
    for (double snr_db = 0; snr_db <= 40.0; snr_db += 1.0) {
        double snr_linear = pow(10.0, snr_db / 10.0);
        double c_val = (M * M - 1.0) / (3.0 * log2(M) * snr_linear);
        
        double exact_ber_sum_ssd = 0;
        
        for (const auto& dv : valid_vectors) {
            // SSD (Rotated)
            vector<double> G2_ssd(K, 0.0);
            for (int k = 0; k < K; ++k) {
                double rk = 0;
                for (int j = 0; j < K; ++j) {
                    rk += R_ssd[k][j] * dv.v[j];
                }
                G2_ssd[k] = rk * rk;
            }
            
            int steps = 40;
            double dphi = (M_PI / 2.0) / steps;
            double integral_ssd = 0.0;
            for (int i = 1; i <= steps; ++i) {
                double phi = (i - 0.5) * dphi;
                double sin2_phi = sin(phi) * sin(phi);
                
                double prod_ssd = 1.0;
                for (int k = 0; k < K; ++k) {
                    prod_ssd *= sin2_phi / (sin2_phi + G2_ssd[k] / c_val);
                }
                integral_ssd += prod_ssd;
            }
            
            exact_ber_sum_ssd += dv.kappa * (integral_ssd * dphi / M_PI);
        }
        
        double BER_bound_ssd = exact_ber_sum_ssd / pow(M, K);
        out << fixed << setprecision(1) << snr_db << "," 
            << scientific << setprecision(8) << BER_bound_ssd << "\n";
    }
    
    out.close();
    cout << "Finished SNR sweep for " << name << ". Results saved." << endl;
}

int main() {
    // 実行したいKをここで指定
    int K = 2;
    
    // 以前導出した最適角をセット
    if (K == 2) {
        run_snr_simulation(K, 2, "BPSK", 27.37);
        run_snr_simulation(K, 4, "16QAM", 21.01);
        run_snr_simulation(K, 8, "64QAM", 15.69);
        run_snr_simulation(K, 16, "256QAM", 15.47);
    } else if (K == 3) {
        run_snr_simulation(K, 2, "BPSK", 30.25);
        run_snr_simulation(K, 4, "16QAM", 42.55);
        run_snr_simulation(K, 8, "64QAM", 22.38);
        run_snr_simulation(K, 16, "256QAM", 22.30);
    }
    
    return 0;
}
