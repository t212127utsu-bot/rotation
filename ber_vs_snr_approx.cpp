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

// SNRに対する【近似BER（式70）】を計算する関数
void run_snr_approx_simulation(int K, int M, string name, double opt_angle_deg) {
    cout << "Starting Approx SNR sweep for K=" << K << ", " << name << " at angle " << opt_angle_deg << " deg..." << endl;
    
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
    
    char buf[256];
    sprintf(buf, "C:\\Users\\Ide Nanako\\Desktop\\result3\\appro\\SNR_vs_BER_Approx_K%d_%s_%.2fdeg.csv", K, name.c_str(), opt_angle_deg);
    string filepath(buf);
    ofstream out(filepath);
    out << "SNR(dB),BER_Approx\n";
    
    // 回転行列の準備
    double th_ssd = opt_angle_deg * M_PI / 180.0;
    vector<vector<double>> R_ssd = get_rotation_matrix(K, th_ssd);
    
    // 近似式（積分の漸近解）における係数 C_K を計算 = (2K-1)!! / (2 * (2K)!!)
    double C_K = 1.0;
    for(int i = 1; i <= K; ++i) {
        C_K *= (double)(2 * i - 1) / (2 * i);
    }
    C_K *= 0.5;
    
    // SNRを0dBから40dBまで1dB刻みでループ
    for (double snr_db = 0; snr_db <= 40.0; snr_db += 1.0) {
        double snr_linear = pow(10.0, snr_db / 10.0);
        double c_val = (M * M - 1.0) / (3.0 * log2(M) * snr_linear);
        
        double approx_ber_sum = 0;
        
        for (const auto& dv : valid_vectors) {
            double P = 1.0;
            for (int k = 0; k < K; ++k) {
                double rk = 0;
                for (int j = 0; j < K; ++j) {
                    rk += R_ssd[k][j] * dv.v[j];
                }
                P *= (rk * rk);
            }
            
            // 式(70)に基づく高SNR近似式の計算
            // 距離が0（Pが極小）になる特異点では近似式は発散するため除外処理
            if (P > 1e-15) {
                // 【修正】c_valはSNRの逆数に比例するため、乗算すべきは(c_val)^Kです
                double term = (1.0 / P) * pow(c_val, K);
                approx_ber_sum += dv.kappa * term * C_K;
            } else {
                approx_ber_sum += 1e30; // 発散
            }
        }
        
        double BER_approx = approx_ber_sum / pow(M, K);
        out << fixed << setprecision(1) << snr_db << "," 
            << scientific << setprecision(8) << BER_approx << "\n";
    }
    
    out.close();
    cout << "Finished Approx SNR sweep for " << name << ". Results saved." << endl;
}

int main() {
    int K = 2;
    
    // 比較用に厳密式で出したのと同じ角度の近似値を算出
    if (K == 2) {
        // 真の最適角 (タイポ修正)
        run_snr_approx_simulation(K, 2, "BPSK", 27.37);
        run_snr_approx_simulation(K, 4, "16QAM", 21.06);
        run_snr_approx_simulation(K, 8, "64QAM", 15.49);
        run_snr_approx_simulation(K, 16, "256QAM", 15.48);
        
        // 従来良いとされていた20.7度
        run_snr_approx_simulation(K, 8, "64QAM", 20.7);
        run_snr_approx_simulation(K, 16, "256QAM", 20.7);
    }
    
    return 0;
}
