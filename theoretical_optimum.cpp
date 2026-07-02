#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include <algorithm>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Helper function: base-2 popcount
int popcount(int v) {
    int count = 0;
    while (v) {
        count += v & 1;
        v >>= 1;
    }
    return count;
}

// Recursive Gray code for arbitrary M (must be power of 2)
int gray_code(int m, int M) {
    if (M == 2) return m;
    int half = M / 2;
    if (m < half) {
        return gray_code(m, half);
    } else {
        return half + gray_code(M - 1 - m, half);
    }
}

// Calculate total bit differences for one dimension, over all pairs that differ by v_k
double D_vals(int vk, int M) {
    double sum = 0;
    for (int m = 0; m < M; ++m) {
        int m_hat = m + vk;
        if (m_hat >= 0 && m_hat < M) {
            sum += popcount(gray_code(m, M) ^ gray_code(m_hat, M));
        }
    }
    return sum;
}

// Number of valid pairs that differ by v_k
double N_vals(int vk, int M) {
    return M - std::abs(vk);
}

struct VData {
    vector<double> v;
    double weight;
};

void generate_v(int k, int K, int M, vector<int>& current_v, vector<VData>& v_list) {
    if (k == K) {
        bool all_zero = true;
        int first_non_zero = 0;
        for (int x : current_v) {
            if (x != 0) {
                if (all_zero) first_non_zero = x;
                all_zero = false;
            }
        }
        if (all_zero || first_non_zero < 0) return;

        double W = 0;
        for (int d = 0; d < K; ++d) {
            double term = D_vals(current_v[d], M);
            for (int j = 0; j < K; ++j) {
                if (j != d) term *= N_vals(current_v[j], M);
            }
            W += term;
        }
        W /= (2.0 * K * log2(M));
        W *= 2.0;

        vector<double> vd(K);
        for(int i=0; i<K; ++i) vd[i] = current_v[i];
        
        v_list.push_back({vd, W});
        return;
    }
    for (int vk = -(M - 1); vk <= M - 1; ++vk) {
        current_v.push_back(vk);
        generate_v(k + 1, K, M, current_v, v_list);
        current_v.pop_back();
    }
}

// Generate rotation matrix depending on K and theta (in radians)
vector<vector<double>> get_rotation_matrix(int K, double th) {
    vector<vector<double>> R(K, vector<double>(K, 0.0));
    for (int i = 0; i < K; ++i) R[i][i] = 1.0;

    if (K == 1) {
        // No rotation
    } else if (K == 2) {
        R[0][0] = cos(th); R[0][1] = -sin(th);
        R[1][0] = sin(th); R[1][1] = cos(th);
    } else if (K == 3) {
        double c = cos(th), s = sin(th);
        vector<vector<double>> Ry = {{c, 0, s}, {0, 1, 0}, {-s, 0, c}};
        vector<vector<double>> Rx = {{1, 0, 0}, {0, c, -s}, {0, s, c}};
        vector<vector<double>> Rz = {{c, -s, 0}, {s, c, 0}, {0, 0, 1}};
        
        // R = Ry * Rx * Rz
        vector<vector<double>> R_temp(K, vector<double>(K, 0.0));
        for(int i=0; i<3; ++i)
            for(int j=0; j<3; ++j)
                for(int k=0; k<3; ++k)
                    R_temp[i][j] += Rx[i][k] * Rz[k][j];
        
        vector<vector<double>> R_final(K, vector<double>(K, 0.0));
        for(int i=0; i<3; ++i)
            for(int j=0; j<3; ++j)
                for(int k=0; k<3; ++k)
                    R_final[i][j] += Ry[i][k] * R_temp[k][j];
        R = R_final;
    } else if (K == 4) {
        for (int i = 0; i < K; ++i) {
            for (int j = i + 1; j < K; ++j) {
                vector<vector<double>> G(K, vector<double>(K, 0.0));
                for (int k = 0; k < K; ++k) G[k][k] = 1.0;
                G[i][i] = cos(th);
                G[j][j] = cos(th);
                G[i][j] = -sin(th);
                G[j][i] = sin(th);
                
                vector<vector<double>> nextR(K, vector<double>(K, 0.0));
                for (int r = 0; r < K; ++r) {
                    for (int c = 0; c < K; ++c) {
                        for (int k = 0; k < K; ++k) {
                            nextR[r][c] += R[r][k] * G[k][c];
                        }
                    }
                }
                R = nextR;
            }
        }
    }
    return R;
}

double nCr(int n, int r) {
    if (r > n) return 0;
    if (r * 2 > n) r = n - r;
    if (r == 0) return 1;
    double result = 1;
    for (int i = 1; i <= r; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

int main() {
    int modulations[] = {2, 4, 16, 64, 256};
    int pam_sizes[] = {2, 2, 4, 8, 16};
    
    double snr_db = 30.0;
    double snr_linear = pow(10.0, snr_db / 10.0);
    
    for (int K = 1; K <= 4; ++K) {
        for (int m_idx = 0; m_idx < 5; ++m_idx) {
            int QAM = modulations[m_idx];
            int M = pam_sizes[m_idx];
            
            string filename = "C:\\Users\\Ide Nanako\\Desktop\\result3\\BER_vs_Angle_K" + to_string(K) + "_" + (QAM==2 ? "BPSK" : to_string(QAM) + "QAM") + ".csv";
            ofstream out(filename);
            if (!out.is_open()) {
                cerr << "Failed to open " << filename << endl;
                continue;
            }
            out << "Angle(deg),Theoretical_BER\n";
            
            vector<VData> v_list;
            vector<int> current_v;
            generate_v(0, K, M, current_v, v_list);
            
            double coeff = (2.0 / (pow(M, K) * log2(M)));
            double term2 = pow( ((M*M - 1.0) * 1.0) / (12.0 * log2(M) * snr_linear), K );
            double comb = nCr(2 * K - 1, K);
            double C = coeff * term2 * comb;
            
            for (double angle_deg = 0.0; angle_deg <= 45.0001; angle_deg += 0.01) {
                double th = angle_deg * M_PI / 180.0;
                vector<vector<double>> R = get_rotation_matrix(K, th);
                
                double fR = 0;
                for (const auto& vd : v_list) {
                    double P = 1.0;
                    for (int k = 0; k < K; ++k) {
                        double dot = 0;
                        for (int d = 0; d < K; ++d) {
                            dot += R[d][k] * vd.v[d];
                        }
                        double denom = dot * dot;
                        if (denom < 1e-15) {
                            P = 1e30; 
                            break;
                        }
                        P /= denom;
                    }
                    fR += vd.weight * P;
                }
                
                double BER = C * fR;
                out << fixed << setprecision(2) << angle_deg << "," << scientific << setprecision(8) << BER << "\n";
            }
            out.close();
            cout << "Finished " << filename << endl;
        }
    }
    return 0;
}
