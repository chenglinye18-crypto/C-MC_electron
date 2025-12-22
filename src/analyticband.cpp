#include "Band.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include "mcmodel.h"

// branch index helpers for phonon table
static const int PH_LA = 0;
static const int PH_TA = 1;
static const int PH_LO = 2;
static const int PH_TO = 3;

static const double Q_SI = 1.602176634e-19;
static const double HBAR_SI = 1.054571817e-34;
static const double KB_SI = 1.380649e-23;
static const double M0_SI = 9.1093837015e-31;
static const double PI_SI = 3.1415926535897932;

// -----------------------------------------------------------------------------
// 辅助工具函数
// -----------------------------------------------------------------------------

// 非均匀 K 轴刻度生成（解析能带）
std::vector<double> Band::GenerateNonUniformTicks() {
    std::vector<double> ticks;
    std::vector<double> raw_ticks;

    // ================== 参数配置 ==================
    double k_max = 2.15;        // 覆盖布里渊区边界（单位：pi/a）
    double step_coarse = 0.05;  // 粗网格步长 (背景)

    // ================== 分支逻辑 ==================
    if (this->igzofl) {
        // ##################################################
        // #                 IGZO 模式                      #
        // ##################################################
        // 特点：单能谷，位于 Gamma 点 (0,0,0)
        // 策略：只在 0 附近加密

        double step_gamma = 0.002;  // 极细步长
        double width_gamma = 0.15;  // 加密范围 +/- 0.15

        // 1. Gamma 点附近极细网格
        for (double k = 0.5 * step_gamma; k <= width_gamma; k += step_gamma) {
            raw_ticks.push_back(k);
            raw_ticks.push_back(-k);
        }

        // 2. 背景粗网格 (全范围)
        for (double k = 0.5 * step_coarse; k <= k_max; k += step_coarse) {
            raw_ticks.push_back(k);
            raw_ticks.push_back(-k);
        }

        if (mpi_rank == 0)
            cout << "   [BandGrid] Generating IGZO Mesh (Gamma Focused, Step=0.002)" << endl;

    } else {
        // ##################################################
        // #              Silicon (Si) 模式                 #
        // ##################################################
        // 特点：多能谷，分别位于 0 (Transverse) 和 1.7 (Longitudinal)

        double k_valley_loc = 1.7;
        double step_gamma = 0.002;  // Gamma 点附近
        double width_gamma = 0.06;

        double step_valley = 0.005; // Valley 附近
        double width_valley = 0.25;

        // 1. Gamma 点附近
        for (double k = 0.5 * step_gamma; k <= width_gamma; k += step_gamma) {
            raw_ticks.push_back(k);
            raw_ticks.push_back(-k);
        }

        // 2. Valley 点 (+/- 1.7) 附近
        for (double k = 0.5 * step_valley; k <= width_valley; k += step_valley) {
            raw_ticks.push_back(k_valley_loc + k);
            raw_ticks.push_back(k_valley_loc - k);
            raw_ticks.push_back(-k_valley_loc + k);
            raw_ticks.push_back(-k_valley_loc - k);
        }

        // 3. 背景粗网格
        for (double k = 0.5 * step_coarse; k <= k_max; k += step_coarse) {
            raw_ticks.push_back(k);
            raw_ticks.push_back(-k);
        }
    }

    // ================== 后处理 (排序去重) ==================
    std::sort(raw_ticks.begin(), raw_ticks.end());

    if (raw_ticks.empty()) return ticks;

    // 压入第一个点 (需在范围 -k_max 内)
    if (raw_ticks[0] >= -k_max) {
        ticks.push_back(raw_ticks[0]);
    }

    // 过滤过近的点 (防止粗细网格重叠导致极小步长)
    double min_separation = 0.0008; // 0.4 * 0.002

    for (size_t i = 1; i < raw_ticks.size(); ++i) {
        double curr = raw_ticks[i];
        if (curr < -k_max || curr > k_max) continue;

        if (ticks.empty()) {
            ticks.push_back(curr);
        } else {
            if (curr - ticks.back() > min_separation) {
                ticks.push_back(curr);
            }
        }
    }

    return ticks;
}

// 计算解析 DOS (物理单位: 1/(eV * m^3))
// E: 能量 (eV)
// alpha: 非抛物线性因子 (1/eV)
// ml, mt: 相对有效质量 (无量纲，需在函数内乘 m0)
double CalculateAnalyticDOS_Real(double E_eV, double alpha_eV, double ml_rel, double mt_rel, int N_valley) {
    // 物理常数 (SI)
    double hbar = 1.0545718e-34; // J*s
    double q    = 1.6021766e-19; // C (J/eV)
    double m0   = 9.1093835e-31; // kg
    
    // 转换单位到 SI (Joule)
    double E_J = E_eV * q;
    double alpha_J = alpha_eV / q; 
    
    // Gamma(E) = E(1 + alpha*E)
    double gamma = E_J * (1.0 + alpha_J * E_J);
    double gamma_prime = 1.0 + 2.0 * alpha_J * E_J; // d(Gamma)/dE
    
    if (gamma < 0) return 0.0;
    
    // 有效质量项 (m_dos)^(3/2) = sqrt(ml) * mt * m0^(3/2)
    // 注意：根据 DOS 公式 g(E) ~ (m_dos)^1.5 * ...
    // m_dos = (ml * mt^2)^(1/3)
    // (m_dos)^(3/2) = (ml * mt^2)^(1/2) = sqrt(ml) * mt
    double mass_factor = std::sqrt(ml_rel * m0) * (mt_rel * m0); 
    
    // DOS = N_v * sqrt(2)/pi^2 * 1/hbar^3 * (m_dos)^3/2 * sqrt(gamma) * gamma'
    double prefactor = N_valley * (std::sqrt(2.0) / (PI * PI * std::pow(hbar, 3)));
    
    double dos_J = prefactor * mass_factor * std::sqrt(gamma) * gamma_prime;
    
    // 转换回 1/(eV * m^3)
    return dos_J * q; 
}

// -----------------------------------------------------------------------------
// 核心初始化函数
// -----------------------------------------------------------------------------

void Band::InitAnalyticBand(double alpha_norm, double ml_rel, double mt_rel, string input_path) {
    cout << "Initializing Analytic Band (Kane's Model - Vector V)..." << endl;

    // 说明（单位一致性）：
    // - 全局 sia0 在 init_phpysical_parameter 中已做归一化：a_code = a_real / spr0
    // - 本函数需要把 ticks(pi/a) 转为 1/m 时，使用 (PI/sia0)/spr0

    // 1. 物理常数 (SI)
    const double m0 = M0_SI;
    const double hbar = HBAR_SI;
    const double q = Q_SI;
    const double a_lattice_code = sia0;

    // 2. 参数：IGZO 设为抛物线 (alpha=0)，Si 使用传入的 alpha_norm
    const double alpha_real = this->igzofl ? 0.0 : (alpha_norm / eV0); // 1/eV
    this->analytic_alpha_real = alpha_real;
    const double ml_kg = ml_rel * m0;
    const double mt_kg = mt_rel * m0;

    // 3. 文件后缀（IGZO 写 *_IGZO.txt，Si 保持原名）
    const string file_suffix = this->igzofl ? "_IGZO" : "";

    // ----------------------------------------------------------
    // 4. 波谷配置
    // ----------------------------------------------------------
    const int num_valleys = this->igzofl ? 1 : 6;
    const double K_valley_norm = this->igzofl ? 0.0 : 1.7;
    double valley_dirs[6][3] = {
        {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}
    };
    int valley_l_axis[6] = {0, 0, 1, 1, 2, 2};
    if (this->igzofl) {
        valley_dirs[0][0] = 0.0; valley_dirs[0][1] = 0.0; valley_dirs[0][2] = 0.0;
        valley_l_axis[0] = 0;
    }

    // ----------------------------------------------------------
    // 5. 生成 E-k-v 表
    // ----------------------------------------------------------
    const string ek_file = input_path + "/analytic_ek" + file_suffix + ".txt";
    ofstream out_ek(ek_file.c_str());
    if (!out_ek.is_open()) {
        cout << "Error: Cannot open " << ek_file << " for writing!" << endl;
        exit(1);
    }
    out_ek << "kx(pi/a) ky(pi/a) kz(pi/a) Energy(eV) vx(m/s) vy(m/s) vz(m/s)" << endl;

    const vector<double> ticks = GenerateNonUniformTicks();
    const int num_ticks = static_cast<int>(ticks.size());
    cout << "  Generating E-k-v table (" << (this->igzofl ? "IGZO" : "Si") << ") with "
         << num_ticks << "^3 points..." << endl;

    const double k_conversion = (PI / a_lattice_code) / spr0; // pi/a -> 1/m
    const double J_to_eV = 1.0 / q;

    for (int i = 0; i < num_ticks; i++) {
        for (int j = 0; j < num_ticks; j++) {
            for (int k = 0; k < num_ticks; k++) {
                const double gx = ticks[i];
                const double gy = ticks[j];
                const double gz = ticks[k];

                // 1) 寻找最近波谷
                int best_valley = 0;
                double min_dist_sq = 1.0e99;
                double dx_best = 0.0, dy_best = 0.0, dz_best = 0.0;
                for (int v = 0; v < num_valleys; v++) {
                    const double cx = valley_dirs[v][0] * K_valley_norm;
                    const double cy = valley_dirs[v][1] * K_valley_norm;
                    const double cz = valley_dirs[v][2] * K_valley_norm;
                    const double dist_sq = (gx-cx)*(gx-cx) + (gy-cy)*(gy-cy) + (gz-cz)*(gz-cz);
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        best_valley = v;
                        dx_best = gx - cx;
                        dy_best = gy - cy;
                        dz_best = gz - cz;
                    }
                }

                // 2) 局部坐标 (SI: 1/m)
                const double dkx_real = dx_best * k_conversion;
                const double dky_real = dy_best * k_conversion;
                const double dkz_real = dz_best * k_conversion;

                double kl = 0.0, kt1 = 0.0, kt2 = 0.0;
                const int l_axis = valley_l_axis[best_valley];
                if (l_axis == 0) { // X
                    kl = dkx_real; kt1 = dky_real; kt2 = dkz_real;
                } else if (l_axis == 1) { // Y
                    kl = dky_real; kt1 = dkx_real; kt2 = dkz_real;
                } else { // Z
                    kl = dkz_real; kt1 = dkx_real; kt2 = dky_real;
                }

                // 3) 能量 (Kane 模型：Gamma(E)=E(1+alphaE)=hbar^2 k^2 / 2m)
                const double Gamma_J = (hbar*hbar/2.0) * ( (kl*kl)/ml_kg + (kt1*kt1 + kt2*kt2)/mt_kg );
                const double Gamma_eV = Gamma_J * J_to_eV;

                double E_val = 0.0;
                if (alpha_real > 1e-12) {
                    E_val = (-1.0 + std::sqrt(1.0 + 4.0 * alpha_real * Gamma_eV)) / (2.0 * alpha_real);
                } else {
                    E_val = Gamma_eV;
                }

                // 4) 速度矢量：v = (1/hbar) * dE/dk = (hbar * k / m) / (1 + 2 alpha E)
                const double dGamma_dE = 1.0 + 2.0 * alpha_real * E_val;
                const double v_prefactor = hbar / dGamma_dE;
                const double vl_vel  = v_prefactor * (kl  / ml_kg);
                const double vt1_vel = v_prefactor * (kt1 / mt_kg);
                const double vt2_vel = v_prefactor * (kt2 / mt_kg);

                double vx = 0.0, vy = 0.0, vz = 0.0;
                if (l_axis == 0) {
                    vx = vl_vel; vy = vt1_vel; vz = vt2_vel;
                } else if (l_axis == 1) {
                    vx = vt1_vel; vy = vl_vel; vz = vt2_vel;
                } else {
                    vx = vt1_vel; vy = vt2_vel; vz = vl_vel;
                }

                out_ek << gx << " " << gy << " " << gz << " "
                       << E_val << " " << vx << " " << vy << " " << vz << endl;
            }
        }
    }
    out_ek.close();
    cout << "  E-k-v table generated: " << ek_file << endl;

    // ----------------------------------------------------------
    // 6. 生成 DOS 表并填充内部数组
    // ----------------------------------------------------------
    const string dos_file = input_path + "/analytic_dos" + file_suffix + ".txt";
    ofstream out_dos(dos_file.c_str());
    if (!out_dos.is_open()) {
        cout << "Error: Cannot open " << dos_file << " for writing!" << endl;
        exit(1);
    }
    out_dos << "Energy(eV) DOS(1/eV/m^3) DOS_Norm(CodeUnits)" << endl;

    const double max_E_eV = 7.5;
    const int N_valley = num_valleys;

    // Kane DOS：g(E) = Nv/(2pi^2) * (2 m_dos / hbar^2)^(3/2) * sqrt(gamma_J) * (1+2 alpha E)
    const double md_dos = std::pow(ml_kg * mt_kg * mt_kg, 1.0/3.0);
    const double dos_prefactor = (static_cast<double>(N_valley) / (2.0 * PI_SI * PI_SI))
                               * std::pow(2.0 * md_dos / (HBAR_SI * HBAR_SI), 1.5);

    for (int itab = 0; itab <= MTAB; itab++) {
        const double E_eV = energy[itab] * eV0;
        if (E_eV > max_E_eV) break;

        double dos_real = 0.0; // 1/(eV*m^3)
        if (E_eV >= 0.0) {
            const double gamma_J = (E_eV * Q_SI) * (1.0 + alpha_real * E_eV);
            if (gamma_J >= 0.0) {
                const double gamma_prime = 1.0 + 2.0 * alpha_real * E_eV;
                const double dos_per_Joule = dos_prefactor * std::sqrt(gamma_J) * gamma_prime; // 1/(J*m^3)
                dos_real = dos_per_Joule * Q_SI; // 1/(eV*m^3)
            }
        }

        const double dos_code = dos_real * eV0 * std::pow(spr0, 3.0);
        dos[bandof[PELEC]][itab] = dos_code;
        sumdos[itab][PELEC] = dos_code;
        out_dos << E_eV << " " << dos_real << " " << dos_code << endl;
    }
    out_dos.close();

    DOSMAX[PELEC] = 0.0;
    for (int itab = 0; itab <= MTAB; itab++) {
        if (sumdos[itab][PELEC] > DOSMAX[PELEC]) DOSMAX[PELEC] = sumdos[itab][PELEC];
    }

    cout << "  DOS table generated: " << dos_file << endl;
}

// -----------------------------------------------------------------------------
// 解析能带读取与建表模块
// -----------------------------------------------------------------------------

// 辅助：根据 k 值确定生成时的步长
static double GetGridStep(double k_pi_val) {
    double valley_loc = 1.7; 
    double fine_region_width = 0.15;
    double step_coarse = 0.1;
    double step_fine = 0.01;
    
    if (std::fabs(std::fabs(k_pi_val) - valley_loc) < fine_region_width) {
        return step_fine;
    }
    return step_coarse;
}

void Band::BuildAnalyticLists(string pathname) {
    cout << "Building Analytic Lists (Indexing)..." << endl;
    if (dlist <= 0) dlist = eV0 / 0.002 ;

    analytic_ntlist.assign(MWLE_ana, 0);
    analytic_ptlist.assign(MWLE_ana, 0);
    analytic_tlist.resize(analytic_k_grid.size());

    for (size_t i = 0; i < analytic_k_grid.size(); ++i) {
        double E = analytic_k_grid[i].energy;
        int itab = (int)((E - emin) * dlist);
        if (itab >= 0 && itab < MWLE_ana) analytic_ntlist[itab]++;
    }

    int current_offset = 0;
    for (int itab = 0; itab < MWLE_ana; ++itab) {
        analytic_ptlist[itab] = current_offset;
        current_offset += analytic_ntlist[itab];
    }

    std::vector<int> temp_counters(MWLE_ana, 0);
    for (size_t i = 0; i < analytic_k_grid.size(); ++i) {
        double E = analytic_k_grid[i].energy;
        int itab = (int)((E - emin) * dlist);
        if (itab >= 0 && itab < MWLE_ana) {
            int pos = analytic_ptlist[itab] + temp_counters[itab];
            analytic_tlist[pos] = (int)i;
            temp_counters[itab]++;
        }
    }
    cout << "  Analytic lists built. Indexed " << analytic_k_grid.size() << " states into " << MWLE_ana << " energy bins." << endl;

    string debug_file = pathname + "/debug_analytic_bins.txt";
    ofstream out(debug_file.c_str());
    out << "Bin_Index Energy_Min(eV) Energy_Max(eV) Count Start_Index" << endl;

    int empty_bins = 0;
    int total_bins = MWLE_ana; // 或者是实际用到的最大 itab

    for (int i = 0; i < total_bins; ++i) {
        int count = analytic_ntlist[i];
        if (count == 0) empty_bins++;
        
        // 仅打印有数据的 Bin，或者前 100 个 Bin 用于检查
        if (count > 0 || i < MWLE_ana) {
            double E_min_eV = (i * (1.0/dlist) + emin) * eV0;
            double E_max_eV = ((i+1) * (1.0/dlist) + emin) * eV0;
            out << i << " " << E_min_eV << " " << E_max_eV << " " 
                << count << " " << analytic_ptlist[i] << endl;
        }
    }
    out.close();
    cout << "  [Debug] Bin statistics dumped to: " << debug_file << endl;
    cout << "  [Debug] Total Bins: " << total_bins << ", Empty Bins: " << empty_bins << endl;
}

void Band::InitPhononSpectrum(string input_path) {
    string filename = input_path + (this->igzofl ? "/phonon_dispersion_IGZO.txt" : "/phonon_dispersion.txt");
    cout << "Reading Phonon Spectrum from: " << filename << endl;

    // reset header-derived scalars (vectors cleared below)
    phonon.a0 = 0.0;
    phonon.qmax = 0.0;
    phonon.dq = 0.0;
    phonon.nq_tab = 0;

    ifstream infile(filename.c_str());
    if (!infile) {
        cerr << "Error: Cannot open phonon dispersion file!" << endl;
        exit(1);
    }

    string line;
    double max_q_found = 0.0;
    // 1. 跳过表头并解析 a0、qmax（如果表头提供）
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] != '#') {
            // 数据开始
            break;
        }
        if (line.find("a0=") != string::npos) {
            size_t pos_a0 = line.find("a0=");
            size_t pos_qmax = line.find("qmax=");
            if (pos_a0 != string::npos) {
                phonon.a0 = std::atof(line.substr(pos_a0 + 3).c_str());
            }
            if (pos_qmax != string::npos) {
                phonon.qmax = std::atof(line.substr(pos_qmax + 5).c_str());
            }
        }
    }

    // 2. 清空表格
    int num_branches = 4;
    for(int i=0; i<num_branches; ++i) {
        phonon.omega_table[i].clear();
        phonon.vg_table[i].clear();
    }

    // 3. 解析数据行（当前 line 已经是第一行数据）
    // 格式: q | w_LA w_TA w_LO w_TO | v_LA v_TA v_LO v_TO
    do {
        if (line.empty() || line[0] == '#') continue;

        stringstream ss(line);
        double q_val;
        double w[4], v[4];

        ss >> q_val;
        for(int i=0; i<4; ++i) ss >> w[i];
        for(int i=0; i<4; ++i) ss >> v[i];

        // 允许存在非数值表头行（例如 "q omega..."），解析失败则跳过
        if (ss.fail()) continue;
        if (!std::isfinite(q_val)) continue;
        bool ok = true;
        for (int i = 0; i < 4; ++i) {
            if (!std::isfinite(w[i]) || !std::isfinite(v[i])) { ok = false; break; }
        }
        if (!ok) continue;
        if (q_val > max_q_found) max_q_found = q_val;

        for(int i=0; i<4; ++i) {
            phonon.omega_table[i].push_back(w[i]);
            phonon.vg_table[i].push_back(v[i]);
        }

    } while (std::getline(infile, line));

    infile.close();

    phonon.nq_tab = phonon.omega_table[0].size();
    // fallback: if header didn't provide qmax, infer from data
    if (phonon.qmax <= 0.0 && max_q_found > 0.0) phonon.qmax = max_q_found;
    if (phonon.nq_tab > 1) {
        phonon.dq = phonon.qmax / (phonon.nq_tab - 1);
    }

    cout << "  Loaded " << phonon.nq_tab << " points for phonon dispersion." << endl;
    cout << "  a0 = " << phonon.a0 << ", qmax = " << phonon.qmax << endl;
}

// ----------------------------------------------------------------------------- 
// 解析散射率：辅助函数
// -----------------------------------------------------------------------------

double Band::GetPhononOmega(int branch, double q) {
    if (phonon.nq_tab < 2 || branch < 0 || branch > 3) return 0.0;

    double dq = phonon.dq;
    if (dq <= 0 && phonon.qmax > 0 && phonon.nq_tab > 1) {
        dq = phonon.qmax / (phonon.nq_tab - 1);
    }
    if (dq <= 0) return 0.0;

    if (q <= 0) return phonon.omega_table[branch].front();
    if (q >= phonon.qmax) return phonon.omega_table[branch].back();

    int idx = static_cast<int>(q / dq);
    if (idx >= phonon.nq_tab - 1) idx = phonon.nq_tab - 2;
    double q1 = idx * dq;
    double t = (q - q1) / dq;
    double w1 = phonon.omega_table[branch][idx];
    double w2 = phonon.omega_table[branch][idx + 1];
    return w1 + (w2 - w1) * t;
}

double Band::GetKaneK_SI(double E_eV) {
    double ml_rel = 0.916;
    double mt_rel = 0.190;
    double ml = ml_rel * M0_SI;
    double mt = mt_rel * M0_SI;
    double md = std::pow(ml * mt * mt, 1.0/3.0);

    double alpha = 0.5;
    double term = E_eV * (1.0 + alpha * E_eV);
    if (term < 0) return 0.0;
    return std::sqrt(2.0 * md * term * Q_SI) / HBAR_SI;
}

double Band::GetKaneDOS_SI(double E_eV) {
    if (E_eV < 0) return 0.0;
    double ml = 0.916 * M0_SI;
    double mt = 0.190 * M0_SI;
    double alpha = 0.5;

    double E_J = E_eV * Q_SI;
    double alpha_J = alpha / Q_SI;
    double gamma = E_J * (1.0 + alpha_J * E_J);
    double gamma_prime = 1.0 + 2.0 * alpha_J * E_J;
    if (gamma < 0) return 0.0;

    double md = std::pow(ml * mt * mt, 1.0/3.0);
    double pre = std::sqrt(2.0) * std::pow(md, 1.5) / (PI_SI * PI_SI * std::pow(HBAR_SI, 3));
    return pre * std::sqrt(gamma) * gamma_prime;
}

double Band::GetOverlapFactor(double q, double Rs) {
    double qRs = q * Rs;
    if (std::fabs(qRs) < 1e-12) return 1.0;
    return 3.0 / (qRs * qRs * qRs) * (std::sin(qRs) - qRs * std::cos(qRs));
}

bool Band::CheckAllowedQ(double E_eV, double ks, double q, double hw_eV, int type,
                         double md_SI, double alpha_eV) {
    if (q < 1e-12) return false;
    if (ks <= 0.0) return false;

    const double hbar = HBAR_SI;
    const double qJ_per_eV = Q_SI;

    if (type == -1) { // emission: E' = E - hw
        if (E_eV - hw_eV <= 0.0) return false;
        const double num = md_SI * hw_eV * qJ_per_eV;
        const double den = (hbar * hbar) * q * ks;
        if (std::fabs(den) < 1e-60) return false;
        const double non_par = 1.0 + alpha_eV * (2.0 * E_eV - hw_eV);
        const double term = (num / den) * non_par;
        const double cos_theta = q / (2.0 * ks) + term;
        return std::fabs(cos_theta) <= 1.0;
    }

    if (type == 1) { // absorption: E' = E + hw
        const double num = md_SI * hw_eV * qJ_per_eV;
        const double den = (hbar * hbar) * q * ks;
        if (std::fabs(den) < 1e-60) return false;
        const double non_par = 1.0 + alpha_eV * (2.0 * E_eV + hw_eV);
        const double term = (num / den) * non_par;
        const double cos_theta = -q / (2.0 * ks) + term;
        return std::fabs(cos_theta) <= 1.0;
    }

    if (type == 0) { // elastic fallback
        return q <= 2.0 * ks;
    }

    return false;
}

/*
void Band::BuildAnalyticScatteringTable() {
    if (this->igzofl) {
        cout << "Building Analytic Scattering Table (IGZO: 3 Processes)..." << endl;

        const double T_lattice = T0;
        const double rho = 6100.0;

        const double E_ac_eV = 5.0; //形变势
        const double D_LA = E_ac_eV * Q_SI;
        const double D_TA = E_ac_eV * Q_SI;

        const double ml = mell * M0_SI;
        const double mt = melt * M0_SI;
        const double md = std::pow(ml * mt * mt, 1.0/3.0);

        const double a0 = phonon.a0;
        const double Rs = (a0 > 0) ? a0 * std::pow(3.0/(16.0*PI_SI), 1.0/3.0) : 0.0;
        const int nq_int = 400;
        const double dq_int = (phonon.qmax > 0 && nq_int > 1) ? phonon.qmax / (nq_int - 1) : 0.0;

        scpre = 3; // 0: acoustic, 1: optical abs, 2: optical em
        const int band_idx = bandof[PELEC];

        // 预清零
        for (int iproc = 0; iproc < scpre; ++iproc) {
            for (int ib = 0; ib < NBE; ++ib) {
                scatte[iproc][ib][ib] = 0.0;
                for (int itab = 0; itab <= MTAB; ++itab) {
                    dose[iproc][ib][itab] = 0.0;
                }
            }
        }
        for (int ib = 0; ib < NBE; ++ib) {
            for (int jb = 0; jb < NBE; ++jb) {
                for (int itab = 0; itab <= MTAB; ++itab) {
                    scattiie[ib][jb][itab] = 0.0;
                }
            }
        }

        for (int iproc = 0; iproc < scpre; ++iproc) {
            scatte[iproc][band_idx][band_idx] = 1.0;
        }

        const double alpha_real = this->analytic_alpha_real;
        const double md_dos = std::pow(ml * mt * mt, 1.0/3.0);
        const double dos_prefactor = (1.0 / (2.0 * PI_SI * PI_SI))
                                   * std::pow(2.0 * md_dos / (HBAR_SI * HBAR_SI), 1.5); // 1/(J^(3/2) m^3)

        auto get_dos_si_from_table = [&](double E_eV_query) -> double {
            if (E_eV_query < 0) return 0.0;

            // if beyond DOS table range, extrapolate using the same analytic DOS form
            const double Emax_eV = energy[MTAB] * eV0;
            if (E_eV_query > Emax_eV) {
                const double gamma_J = (E_eV_query * Q_SI) * (1.0 + (alpha_real > 0.0 ? alpha_real * E_eV_query : 0.0));
                if (gamma_J <= 0.0) return 0.0;
                const double gamma_prime = 1.0 + 2.0 * (alpha_real > 0.0 ? alpha_real * E_eV_query : 0.0);
                const double dos_per_J = dos_prefactor * std::sqrt(gamma_J) * gamma_prime; // 1/(J*m^3)
                return dos_per_J * Q_SI; // 1/(eV*m^3)
            }

            const double E_norm = E_eV_query / eV0;
            int itab_q = static_cast<int>(((E_norm - emin) / dtable) + 0.5);
            if (itab_q < 0) itab_q = 0;
            if (itab_q > MTAB) itab_q = MTAB;
            return sumdos[itab_q][PELEC] / (eV0 * std::pow(spr0, 3.0));
        };

        // 光学声子能量：取 LO 分支在 Gamma 的值；若过小则兜底为 60 meV
        double omega_LO = 0.0;
        if (phonon.nq_tab > 0) omega_LO = phonon.omega_table[PH_LO].front();
        if (!(omega_LO > 0.0) && phonon.nq_tab > 1) omega_LO = phonon.omega_table[PH_LO][1];
        if (!(omega_LO > 0.0) && phonon.nq_tab > 2) omega_LO = phonon.omega_table[PH_LO][phonon.nq_tab / 2];
        double hw_LO_eV = (omega_LO > 0.0) ? (HBAR_SI * omega_LO / Q_SI) : 0.0;
        if (hw_LO_eV < 0.02) {
            hw_LO_eV = 0.06;
            if (mpi_rank == 0) {
                cout << "  [Warning] LO energy too low in table, using default 60 meV for IGZO." << endl;
            }
        }

        const double w0_LO = hw_LO_eV * Q_SI / HBAR_SI;
        const double Nq_LO = 1.0 / (std::exp(hw_LO_eV * Q_SI / (KB_SI * T_lattice)) - 1.0);
        const double Dopt_eVm = 5e10;
        const double D_Jm = Dopt_eVm * Q_SI;
        const double C_LO = (PI_SI * D_Jm * D_Jm) / (2.0 * rho * w0_LO);

        for (int itab = 0; itab <= MTAB; ++itab) {
            const double E_eV = energy[itab] * eV0;
            sumscatt[itab][band_idx] = 0.0;

            // -------- Acoustic (elastic; phonon spectrum integral) --------
            double Rate_AC_SI = 0.0;
            double term = E_eV * (1.0 + alpha_real * E_eV);
            if (term < 0.0) term = 0.0;
            const double ks = (term > 0.0) ? std::sqrt(2.0 * md * term * Q_SI) / HBAR_SI : 0.0;

            if (ks > 1e-30 && dq_int > 0 && phonon.nq_tab > 1) {
                double integ_LA = 0.0;
                double integ_TA = 0.0;

                for (int iq = 0; iq < nq_int; ++iq) {
                    const double q = iq * dq_int;
                    if (q < 1e-12) continue;

                    const double w_LA = GetPhononOmega(PH_LA, q);
                    const double w_TA = GetPhononOmega(PH_TA, q);
                    if (w_LA <= 0.0 || w_TA <= 0.0) continue;

                    const double N_LA = 1.0 / (std::exp(HBAR_SI * w_LA / (KB_SI * T_lattice)) - 1.0);
                    const double N_TA = 1.0 / (std::exp(HBAR_SI * w_TA / (KB_SI * T_lattice)) - 1.0);

                    const double Iq = (Rs > 0) ? GetOverlapFactor(q, Rs) : 1.0;
                    const double q3_I2 = q * q * q * Iq * Iq;

                    const double hw_LA_eV = HBAR_SI * w_LA / Q_SI;
                    const double hw_TA_eV = HBAR_SI * w_TA / Q_SI;

                    // absorption
                    if (CheckAllowedQ(E_eV, ks, q, hw_LA_eV, 1, md, alpha_real)) {
                        integ_LA += (1.0 / w_LA) * N_LA * q3_I2;
                    }
                    // emission
                    if (E_eV > hw_LA_eV && CheckAllowedQ(E_eV, ks, q, hw_LA_eV, -1, md, alpha_real)) {
                        integ_LA += (1.0 / w_LA) * (N_LA + 1.0) * q3_I2;
                    }

                    if (CheckAllowedQ(E_eV, ks, q, hw_TA_eV, 1, md, alpha_real)) {
                        integ_TA += (1.0 / w_TA) * N_TA * q3_I2;
                    }
                    if (E_eV > hw_TA_eV && CheckAllowedQ(E_eV, ks, q, hw_TA_eV, -1, md, alpha_real)) {
                        integ_TA += (1.0 / w_TA) * (N_TA + 1.0) * q3_I2;
                    }
                }

                const double pre = md / (4.0 * PI_SI * rho * HBAR_SI * HBAR_SI * ks);
                Rate_AC_SI = pre * (D_LA * D_LA * integ_LA + D_TA * D_TA * integ_TA) * dq_int;
            }

            dose[0][band_idx][itab] = Rate_AC_SI * time0;
            sumscatt[itab][band_idx] += dose[0][band_idx][itab];

            // -------- Optical (equivalent POP via zero-order optical DP) --------
            const double g_abs = get_dos_si_from_table(E_eV + hw_LO_eV);
            const double Rate_Abs_SI = C_LO * Nq_LO * g_abs;
            dose[1][band_idx][itab] = Rate_Abs_SI * time0;
            sumscatt[itab][band_idx] += dose[1][band_idx][itab];

            double Rate_Em_SI = 0.0;
            if (E_eV > hw_LO_eV) {
                const double g_em = get_dos_si_from_table(E_eV - hw_LO_eV);
                Rate_Em_SI = C_LO * (Nq_LO + 1.0) * g_em;
            }
            dose[2][band_idx][itab] = Rate_Em_SI * time0;
            sumscatt[itab][band_idx] += dose[2][band_idx][itab];

            scattiie[band_idx][band_idx][itab] = 0.0;
        }

        double max_gamma = 0.0;
        for (int itab = 0; itab <= MTAB; ++itab) {
            if (sumscatt[itab][band_idx] > max_gamma) {
                max_gamma = sumscatt[itab][band_idx];
            }
        }

        if (nt <= 0) nt = 1;
        int fill_nt = nt;
        if (fill_nt > MNTet) fill_nt = MNTet;
        for(int it = 0; it < fill_nt; it++) {
            gamtet[it] = max_gamma;
        }
        gamma[PELEC] = max_gamma;

        // =========================================================
        // [DEBUG] 导出散射率与能量关系 (Export Scattering Rates)
        // =========================================================
        if (mpi_rank == 0) {
            const string base_dir = pathname.empty() ? string("input") : pathname;
            string dump_name = this->igzofl ? (base_dir + "/scattering_rates_IGZO.txt")
                                            : (base_dir + "/scattering_rates_Si.txt");
            ofstream out_scat(dump_name.c_str());
            if (!out_scat.is_open()) {
                cout << "  [Debug] Failed to open scattering dump file: " << dump_name << endl;
            } else {
                out_scat << "Energy(eV) Total(1/s)";
                for(int i=0; i<scpre; i++) {
                    out_scat << " Mech_" << i;
                }
                out_scat << endl;

                for (int itab = 0; itab <= MTAB; ++itab) {
                    double E_real = energy[itab] * eV0;
                    if (E_real > 3.0) break;

                    double rate_total = sumscatt[itab][band_idx] / time0;
                    out_scat << E_real << " " << rate_total;

                    for (int i = 0; i < scpre; ++i) {
                        double rate_proc = dose[i][band_idx][itab] / time0;
                        out_scat << " " << rate_proc;
                    }
                    out_scat << endl;
                }
                out_scat.close();
                cout << "  [Debug] Scattering rates saved to: " << dump_name << endl;
            }
        }
        // =========================================================

        cout << "  Analytic scattering table built (IGZO). Max Rate (norm) = " << max_gamma << endl;
        return;
    }

    cout << "Building Analytic Scattering Table (14 Processes)..." << endl;

    double T_lattice = T0;
    double rho = 2330.0;

    double D_LA = 6.39 * Q_SI;
    double D_TA = 3.01 * Q_SI;
    double ml = 0.916 * M0_SI;
    double mt = 0.190 * M0_SI;
    double md = std::pow(ml * mt * mt, 1.0/3.0);

    double a0 = phonon.a0;
    double Rs = (a0 > 0) ? a0 * std::pow(3.0/(16.0*PI_SI), 1.0/3.0) : 0.0;
    int nq_int = 400;
    double dq_int = (phonon.qmax > 0 && nq_int > 1) ? phonon.qmax / (nq_int - 1) : 0.0;

    struct IvParam { double E_meV; double D_1e8; int Z; };
    IvParam iv_params[6] = {
        {10.0, 0.3, 1},
        {19.0, 1.5, 1},
        {62.0, 6.0, 1},
        {19.0, 0.5, 4},
        {51.0, 3.5, 4},
        {57.0, 1.5, 4}
    };

    scpre = 14;
    int band_idx = bandof[PELEC];

    // 预清零
    for (int iproc = 0; iproc < scpre; ++iproc) {
        for (int ib = 0; ib < NBE; ++ib) {
            scatte[iproc][ib][ib] = 0.0;
            for (int itab = 0; itab <= MTAB; ++itab) {
                dose[iproc][ib][itab] = 0.0;
            }
        }
    }
    for (int ib = 0; ib < NBE; ++ib) {
        for (int jb = 0; jb < NBE; ++jb) {
            for (int itab = 0; itab <= MTAB; ++itab) {
                scattiie[ib][jb][itab] = 0.0;
            }
        }
    }

    for (int iproc = 0; iproc < scpre; ++iproc) {
        scatte[iproc][band_idx][band_idx] = 1.0;
    }

    for (int itab = 0; itab <= MTAB; ++itab) {
        double E_norm = energy[itab];
        double E_eV = E_norm * eV0;
        sumscatt[itab][band_idx] = 0.0;

        double Rate_AC_SI = 0.0;
        double ks = GetKaneK_SI(E_eV);

        if (ks > 1e-30 && dq_int > 0 && phonon.nq_tab > 1) {
            double integ_LA = 0.0;
            double integ_TA = 0.0;
            const double alpha_si = 0.5;

            for (int iq = 0; iq < nq_int; ++iq) {
                double q = iq * dq_int;
                if (q < 1e-12) continue;

                double w_LA = GetPhononOmega(PH_LA, q);
                double w_TA = GetPhononOmega(PH_TA, q);

                if (w_LA <= 0 || w_TA <= 0) continue;

                double N_LA = 1.0 / (std::exp(HBAR_SI * w_LA / (KB_SI * T_lattice)) - 1.0);
                double N_TA = 1.0 / (std::exp(HBAR_SI * w_TA / (KB_SI * T_lattice)) - 1.0);

                double Iq = (Rs > 0) ? GetOverlapFactor(q, Rs) : 1.0;
                double q3_I2 = q * q * q * Iq * Iq;

                double hw_LA_eV = HBAR_SI * w_LA / Q_SI;
                double hw_TA_eV = HBAR_SI * w_TA / Q_SI;

                if (CheckAllowedQ(E_eV, ks, q, hw_LA_eV, 1, md, alpha_si)) {
                    integ_LA += (1.0/w_LA) * N_LA * q3_I2;
                }
                if (E_eV > hw_LA_eV && CheckAllowedQ(E_eV, ks, q, hw_LA_eV, -1, md, alpha_si)) {
                    integ_LA += (1.0/w_LA) * (N_LA + 1.0) * q3_I2;
                }

                if (CheckAllowedQ(E_eV, ks, q, hw_TA_eV, 1, md, alpha_si)) {
                    integ_TA += (1.0/w_TA) * N_TA * q3_I2;
                }
                if (E_eV > hw_TA_eV && CheckAllowedQ(E_eV, ks, q, hw_TA_eV, -1, md, alpha_si)) {
                    integ_TA += (1.0/w_TA) * (N_TA + 1.0) * q3_I2;
                }
            }
            double pre = md / (4.0 * PI_SI * rho * HBAR_SI * HBAR_SI * ks);
            Rate_AC_SI = pre * (D_LA * D_LA * integ_LA + D_TA * D_TA * integ_TA) * dq_int;
        }

        dose[0][band_idx][itab] = Rate_AC_SI * time0;
        sumscatt[itab][band_idx] += dose[0][band_idx][itab];

        for (int i = 0; i < 6; ++i) {
            double hw_eV = iv_params[i].E_meV * 1e-3;
            double w0 = hw_eV * Q_SI / HBAR_SI;
            double Nq = 1.0 / (std::exp(hw_eV * Q_SI / (KB_SI * T_lattice)) - 1.0);
            double D_Jm = iv_params[i].D_1e8 * 1e10 * Q_SI;
            double C = (PI_SI * D_Jm * D_Jm * iv_params[i].Z) / (2.0 * rho * w0);

            int idx_abs = 1 + 2 * i;
            double Rate_Abs_SI = 0.0;
            double g_ab = GetKaneDOS_SI(E_eV + hw_eV);
            Rate_Abs_SI = C * Nq * g_ab;

            dose[idx_abs][band_idx][itab] = Rate_Abs_SI * time0;
            sumscatt[itab][band_idx] += dose[idx_abs][band_idx][itab];

            int idx_em = 2 + 2 * i;
            double Rate_Em_SI = 0.0;
            if (E_eV > hw_eV) {
                double g_em = GetKaneDOS_SI(E_eV - hw_eV);
                Rate_Em_SI = C * (Nq + 1.0) * g_em;
            }
            dose[idx_em][band_idx][itab] = Rate_Em_SI * time0;
            sumscatt[itab][band_idx] += dose[idx_em][band_idx][itab];
        }

        scattiie[band_idx][band_idx][itab] = 0.0;
    }

    double max_gamma = 0.0;
    for (int itab = 0; itab <= MTAB; ++itab) {
        if (sumscatt[itab][band_idx] > max_gamma) {
            max_gamma = sumscatt[itab][band_idx];
        }
    }

    // [修复] 确保 gamtet 分配并填充最大散射率
    if (nt <= 0) nt = 1;
    int fill_nt = nt;
    if (fill_nt > MNTet) fill_nt = MNTet;
    for(int it = 0; it < fill_nt; it++) {
        gamtet[it] = max_gamma;
    }
    gamma[PELEC] = max_gamma;

    // =========================================================
    // [DEBUG] 导出散射率与能量关系 (Export Scattering Rates)
    // =========================================================
    if (mpi_rank == 0) {
        const string base_dir = pathname.empty() ? string("input") : pathname;
        string dump_name = this->igzofl ? (base_dir + "/scattering_rates_IGZO.txt")
                                        : (base_dir + "/scattering_rates_Si.txt");
        ofstream out_scat(dump_name.c_str());
        if (!out_scat.is_open()) {
            cout << "  [Debug] Failed to open scattering dump file: " << dump_name << endl;
        } else {
            out_scat << "Energy(eV) Total(1/s)";
            for(int i=0; i<scpre; i++) {
                out_scat << " Mech_" << i;
            }
            out_scat << endl;

            for (int itab = 0; itab <= MTAB; ++itab) {
                double E_real = energy[itab] * eV0;
                if (E_real > 3.0) break;

                double rate_total = sumscatt[itab][band_idx] / time0;
                out_scat << E_real << " " << rate_total;

                for (int i = 0; i < scpre; ++i) {
                    double rate_proc = dose[i][band_idx][itab] / time0;
                    out_scat << " " << rate_proc;
                }
                out_scat << endl;
            }
            out_scat.close();
            cout << "  [Debug] Scattering rates saved to: " << dump_name << endl;
        }
    }
    // =========================================================

    cout << "  Analytic scattering table built. Max Rate (norm) = " << max_gamma << endl;
}
*/


void Band::BuildAnalyticScatteringTable() {
    // =========================================================================
    // Lambda: 从表中获取态密度 (用于 IGZO POP)
    // =========================================================================
    auto get_dos_si_from_table = [&](double E_eV_query) -> double {
        if (E_eV_query < 0) return 0.0;
        const double E_norm = E_eV_query / eV0;
        int itab_q = static_cast<int>(((E_norm - emin) / dtable) + 0.5);
        if (itab_q >= 0 && itab_q <= MTAB) {
            return sumdos[itab_q][PELEC] / (eV0 * std::pow(spr0, 3.0));
        } else {
            return sumdos[MTAB][PELEC] / (eV0 * std::pow(spr0, 3.0));
        }
    };

    // =========================================================================
    // Lambda: 解析计算单谷 Kane DOS (用于 Si 谷间散射 - 复刻旧版逻辑)
    // =========================================================================
    auto calc_kane_dos_si_analytical = [&](double E_eV_val) -> double {
        if (E_eV_val <= 0.0) return 0.0;
        
        // Si 参数 (硬编码，确保与旧版一致)
        const double ml_si = 0.916 * M0_SI;
        const double mt_si = 0.190 * M0_SI;
        const double alpha_si = 0.5; // 1/eV
        const double md_si = std::pow(ml_si * mt_si * mt_si, 1.0/3.0);
        
        double E_J = E_eV_val * Q_SI;
        double alpha_J = alpha_si / Q_SI;
        
        // g(E) = sqrt(2) * md^1.5 / (pi^2 * hbar^3) * sqrt(E(1+aE)) * (1+2aE)
        // 这里的单位是 J^-1 m^-3
        double prefactor = std::sqrt(2.0) * std::pow(md_si, 1.5) / (PI_SI * PI_SI * std::pow(HBAR_SI, 3.0));
        double gamma = E_J * (1.0 + alpha_J * E_J);
        double dgamma = 1.0 + 2.0 * alpha_J * E_J;
        
        return prefactor * std::sqrt(gamma) * dgamma;
    };

    // -------------------------------------------------------------------------
    // 准备积分网格 (用于声学支周期性延展积分)
    // -------------------------------------------------------------------------
    const double q_grid_limit = 6.0e10; 
    const int nq_int = 1000; 
    const double dq_int = q_grid_limit / (nq_int - 1);
    vector<double> q_grid(nq_int);
    for(int i=0; i<nq_int; ++i) q_grid[i] = i * dq_int;

    // -------------------------------------------------------------------------
    // Helper: 周期性对称延展（Zone Folding / triangle-wave folding）
    // 将任意 q 映射回 [0, qmax]，用于查表 omega(q)
    // -------------------------------------------------------------------------
    auto fold_q_to_0_qmax = [&](double q, double qmax) -> double {
        if (!(qmax > 0.0)) return q;
        const double q_period = 2.0 * qmax;
        double q_mod = std::fmod(q, q_period);
        if (q_mod < 0.0) q_mod += q_period;
        return qmax - std::abs(q_mod - qmax);
    };

    // -------------------------------------------------------------------------
    // IGZO 分支 (Acoustic 积分 + POP 查表)
    // -------------------------------------------------------------------------
    if (this->igzofl) {
        cout << "Building Analytic Scattering Table (IGZO: Periodic Extension)..." << endl;

        const double T_lattice = T0;
        const double rho = 6100.0;
        const double E_ac_eV = 5.0; 
        const double D_LA = E_ac_eV * Q_SI; 
        const double D_TA = E_ac_eV * Q_SI;

        const double ml = mell * M0_SI;
        const double mt = melt * M0_SI;
        const double md_SI = std::pow(ml * mt * mt, 1.0/3.0);
        const double alpha_val = 0.0; 

        // ---------------------------------------------------------------------
        // [IGZO-Amorphous] 非晶无序修正因子：低能区指数增强散射
        // delta_E(E) = E_tail * (1 - E/E_max_corr), for E < E_max_corr
        // S_disorder = exp(delta_E / kBT)
        // ---------------------------------------------------------------------
        const double E_tail_eV = 0.12;
        const double E_max_corr_eV = 3.0;
        const double kBT_eV = (KB_SI * T_lattice) / Q_SI;

        // Optical Parameters (NPOP via effective optical deformation potential)
        const double Dopt_eVm = 5e10;
        const double D_LO = Dopt_eVm * Q_SI; // J/m

        scpre = 3; 
        const int band_idx = bandof[PELEC];

        // 清零
        for (int iproc = 0; iproc < scpre; ++iproc) {
            for (int ib = 0; ib < NBE; ++ib) {
                scatte[iproc][ib][ib] = 0.0;
                for (int itab = 0; itab <= MTAB; ++itab) dose[iproc][ib][itab] = 0.0;
            }
            scatte[iproc][band_idx][band_idx] = 1.0;
        }
        for (int ib = 0; ib < NBE; ++ib) for (int jb = 0; jb < NBE; ++jb) for (int itab = 0; itab <= MTAB; ++itab) scattiie[ib][jb][itab] = 0.0;

        // --- 能量循环 ---
        for (int itab = 0; itab <= MTAB; ++itab) {
            const double E_eV = energy[itab] * eV0;
            sumscatt[itab][band_idx] = 0.0;
            dose[0][band_idx][itab] = 0.0;
            dose[1][band_idx][itab] = 0.0;
            dose[2][band_idx][itab] = 0.0;

            // [IGZO-Amorphous] 计算能量相关无序增强因子
            double S_disorder = 1.0;
            if (kBT_eV > 0.0 && E_eV < E_max_corr_eV) {
                const double delta_E = E_tail_eV * (1.0 - (E_eV / E_max_corr_eV));
                const double expo = delta_E / kBT_eV;
                S_disorder = std::exp(std::min(expo, 700.0)); // 防止溢出；常规参数下 expo ~ O(1-10)
            }

            // 1. Acoustic Scattering (使用周期性延展积分)
            double Rate_AC = 0.0;
            double term_k = E_eV * (1.0 + alpha_val * E_eV);
            if (term_k < 0) term_k = 0.0;
            double ks = std::sqrt(2.0 * md_SI * term_k * Q_SI) / HBAR_SI;

            if (ks > 1e-10) {
                double sum_integ_LA = 0.0;
                double sum_integ_TA = 0.0;
                double sum_LO_abs = 0.0;
                double sum_LO_ems = 0.0;

                for (int iq = 0; iq < nq_int; ++iq) {
                    double q = q_grid[iq];
                    if (q < 1e-12) continue;

                    const double q_mapped = fold_q_to_0_qmax(q, phonon.qmax);
                    const double q_strength = (q_mapped > 1e-12) ? q_mapped : 1e-12;

                    // 查表必须用映射后的 q；强度也必须使用 q_mapped，避免高能区奇点
                    double w_LA = GetPhononOmega(PH_LA, q_strength);
                    double w_TA = GetPhononOmega(PH_TA, q_strength);
                    double w_LO = GetPhononOmega(PH_LO, q_strength);
                    if (!(w_LA > 0.0) && phonon.nq_tab > 1 && phonon.dq > 0.0 && phonon.omega_table[PH_LA].size() > 1) {
                        w_LA = (phonon.omega_table[PH_LA][1] / phonon.dq) * q_strength;
                    }
                    if (!(w_TA > 0.0) && phonon.nq_tab > 1 && phonon.dq > 0.0 && phonon.omega_table[PH_TA].size() > 1) {
                        w_TA = (phonon.omega_table[PH_TA][1] / phonon.dq) * q_strength;
                    }
                    if (!(w_LA > 0.0) || !(w_TA > 0.0) || !(w_LO > 0.0)) continue;

                    double hw_LA = HBAR_SI * w_LA / Q_SI; 
                    double hw_TA = HBAR_SI * w_TA / Q_SI; 
                    double hw_LO = HBAR_SI * w_LO / Q_SI;

                    double N_LA = 1.0 / (std::exp(w_LA * HBAR_SI / (KB_SI * T_lattice)) - 1.0);
                    double N_TA = 1.0 / (std::exp(w_TA * HBAR_SI / (KB_SI * T_lattice)) - 1.0);
                    double N_LO = 1.0 / (std::exp(w_LO * HBAR_SI / (KB_SI * T_lattice)) - 1.0);

                    const double base_term = q_mapped * q_mapped * q_mapped;
                    const double base_opt = q_strength; // optical strength term uses q_mapped (mapped into 1st BZ)
                    
                    if (CheckAllowedQ(E_eV, ks, q, hw_LA, -1, md_SI, alpha_val)) { 
                        if (E_eV > hw_LA) sum_integ_LA += (1.0/w_LA) * (N_LA + 1.0) * base_term;
                    }
                    if (CheckAllowedQ(E_eV, ks, q, hw_LA, 1, md_SI, alpha_val)) { 
                        sum_integ_LA += (1.0/w_LA) * N_LA * base_term;
                    }
                    if (CheckAllowedQ(E_eV, ks, q, hw_TA, -1, md_SI, alpha_val)) {
                        if (E_eV > hw_TA) sum_integ_TA += (1.0/w_TA) * (N_TA + 1.0) * base_term;
                    }
                    if (CheckAllowedQ(E_eV, ks, q, hw_TA, 1, md_SI, alpha_val)) {
                        sum_integ_TA += (1.0/w_TA) * N_TA * base_term;
                    }

                    // LO optical: absorption / emission (Mech 1/2)
                    if (CheckAllowedQ(E_eV, ks, q, hw_LO, 1, md_SI, alpha_val)) {
                        sum_LO_abs += (1.0 / w_LO) * N_LO * base_opt;
                    }
                    if (E_eV > hw_LO && CheckAllowedQ(E_eV, ks, q, hw_LO, -1, md_SI, alpha_val)) {
                        sum_LO_ems += (1.0 / w_LO) * (N_LO + 1.0) * base_opt;
                    }
                }
                
                double pre = md_SI / (4.0 * PI_SI * rho * HBAR_SI * HBAR_SI * ks);
                Rate_AC = pre * (D_LA*D_LA * sum_integ_LA + D_TA*D_TA * sum_integ_TA) * dq_int;

                const double Rate_LO_abs = pre * (D_LO * D_LO) * sum_LO_abs * dq_int;
                const double Rate_LO_ems = pre * (D_LO * D_LO) * sum_LO_ems * dq_int;

                dose[0][band_idx][itab] = (Rate_AC * S_disorder) * time0;
                dose[1][band_idx][itab] = (Rate_LO_abs * S_disorder) * time0;
                dose[2][band_idx][itab] = (Rate_LO_ems * S_disorder) * time0;
                sumscatt[itab][band_idx] = dose[0][band_idx][itab] + dose[1][band_idx][itab] + dose[2][band_idx][itab];
            }
        }

    } else {
        // ---------------------------------------------------------------------
        // Silicon 分支
        // ---------------------------------------------------------------------
        cout << "Building Analytic Scattering Table (Si: Periodic Acoustic + Analytical Intervalley)..." << endl;

        double T_lattice = T0;
        double rho = 2330.0;
        double D_LA_J = 6.39 * Q_SI;
        double D_TA_J = 3.01 * Q_SI;
        
        double ml = 0.916 * M0_SI;
        double mt = 0.190 * M0_SI;
        double md_SI = std::pow(ml * mt * mt, 1.0/3.0);
        double alpha_val = 0.5;

        double a0 = phonon.a0;
        double Rs = (a0 > 0) ? a0 * std::pow(3.0/(16.0*PI_SI), 1.0/3.0) : 0.0;

        struct IvParam { double E_meV; double D_1e8; int Z; };
        IvParam iv_params[6] = {
            {10.0, 0.3, 1}, {19.0, 1.5, 1}, {62.0, 6.0, 1},
            {19.0, 0.5, 4}, {51.0, 3.5, 4}, {57.0, 1.5, 4}
        };

        scpre = 14; 
        int band_idx = bandof[PELEC];

        // 清零
        for (int iproc = 0; iproc < scpre; ++iproc) {
            for (int ib = 0; ib < NBE; ++ib) {
                scatte[iproc][ib][ib] = 0.0;
                for (int itab = 0; itab <= MTAB; ++itab) dose[iproc][ib][itab] = 0.0;
            }
            scatte[iproc][band_idx][band_idx] = 1.0;
        }
        for (int ib = 0; ib < NBE; ++ib) for (int jb = 0; jb < NBE; ++jb) for (int itab = 0; itab <= MTAB; ++itab) scattiie[ib][jb][itab] = 0.0;

        for (int itab = 0; itab <= MTAB; ++itab) {
            double E_norm = energy[itab];
            double E_eV = E_norm * eV0;
            sumscatt[itab][band_idx] = 0.0;
            
            // 1. Acoustic: Periodic Extension
            double Rate_AC = 0.0;
            double term_k = E_eV * (1.0 + alpha_val * E_eV);
            if (term_k < 0) term_k = 0.0;
            double ks = std::sqrt(2.0 * md_SI * term_k * Q_SI) / HBAR_SI;

            if (ks > 1e-10) {
                double sum_LA = 0.0, sum_TA = 0.0;
                for (int iq = 0; iq < nq_int; ++iq) {
                    double q = q_grid[iq];
                    if (q < 1e-12) continue;

                    const double q_mapped = fold_q_to_0_qmax(q, phonon.qmax);
                    const double q_strength = (q_mapped > 1e-12) ? q_mapped : 1e-12;

                    double w_LA = GetPhononOmega(PH_LA, q_strength);
                    double w_TA = GetPhononOmega(PH_TA, q_strength);
                    if (!(w_LA > 0.0) && phonon.nq_tab > 1 && phonon.dq > 0.0 && phonon.omega_table[PH_LA].size() > 1) {
                        w_LA = (phonon.omega_table[PH_LA][1] / phonon.dq) * q_strength;
                    }
                    if (!(w_TA > 0.0) && phonon.nq_tab > 1 && phonon.dq > 0.0 && phonon.omega_table[PH_TA].size() > 1) {
                        w_TA = (phonon.omega_table[PH_TA][1] / phonon.dq) * q_strength;
                    }
                    if (!(w_LA > 0.0) || !(w_TA > 0.0)) continue;

                    double hw_LA = HBAR_SI * w_LA / Q_SI;
                    double hw_TA = HBAR_SI * w_TA / Q_SI;

                    double N_LA = 1.0 / (std::exp(w_LA * HBAR_SI / (KB_SI * T_lattice)) - 1.0);
                    double N_TA = 1.0 / (std::exp(w_TA * HBAR_SI / (KB_SI * T_lattice)) - 1.0);

                    const double Iq = (Rs > 0) ? GetOverlapFactor(q_strength, Rs) : 1.0;
                    const double base_term = (Iq * Iq) * (q_mapped * q_mapped * q_mapped);

                    // LA
                    if (CheckAllowedQ(E_eV, ks, q, hw_LA, -1, md_SI, alpha_val)) {
                        if (E_eV > hw_LA) sum_LA += (1.0/w_LA) * (N_LA + 1.0) * base_term;
                    }
                    if (CheckAllowedQ(E_eV, ks, q, hw_LA, 1, md_SI, alpha_val)) {
                        sum_LA += (1.0/w_LA) * N_LA * base_term;
                    }
                    // TA
                    if (CheckAllowedQ(E_eV, ks, q, hw_TA, -1, md_SI, alpha_val)) {
                        if (E_eV > hw_TA) sum_TA += (1.0/w_TA) * (N_TA + 1.0) * base_term;
                    }
                    if (CheckAllowedQ(E_eV, ks, q, hw_TA, 1, md_SI, alpha_val)) {
                        sum_TA += (1.0/w_TA) * N_TA * base_term;
                    }
                }
                double pre = md_SI / (4.0 * PI_SI * rho * HBAR_SI * HBAR_SI * ks);
                Rate_AC = pre * (D_LA_J*D_LA_J * sum_LA + D_TA_J*D_TA_J * sum_TA) * dq_int;
            }

            dose[0][band_idx][itab] = Rate_AC * time0;
            sumscatt[itab][band_idx] += dose[0][band_idx][itab];

            // 2. Intervalley Scattering (Fix: 使用 calc_kane_dos_si_analytical 替代查表)
            for (int i = 0; i < 6; ++i) {
                double hw_eV = iv_params[i].E_meV * 1e-3;
                double w0 = hw_eV * Q_SI / HBAR_SI;
                double Nq = 1.0 / (std::exp(hw_eV * Q_SI / (KB_SI * T_lattice)) - 1.0);
                double D_Jm = iv_params[i].D_1e8 * 1e10 * Q_SI;
                double C = (PI_SI * D_Jm * D_Jm * iv_params[i].Z) / (2.0 * rho * w0);

                int idx_abs = 1 + 2 * i;
                // [FIX] 计算单谷 DOS
                double g_ab = calc_kane_dos_si_analytical(E_eV + hw_eV); 
                double Rate_Abs_SI = C * Nq * g_ab;
                dose[idx_abs][band_idx][itab] = Rate_Abs_SI * time0;
                sumscatt[itab][band_idx] += dose[idx_abs][band_idx][itab];

                int idx_em = 2 + 2 * i;
                double Rate_Em_SI = 0.0;
                if (E_eV > hw_eV) {
                    // [FIX] 计算单谷 DOS
                    double g_em = calc_kane_dos_si_analytical(E_eV - hw_eV);
                    Rate_Em_SI = C * (Nq + 1.0) * g_em;
                }
                dose[idx_em][band_idx][itab] = Rate_Em_SI * time0;
                sumscatt[itab][band_idx] += dose[idx_em][band_idx][itab];
            }
            
            scattiie[band_idx][band_idx][itab] = 0.0; 
        }
    }

    // 后处理
    double max_gamma = 0.0;
    for (int itab = 0; itab <= MTAB; ++itab) {
        if (sumscatt[itab][bandof[PELEC]] > max_gamma) {
            max_gamma = sumscatt[itab][bandof[PELEC]];
        }
    }
    if (nt <= 0) nt = 1;
    int fill_nt = (nt > MNTet) ? MNTet : nt;
    for(int it = 0; it < fill_nt; it++) gamtet[it] = max_gamma;
    gamma[PELEC] = max_gamma;

    // [DEBUG] 导出
    if (mpi_rank == 0) {
        const string base_dir = pathname.empty() ? string("input") : pathname;
        string dump_name = this->igzofl ? (base_dir + "/scattering_rates_IGZO.txt")
                                        : (base_dir + "/scattering_rates_Si.txt");
        ofstream out_scat(dump_name.c_str());
        if (out_scat.is_open()) {
            out_scat << "Energy(eV) Total(1/s)";
            for(int i=0; i<scpre; i++) out_scat << " Mech_" << i;
            out_scat << endl;
            for (int itab = 0; itab <= MTAB; ++itab) {
                double E_real = energy[itab] * eV0;
                if (E_real > 3.0) break;
                out_scat << E_real << " " << sumscatt[itab][bandof[PELEC]] / time0;
                for (int i = 0; i < scpre; ++i) out_scat << " " << dose[i][bandof[PELEC]][itab] / time0;
                out_scat << endl;
            }
            out_scat.close();
        }
    }
    cout << "  Analytic scattering table built. Max Rate (norm) = " << max_gamma << endl;
}



int Band::GetAxisIndex(double k_norm) {
    int idx = static_cast<int>((k_norm + 2.0) / 0.01 + 0.5);
    if (idx < 0) idx = 0;
    return idx;
}

void Band::InitAxisLookupTable() {
    cout << "Initializing Axis Lookup Table (O(1) Map)..." << endl;

    std::vector<double> ticks = GenerateNonUniformTicks();
    num_ticks_axis = static_cast<int>(ticks.size());

    // k 刻度（代码单位）
    double a_lattice = sia0;
    double conversion = (PI / a_lattice);
    k_ticks_code.resize(num_ticks_axis);
    for (int i = 0; i < num_ticks_axis; ++i) {
        k_ticks_code[i] = ticks[i] * conversion;
    }

    // 计算 K 空间边界（Voronoi 中点）
    k_boundaries.clear();
    if (num_ticks_axis > 1) {
        k_boundaries.resize(num_ticks_axis - 1);
        for (int i = 0; i < num_ticks_axis - 1; ++i) {
            k_boundaries[i] = 0.5 * (k_ticks_code[i] + k_ticks_code[i + 1]);
        }
    }
    cout << "  [Info] K-Space Boundaries computed. Count: " << k_boundaries.size() << endl;

    k_map_min = -2.15;
    k_map_max = 2.15;
    double resolution = 0.001;
    k_map_scale = 1.0 / resolution;

    int map_size = static_cast<int>((k_map_max - k_map_min) * k_map_scale) + 1;
    k_axis_map.resize(map_size);

    for (int i = 0; i < map_size; ++i) {
        double k_val = k_map_min + i * resolution;

        std::vector<double>::iterator it = std::lower_bound(ticks.begin(), ticks.end(), k_val);

        int idx = 0;
        if (it == ticks.begin()) {
            idx = 0;
        } else if (it == ticks.end()) {
            idx = num_ticks_axis - 1;
        } else {
            double val_upper = *it;
            double val_lower = *(it - 1);

            if ((val_upper - k_val) < (k_val - val_lower)) {
                idx = static_cast<int>(it - ticks.begin());
            } else {
                idx = static_cast<int>(it - ticks.begin()) - 1;
            }
        }
        k_axis_map[i] = idx;
    }

    cout << "  Axis Map built. Size: " << map_size << ", Resolution: " << resolution << endl;
    cout << "  Grid dimensions: " << num_ticks_axis << "^3" << endl;
}

int Band::GetAxisIndex_O1(double k_val) {
    int map_idx = static_cast<int>((k_val - k_map_min) * k_map_scale);
    if (map_idx < 0) return 0;
    if (map_idx >= static_cast<int>(k_axis_map.size())) return num_ticks_axis - 1;
    return k_axis_map[map_idx];
}

double Band::GetAnalyticGridTime(Particle* p, double Fx, double Fy, double Fz) {
    double min_dt = 1.0e99;
    last_k_col_dir = -1;

    // 使用粒子已有索引，避免重复查找
    auto update_dt = [&](double k_val, int idx, double F, int axis_id) {
        if (std::fabs(F) < 1.0e-20) return;

        double wall;
        int dir_offset = (F > 0) ? 0 : 1; // 0:+, 1:-

        if (F > 0) {
            if (idx >= num_ticks_axis - 1) wall = 1.0e99;
            else wall = k_boundaries[idx];
        } else {
            if (idx <= 0) wall = -1.0e99;
            else wall = k_boundaries[idx - 1];
        }

        double dt = (wall - k_val) / F;
        if (dt > 1.0e-14 && dt < min_dt) {
            min_dt = dt;
            last_k_col_dir = axis_id * 2 + dir_offset;
        }
    };

    update_dt(p->kx, p->kx_idx, Fx, 0); // X
    update_dt(p->ky, p->ky_idx, Fy, 1); // Y
    update_dt(p->kz, p->kz_idx, Fz, 2); // Z

    if (min_dt < 1.0e-8) min_dt = 1.0e-8;
    return min_dt;
}

void Band::GetAnalyticV_By_Index(Particle* p) {
    if (p->kx_idx < 0) p->kx_idx = 0;
    if (p->kx_idx >= num_ticks_axis) p->kx_idx = num_ticks_axis - 1;
    if (p->ky_idx < 0) p->ky_idx = 0;
    if (p->ky_idx >= num_ticks_axis) p->ky_idx = num_ticks_axis - 1;
    if (p->kz_idx < 0) p->kz_idx = 0;
    if (p->kz_idx >= num_ticks_axis) p->kz_idx = num_ticks_axis - 1;

    int N = num_ticks_axis;
    int flat_idx = p->kx_idx * (N * N) + p->ky_idx * N + p->kz_idx;
    if (flat_idx >= 0 && flat_idx < static_cast<int>(analytic_k_grid.size())) {
        const AnalyticKPoint& pt = analytic_k_grid[flat_idx];
        analytic_vx = pt.vx;
        analytic_vy = pt.vy;
        analytic_vz = pt.vz;
        p->energy = pt.energy;
    } else {
        analytic_vx = analytic_vy = analytic_vz = 0.0;
    }
}

double Band::GetAnalyticImpurityRate(double E, double DA, double Rho, double eps_si, double frickel) {
    // 阈值判断，与 GetImpScRate 中保持一致
    double rvscrt = 1.0/scrt0;
    if ((E * eV0 >= 0.120) || (DA * conc0 < 1e22)) {
    //if (false) {
        return rvscrt;
    }

    double eebeta = frickel * 4.0 * PI * Rho / (2.0 * meld * eps_si);
    if (std::abs(eebeta) < 1.0e-20) eebeta = 1.0e-20;

    double denom_core = PI * std::pow(frickel * Rho, 2.0);
    if (denom_core < 1.0e-30) return 0.0;

    double term_in_bracket = 4.0 * E * melt / (eebeta * meld);
    double denom_full = denom_core * (1.0 + term_in_bracket);
    if (denom_full < 1.0e-30) return 0.0;

    rvscrt = DA * meld * std::sqrt(2.0 * meld * E) / denom_full;
    rvscrt = Max(rvscrt,1.0/scrt0);
    return rvscrt;
}

void Band::AnalyticImpurityScatter(Particle* p, double DA, double Rho, double eps_si, double frickel, double ImpScGamma_Max) {
    analytic_self_scatter = false;

    double gamimp = GetAnalyticImpurityRate(p->energy, DA, Rho, eps_si, frickel);
    if (gamimp < Random() * ImpScGamma_Max) {
        analytic_self_scatter = true;
        return;
    }

    double k_conv_real = 1.0 / spr0;
    double kx_real = p->kx * k_conv_real;
    double ky_real = p->ky * k_conv_real;
    double kz_real = p->kz * k_conv_real;

    double a_lattice = sia0;
    double K0_code = 0.0;
    if (!this->igzofl) {
        double K_valley_mag = 0.85 * (2.0 * PI / a_lattice);
        K0_code = K_valley_mag;
    }

    double kl = 0.0, kt1 = 0.0, kt2 = 0.0;
    int axis = 0;
    if (std::fabs(p->kx) >= std::fabs(p->ky) && std::fabs(p->kx) >= std::fabs(p->kz)) {
        axis = 0;
        kl = (p->kx > 0) ? (p->kx - K0_code) : (p->kx + K0_code);
        kt1 = p->ky; kt2 = p->kz;
    } else if (std::fabs(p->ky) >= std::fabs(p->kx) && std::fabs(p->ky) >= std::fabs(p->kz)) {
        axis = 1;
        kl = (p->ky > 0) ? (p->ky - K0_code) : (p->ky + K0_code);
        kt1 = p->kx; kt2 = p->kz;
    } else {
        axis = 2;
        kl = (p->kz > 0) ? (p->kz - K0_code) : (p->kz + K0_code);
        kt1 = p->kx; kt2 = p->ky;
    }

    double scale_l = std::sqrt(meld / mell);
    double scale_t = std::sqrt(meld / melt);
    double x_bh = kl * scale_l;
    double y_bh = kt1 * scale_t;
    double z_bh = kt2 * scale_t;

    double betaq = frickel * 4.0 * PI * Rho / eps_si;
    double eebeta = betaq / (2.0 * meld);
    double alfa = eebeta * meld / (2.0 * p->energy * melt);
    double pr = Random() / (1.0 + 0.5 * alfa);
    double costr = 1.0 - alfa * pr / (1.0 - pr);
    if (costr > 1.0) costr = 1.0;
    if (costr < -1.0) costr = -1.0;
    double sintr = std::sqrt(std::max(0.0, 1.0 - costr * costr));
    double phi = 2.0 * PI * Random();
    double cospr = std::cos(phi);
    double sinpr = std::sin(phi);

    double k_mag = std::sqrt(x_bh * x_bh + y_bh * y_bh + z_bh * z_bh);
    if (k_mag < 1e-20) k_mag = 1e-20;
    double cos_theta_old = z_bh / k_mag;
    double sin_theta_old = std::sqrt(std::max(0.0, 1.0 - cos_theta_old * cos_theta_old));
    double cos_phi_old = 1.0, sin_phi_old = 0.0;
    if (sin_theta_old > 1e-9) {
        cos_phi_old = x_bh / (k_mag * sin_theta_old);
        sin_phi_old = y_bh / (k_mag * sin_theta_old);
    }

    double x_new_prime = k_mag * sintr * cospr;
    double y_new_prime = k_mag * sintr * sinpr;
    double z_new_prime = k_mag * costr;

    double x_bh_new = cos_phi_old * cos_theta_old * x_new_prime - sin_phi_old * y_new_prime + cos_phi_old * sin_theta_old * z_new_prime;
    double y_bh_new = sin_phi_old * cos_theta_old * x_new_prime + cos_phi_old * y_new_prime + sin_phi_old * sin_theta_old * z_new_prime;
    double z_bh_new = -sin_theta_old * x_new_prime + cos_theta_old * z_new_prime;

    double kl_new = x_bh_new * std::sqrt(mell / meld);
    double kt1_new = y_bh_new * std::sqrt(melt / meld);
    double kt2_new = z_bh_new * std::sqrt(melt / meld);

    if (axis == 0) {
        p->kx = (p->kx > 0) ? (K0_code + kl_new) : (-K0_code + kl_new);
        p->ky = kt1_new;
        p->kz = kt2_new;
    } else if (axis == 1) {
        p->ky = (p->ky > 0) ? (K0_code + kl_new) : (-K0_code + kl_new);
        p->kx = kt1_new;
        p->kz = kt2_new;
    } else {
        p->kz = (p->kz > 0) ? (K0_code + kl_new) : (-K0_code + kl_new);
        p->kx = kt1_new;
        p->ky = kt2_new;
    }

    // 更新索引并查表速度
    double to_pi = 1.0 / (PI / a_lattice);
    p->kx_idx = GetAxisIndex_O1(p->kx * to_pi);
    p->ky_idx = GetAxisIndex_O1(p->ky * to_pi);
    p->kz_idx = GetAxisIndex_O1(p->kz * to_pi);

    GetAnalyticV_FromTable(p);
}

// -----------------------------------------------------------------------------
// 解析能带：声子散射与能谷操作模块
// -----------------------------------------------------------------------------

void Band::InitValleyConfiguration() {
    if (this->igzofl) {
        if (mpi_rank == 0) cout << "  [Band] Configuring Single Valley at Gamma (IGZO)..." << endl;
        valley_k0_norm = 0.0;
        for (int i = 0; i < 6; ++i) {
            valley_centers[i][0] = 0.0;
            valley_centers[i][1] = 0.0;
            valley_centers[i][2] = 0.0;
            valley_axis[i] = 0;
        }
        return;
    }

    double a_lattice = sia0;
    double K_code = 0.85 * (2.0 * PI / a_lattice);
    valley_k0_norm = K_code;

    cout << "Initializing Valley Config. K0_norm (Code Unit) = " << valley_k0_norm << endl;

    // +X / -X
    valley_centers[0][0] = valley_k0_norm;  valley_centers[0][1] = 0.0;               valley_centers[0][2] = 0.0; valley_axis[0] = 0;
    valley_centers[1][0] = -valley_k0_norm; valley_centers[1][1] = 0.0;               valley_centers[1][2] = 0.0; valley_axis[1] = 0;
    // +Y / -Y
    valley_centers[2][0] = 0.0;             valley_centers[2][1] = valley_k0_norm;    valley_centers[2][2] = 0.0; valley_axis[2] = 1;
    valley_centers[3][0] = 0.0;             valley_centers[3][1] = -valley_k0_norm;   valley_centers[3][2] = 0.0; valley_axis[3] = 1;
    // +Z / -Z
    valley_centers[4][0] = 0.0;             valley_centers[4][1] = 0.0;               valley_centers[4][2] = valley_k0_norm; valley_axis[4] = 2;
    valley_centers[5][0] = 0.0;             valley_centers[5][1] = 0.0;               valley_centers[5][2] = -valley_k0_norm; valley_axis[5] = 2;
}

int Band::GetValleyID(double kx, double ky, double kz) {
    int best_id = 0;
    double min_d2 = 1.0e99;
    for (int i = 0; i < 6; ++i) {
        double dx = kx - valley_centers[i][0];
        double dy = ky - valley_centers[i][1];
        double dz = kz - valley_centers[i][2];
        double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < min_d2) {
            min_d2 = d2;
            best_id = i;
        }
    }
    return best_id;
}

void Band::AnalyticPhononScatter(Particle* p) {
    int itab = static_cast<int>((p->energy - emin) * dlist);
    if (itab < 0) itab = 0;
    if (itab >= MWLE_ana) itab = MWLE_ana - 1;

    int band_idx = bandof[PELEC];
    double total_rate = sumscatt[itab][band_idx];

    analytic_self_scatter = false;
    if (Random() * gamma[PELEC] > total_rate) {
        analytic_self_scatter = true;
        return;
    }

    double r_mech = Random() * total_rate;
    double acc = 0.0;
    int iscat = -1;
    for (int i = 0; i < scpre; ++i) {
        acc += dose[i][band_idx][itab];
        if (acc >= r_mech) { iscat = i; break; }
    }
    if (iscat < 0) iscat = scpre - 1;

    double delta_E_eV = 0.0;
    int valley_rule = 0; // 0 intra, 1 g-opposite, 2 f-perp
    static const double E_ph_meV[] = {
        0, 10.0,10.0, 19.0,19.0, 62.0,62.0, 19.0,19.0, 51.0,51.0, 57.0,57.0
    };

    if (this->igzofl) {
        // IGZO: 0 acoustic (elastic), 1 optical absorption, 2 optical emission (intravalley only)
        double omega_LO = 0.0;
        if (phonon.nq_tab > 0) omega_LO = phonon.omega_table[PH_LO].front();
        if (!(omega_LO > 0.0) && phonon.nq_tab > 1) omega_LO = phonon.omega_table[PH_LO][1];
        if (!(omega_LO > 0.0) && phonon.nq_tab > 2) omega_LO = phonon.omega_table[PH_LO][phonon.nq_tab / 2];
        double hw = (omega_LO > 0.0) ? (HBAR_SI * omega_LO / Q_SI) : 0.06;
        if (hw < 0.02) hw = 0.06;

        if (iscat == 0) {
            delta_E_eV = 0.0;
        } else if (iscat == 1) {
            delta_E_eV = hw;
        } else if (iscat == 2) {
            delta_E_eV = -hw;
            if (p->energy * eV0 < hw) {
                analytic_self_scatter = true;
                return;
            }
        } else {
            analytic_self_scatter = true;
            return;
        }
        valley_rule = 0;
    } else {
        if (iscat == 0) {
            delta_E_eV = 0.0;
            valley_rule = 0;
        } else if (iscat >= 1 && iscat <= 6) {
            double hw = E_ph_meV[iscat] * 1e-3;
            delta_E_eV = (iscat % 2 != 0) ? hw : -hw;
            valley_rule = 1;
        } else if (iscat >= 7 && iscat <= 12) {
            double hw = E_ph_meV[iscat] * 1e-3;
            delta_E_eV = (iscat % 2 != 0) ? hw : -hw;
            valley_rule = 2;
        } else {
            analytic_self_scatter = true;
            return;
        }
    }

    int current_valley = GetValleyID(p->kx, p->ky, p->kz);
    int target_valley = current_valley;
    if (this->igzofl) {
        current_valley = 0;
        target_valley = 0;
    }
    if (valley_rule == 1) {
        target_valley = current_valley ^ 1;
    } else if (valley_rule == 2) {
        int axis = valley_axis[current_valley];
        int r = static_cast<int>(Random() * 4.0);
        static const int axis_map[3][2] = { {1,2}, {0,2}, {0,1} };
        int next_axis = axis_map[axis][r / 2];
        int next_dir = r % 2;
        target_valley = next_axis * 2 + next_dir;
    }

    double E_final_norm = p->energy + delta_E_eV / eV0;
    if (E_final_norm < 1e-6) E_final_norm = 1e-6;

    Particle temp_p;
    SelectAnalyticKState(&temp_p, E_final_norm);
    int sample_valley = GetValleyID(temp_p.kx, temp_p.ky, temp_p.kz);

    double dkx = temp_p.kx - valley_centers[sample_valley][0];
    double dky = temp_p.ky - valley_centers[sample_valley][1];
    double dkz = temp_p.kz - valley_centers[sample_valley][2];

    double ql, qt1, qt2;
    int s_axis = valley_axis[sample_valley];
    if (s_axis == 0)      { ql = dkx; qt1 = dky; qt2 = dkz; }
    else if (s_axis == 1) { ql = dky; qt1 = dkx; qt2 = dkz; }
    else                  { ql = dkz; qt1 = dkx; qt2 = dky; }

    int t_axis = valley_axis[target_valley];
    double nkx=0.0, nky=0.0, nkz=0.0;
    if (t_axis == 0) {
        nkx = valley_centers[target_valley][0] + ql;
        nky = qt1;
        nkz = qt2;
    } else if (t_axis == 1) {
        nky = valley_centers[target_valley][1] + ql;
        nkx = qt1;
        nkz = qt2;
    } else {
        nkz = valley_centers[target_valley][2] + ql;
        nkx = qt1;
        nky = qt2;
    }

    p->kx = nkx;
    p->ky = nky;
    p->kz = nkz;

    // 更新索引并查表
    double a_lattice = sia0;
    double to_pi = 1.0 / (PI / a_lattice);
    p->kx_idx = GetAxisIndex_O1(p->kx * to_pi);
    p->ky_idx = GetAxisIndex_O1(p->ky * to_pi);
    p->kz_idx = GetAxisIndex_O1(p->kz * to_pi);

    GetAnalyticV_FromTable(p);
}
// 通过索引直接查表，避免浮点抖动
void Band::GetAnalyticV_FromTable(Particle* p) {
    int N = num_ticks_axis;

    // 越界钳位
    if (p->kx_idx < 0) p->kx_idx = 0; else if (p->kx_idx >= N) p->kx_idx = N - 1;
    if (p->ky_idx < 0) p->ky_idx = 0; else if (p->ky_idx >= N) p->ky_idx = N - 1;
    if (p->kz_idx < 0) p->kz_idx = 0; else if (p->kz_idx >= N) p->kz_idx = N - 1;

    int flat_idx = p->kx_idx * (N * N) + p->ky_idx * N + p->kz_idx;

    if (flat_idx >= 0 && flat_idx < static_cast<int>(analytic_k_grid.size())) {
        const AnalyticKPoint& pt = analytic_k_grid[flat_idx];
        analytic_vx = pt.vx;
        analytic_vy = pt.vy;
        analytic_vz = pt.vz;
        p->energy = pt.energy;
    } else {
        analytic_vx = analytic_vy = analytic_vz = 0.0;
        cout << "Error: Index out of bounds in GetAnalyticV_FromTable: " << flat_idx << endl;
    }
}

// 根据粒子索引查表，返回速度与能量（不修改输入粒子）
void Band::GetAnalyticStateByIndex(Particle* p, double &vx_out, double &vy_out, double &vz_out, double &E_out) {
    vx_out = vy_out = vz_out = 0.0;
    E_out = 0.0;
    if (!use_analytic_band || p == nullptr) return;

    Particle tmp = *p; // 不影响原粒子
    GetAnalyticV_FromTable(&tmp);
    vx_out = analytic_vx * velo0;
    vy_out = analytic_vy * velo0;
    vz_out = analytic_vz * velo0;
    E_out = tmp.energy * pot0;
}

void Band::BuildAnalyticInjectionTable() {
    cout << "Building Analytic Injection Table (Output to file)..." << endl;

    double Ef_min_eV = -0.4;
    double Ef_max_eV = 0.4;
    int n_Ef_steps = 10000;

    double T_lattice = T0; 
    double kb_T_eV = T_lattice * 8.617333262e-5; 

    string filename = pathname + "/analytic_density_cross_number.txt";
    ofstream out(filename.c_str());
    if (!out) {
        cerr << "Error: Cannot open output file " << filename << endl;
        return;
    }
    out << n_Ef_steps << endl;
    out << Ef_min_eV << " " << Ef_max_eV << endl;

    double dEf_eV = (Ef_max_eV - Ef_min_eV) / (n_Ef_steps - 1);

    int itab_max_integ = static_cast<int>(1.5 / (dtable * eV0));
    if (itab_max_integ > MTAB) itab_max_integ = MTAB;
    double dE_eV = dtable * eV0; 

    vector<double> dos_SI_cache(itab_max_integ);
    vector<double> v_SI_cache(itab_max_integ);

    double hbar = 1.0545718e-34;
    double ml = 0.916 * 9.1093837e-31;
    double mt = 0.190 * 9.1093837e-31;
    double md = std::pow(ml * mt * mt, 1.0/3.0);
    double alpha = 0.5;

    for (int itab = 0; itab < itab_max_integ; ++itab) {
        double E_norm = energy[itab];
        double E_eV = E_norm * eV0;

        double g_norm = sumdos[itab][PELEC];
        dos_SI_cache[itab] = g_norm / (eV0 * std::pow(spr0, 3)); 

        double k_SI = GetKaneK_SI(E_eV);
        double denom = md * (1.0 + 2.0 * alpha * E_eV);
        if (denom > 0) v_SI_cache[itab] = (hbar * k_SI) / denom;
        else v_SI_cache[itab] = 0.0;
    }

    for (int i = 0; i < n_Ef_steps; ++i) {
        double current_Ef_eV = Ef_min_eV + i * dEf_eV;

        double n_accum_SI = 0.0; 
        double J_accum_SI = 0.0; 

        for (int itab = 0; itab < itab_max_integ; ++itab) {
            double E_norm = energy[itab];
            double E_eV = E_norm * eV0;

            double exponent = (E_eV - current_Ef_eV) / kb_T_eV;
            double f_dist = 0.0;
            if (exponent > 50.0) f_dist = 0.0;
            else if (exponent < -50.0) f_dist = 1.0;
            else f_dist = 1.0 / (1.0 + std::exp(exponent));

            double g = dos_SI_cache[itab];
            double v = v_SI_cache[itab];

            n_accum_SI += g * f_dist * dE_eV;
            J_accum_SI += g * f_dist * (v * 0.25) * dE_eV; // 半球积分 1/4
        }

        out << current_Ef_eV << " " << n_accum_SI << " " << J_accum_SI << endl;
    }

    out.close();
    cout << "  Analytic injection table generated: " << filename << endl;

    // 生成后立即读入，保持内存与文件一致
    ReadAnalyticInjectionTable();
}

void Band::ReadAnalyticInjectionTable() {
    string filename = pathname + "/analytic_density_cross_number.txt";
    cout << "Reading Analytic Injection Table from: " << filename << endl;

    ifstream in(filename.c_str());
    if (!in) {
        cerr << "Error: Cannot open injection table file." << endl;
        exit(1);
    }

    int n_steps;
    double ef_min_v, ef_max_v;
    in >> n_steps;
    in >> ef_min_v >> ef_max_v;

    NEf = n_steps;
    Efmin = ef_min_v / eV0; // 归一化
    Efmax = ef_max_v / eV0;
    deltaEf = (Efmax - Efmin) / (NEf - 1);

    Ef_value.resize(NEf);
    Ef_density.resize(NEf);
    Ef_cross_number.resize(NEf);

    double ef_read, n_read, j_read;
    for (int i = 0; i < NEf; ++i) {
        in >> ef_read >> n_read >> j_read;

        Ef_value[i] = ef_read / eV0;
        Ef_density[i] = n_read / conc0;                      // m^-3 -> norm
        Ef_cross_number[i] = j_read * (spr0 * spr0 * time0); // m^-2 s^-1 -> norm
    }

    in.close();
    cout << "  Injection table loaded. Range: [" << ef_min_v << ", " << ef_max_v << "] eV." << endl;

    int mid = NEf / 2;
    cout << "  [Check] Ef= " << Ef_value[mid]*eV0 << " eV, "
         << "n= " << Ef_density[mid]*conc0*1e-6 << " cm^-3" << endl;
}

// -----------------------------------------------------------------------------
// 运行时核心工具：状态选择与速度计算
// -----------------------------------------------------------------------------

void Band::SelectAnalyticKState(Particle* p, double E_target) {
    int itab = (int)((E_target - emin) * dlist);
    if (itab < 0) itab = 0;
    if (itab >= MWLE_ana) itab = MWLE_ana - 1;

    int start_index = analytic_ptlist[itab];
    int count = analytic_ntlist[itab];

    if (count <= 0) {
        return;
    }

    double total_weight = 0.0;
    for (int i = 0; i < count; ++i) {
        int grid_idx = analytic_tlist[start_index + i];
        total_weight += analytic_k_grid[grid_idx].weight;
    }

    double r = Random() * total_weight;
    double current_weight = 0.0;
    int selected_grid_idx = -1;

    for (int i = 0; i < count; ++i) {
        int grid_idx = analytic_tlist[start_index + i];
        current_weight += analytic_k_grid[grid_idx].weight;
        if (current_weight >= r) {
            selected_grid_idx = grid_idx;
            break;
        }
    }

    if (selected_grid_idx == -1) {
        selected_grid_idx = analytic_tlist[start_index + count - 1];
    }

    const AnalyticKPoint& pt = analytic_k_grid[selected_grid_idx];
    p->kx = pt.kx;
    p->ky = pt.ky;
    p->kz = pt.kz;

    // 直接反推索引（flat_idx -> kx/ky/kz_idx）
    int N = num_ticks_axis;
    p->kz_idx = selected_grid_idx % N;
    int temp = selected_grid_idx / N;
    p->ky_idx = temp % N;
    p->kx_idx = temp / N;

    // 同步速度与能量
    analytic_vx = pt.vx;
    analytic_vy = pt.vy;
    analytic_vz = pt.vz;
    p->energy = pt.energy;
}

/*
void Band::GetAnalyticV(Particle* p) {
    static const double HBAR = 1.0545718e-34;
    static const double M0 = 9.1093837e-31;

    double ml = 0.916 * M0;
    double mt = 0.190 * M0;
    double alpha_eV = 0.5;

    double k_scale = 1.0 / spr0;
    double kx_real = p->kx * k_scale;
    double ky_real = p->ky * k_scale;
    double kz_real = p->kz * k_scale;

    double E_eV = p->energy * eV0;

    double abs_kx = std::fabs(kx_real);
    double abs_ky = std::fabs(ky_real);
    double abs_kz = std::fabs(kz_real);

    double K_valley = 0.85 * (2.0 * PI / 5.431e-10);

    double kl = 0.0;
    double kt_vec[2] = {0.0, 0.0};
    int axis = 0;

    if (abs_kx >= abs_ky && abs_kx >= abs_kz) {
        axis = 0;
        kl = (kx_real > 0) ? (kx_real - K_valley) : (kx_real + K_valley);
        kt_vec[0] = ky_real;
        kt_vec[1] = kz_real;
    } else if (abs_ky >= abs_kx && abs_ky >= abs_kz) {
        axis = 1;
        kl = (ky_real > 0) ? (ky_real - K_valley) : (ky_real + K_valley);
        kt_vec[0] = kx_real;
        kt_vec[1] = kz_real;
    } else {
        axis = 2;
        kl = (kz_real > 0) ? (kz_real - K_valley) : (kz_real + K_valley);
        kt_vec[0] = kx_real;
        kt_vec[1] = ky_real;
    }

    double prefactor = HBAR / (1.0 + 2.0 * alpha_eV * E_eV);

    double vl = prefactor * (kl / ml);
    double vt1 = prefactor * (kt_vec[0] / mt);
    double vt2 = prefactor * (kt_vec[1] / mt);

    double v_norm_scale = 1.0 / velo0;

    if (axis == 0) {
        analytic_vx = vl * v_norm_scale;
        analytic_vy = vt1 * v_norm_scale;
        analytic_vz = vt2 * v_norm_scale;
    } else if (axis == 1) {
        analytic_vx = vt1 * v_norm_scale;
        analytic_vy = vl * v_norm_scale;
        analytic_vz = vt2 * v_norm_scale;
    } else {
        analytic_vx = vt1 * v_norm_scale;
        analytic_vy = vt2 * v_norm_scale;
        analytic_vz = vl * v_norm_scale;
    }
}
*/

void Band::ReadAnalyticData(string input_path) {
    cout << "Reading Analytic Band Data (Direct Lookup Mode)..." << endl;

    // 1) 读 DOS 表
    const string file_suffix = this->igzofl ? "_IGZO" : "";
    string dos_file = input_path + "/analytic_dos" + file_suffix + ".txt";
    ifstream in_dos(dos_file.c_str());
    if (!in_dos) {
        cerr << "Error: Cannot open " << dos_file << endl;
        exit(1);
    }
    if (mpi_rank == 0) cout << "  [Band] Reading DOS data from: " << dos_file << endl;
    char buffer[256];
    in_dos.getline(buffer, 256); // 跳过表头

    double E_eV, dos_real, dos_norm_val;
    int band_idx = bandof[PELEC];
    for(int i=0; i<=MTAB; i++) {
        dos[band_idx][i] = 0.0;
        sumdos[i][PELEC] = 0.0;
    }
    DOSMAX[PELEC] = 0.0;

    while(in_dos >> E_eV >> dos_real >> dos_norm_val) {
        double E_norm = E_eV / eV0;
        int itab = (int)((E_norm - emin) / dtable + 0.5);
        if (itab >= 0 && itab <= MTAB) {
            dos[band_idx][itab] = dos_norm_val;
            sumdos[itab][PELEC] = dos_norm_val;
            if (dos_norm_val > DOSMAX[PELEC]) DOSMAX[PELEC] = dos_norm_val;
        }
    }
    in_dos.close();
    cout << "  DOS table loaded successfully." << endl;

    // 2) 读 E-k-v 表，直接使用文件中给出的速度矢量
    string filename = input_path + "/analytic_ek" + file_suffix + ".txt";
    if (mpi_rank == 0) cout << "  [Band] Reading E-k data from: " << filename << endl;
    ifstream in_ek(filename.c_str());
    if (!in_ek) {
        cerr << "Error: Cannot open " << filename << endl;
        exit(1);
    }
    in_ek.getline(buffer, 256); // 跳过表头

    analytic_k_grid.clear();
    analytic_k_grid.reserve(8000000);

    // 归一化转换
    double a_lattice = sia0;
    double k_pi_to_internal = (PI / a_lattice);        // k_internal = k_pi * (pi/a_code)
    double v_real_to_internal = 1.0 / velo0;           // v_internal = v_real / velo0

    double kx_pi, ky_pi, kz_pi, E_eV_in;
    double vx_si, vy_si, vz_si;

    while(in_ek >> kx_pi >> ky_pi >> kz_pi >> E_eV_in >> vx_si >> vy_si >> vz_si) {
        AnalyticKPoint pt;

        // k、能量归一化
        pt.kx = kx_pi * k_pi_to_internal;
        pt.ky = ky_pi * k_pi_to_internal;
        pt.kz = kz_pi * k_pi_to_internal;
        pt.energy = E_eV_in / eV0;

        // 速度直接表读 -> 归一化
        pt.vx = vx_si * v_real_to_internal;
        pt.vy = vy_si * v_real_to_internal;
        pt.vz = vz_si * v_real_to_internal;
        pt.velocity = std::sqrt(pt.vx * pt.vx + pt.vy * pt.vy + pt.vz * pt.vz);

        // 权重
        double step_x = GetGridStep(kx_pi);
        double step_y = GetGridStep(ky_pi);
        double step_z = GetGridStep(kz_pi);
        pt.weight = step_x * step_y * step_z;

        // 简单谷索引（仅作标记）
        if (std::fabs(kx_pi) > 1.5) pt.valley_index = 0;
        else if (std::fabs(ky_pi) > 1.5) pt.valley_index = 1;
        else pt.valley_index = 2;

        analytic_k_grid.push_back(pt);
    }
    in_ek.close();

    cout << "  E-k-v table loaded. Total points: " << analytic_k_grid.size() << endl;
}
