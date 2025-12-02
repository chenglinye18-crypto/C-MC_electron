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

// 生成非均匀网格刻度 (单位: pi/a)
// 范围: [-2, 2], 重点加密波谷位置 +/- 1.7 (即 0.85 * 2)
std::vector<double> GenerateNonUniformTicks() {
    std::vector<double> ticks;
    double k_start = -2.0;
    double k_end = 2.0;
    double current = k_start;
    
    // 波谷中心位置 (归一化为 pi/a 单位)
    double valley_loc1 = 1.7; 
    double valley_loc2 = 0;
    double fine_region_width = 0.3; // 加密区域半径
    
    double step_coarse = 0.2; // 粗网格步长
    double step_fine = 0.02;  // 细网格步长

    while (current <= k_end + 1e-5) {
        ticks.push_back(current);
        
        double abs_k = std::fabs(current);
        // 判断是否在波谷附近
        bool near_valley1 = std::fabs(abs_k - valley_loc1) < fine_region_width;
        bool near_valley2 = std::fabs(abs_k - valley_loc2) < fine_region_width;
        
        if (near_valley1 || near_valley2) {
            current += step_fine;
        } else {
            current += step_coarse;
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

    // 1. 物理常数 (SI)
    double m0 = 9.10938356e-31; 
    double hbar = 1.0545718e-34; 
    double q = 1.60217662e-19;   
    double a_lattice = 5.43e-10; 

    // 参数
    double alpha_real = 0.5; // 1/eV
    double ml_kg = ml_rel * m0;
    double mt_kg = mt_rel * m0;

    // 波谷配置
    double K_valley_norm = 1.7; 
    double valley_dirs[6][3] = {
        {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}
    };
    int valley_l_axis[6] = {0, 0, 1, 1, 2, 2};

    // ----------------------------------------------------------
    // 2. 生成 E-k-v 表
    // ----------------------------------------------------------
    string ek_file = input_path + "/analytic_ek.txt";
    ofstream out_ek(ek_file.c_str());
    out_ek << "kx(pi/a) ky(pi/a) kz(pi/a) Energy(eV) vx(m/s) vy(m/s) vz(m/s)" << endl;

    vector<double> ticks = GenerateNonUniformTicks();
    int num_ticks = ticks.size();
    cout << "  Generating E-k-v table with " << num_ticks << "^3 points..." << endl;

    double k_conversion = PI / a_lattice; // pi/a -> 1/m
    double J_to_eV = 1.0 / q;

    for (int i = 0; i < num_ticks; i++) {
        for (int j = 0; j < num_ticks; j++) {
            for (int k = 0; k < num_ticks; k++) {
                double gx = ticks[i];
                double gy = ticks[j];
                double gz = ticks[k];

                // 1) 最近波谷
                int best_valley = 0;
                double min_dist_sq = 1.0e99;
                double dx_best = 0.0, dy_best = 0.0, dz_best = 0.0;
                for (int v = 0; v < 6; v++) {
                    double cx = valley_dirs[v][0] * K_valley_norm;
                    double cy = valley_dirs[v][1] * K_valley_norm;
                    double cz = valley_dirs[v][2] * K_valley_norm;
                    double dist_sq = (gx-cx)*(gx-cx) + (gy-cy)*(gy-cy) + (gz-cz)*(gz-cz);
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        best_valley = v;
                        dx_best = gx - cx;
                        dy_best = gy - cy;
                        dz_best = gz - cz;
                    }
                }

                // 2) 局部坐标 (真实单位)
                double dkx_real = dx_best * k_conversion;
                double dky_real = dy_best * k_conversion;
                double dkz_real = dz_best * k_conversion;

                double kl = 0.0, kt1 = 0.0, kt2 = 0.0;
                int l_axis = valley_l_axis[best_valley];
                if (l_axis == 0) { // X
                    kl = dkx_real; kt1 = dky_real; kt2 = dkz_real;
                } else if (l_axis == 1) { // Y
                    kl = dky_real; kt1 = dkx_real; kt2 = dkz_real;
                } else { // Z
                    kl = dkz_real; kt1 = dkx_real; kt2 = dky_real;
                }

                // 3) 能量
                double Gamma_J = (hbar*hbar/2.0) * ( (kl*kl)/ml_kg + (kt1*kt1 + kt2*kt2)/mt_kg );
                double Gamma_eV = Gamma_J * J_to_eV;
                double E_val = (-1.0 + sqrt(1.0 + 4.0 * alpha_real * Gamma_eV)) / (2.0 * alpha_real);

                // 4) 速度矢量
                double v_prefactor = hbar / (1.0 + 2.0 * alpha_real * E_val);
                double vl_vel  = v_prefactor * (kl  / ml_kg);
                double vt1_vel = v_prefactor * (kt1 / mt_kg);
                double vt2_vel = v_prefactor * (kt2 / mt_kg);

                double vx = 0.0, vy = 0.0, vz = 0.0;
                if (l_axis == 0) {       // X Valley
                    vx = vl_vel; vy = vt1_vel; vz = vt2_vel;
                } else if (l_axis == 1) { // Y Valley
                    vx = vt1_vel; vy = vl_vel; vz = vt2_vel;
                } else {                 // Z Valley
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
    // 3. 生成 DOS 表 (analytic_dos.txt) 并填充内部数组
    // ----------------------------------------------------------
    string dos_file = input_path + "/analytic_dos.txt";
    ofstream out_dos(dos_file.c_str());
    out_dos << "Energy(eV) DOS(1/eV/m^3) DOS_Norm(CodeUnits)" << endl;
    
    double max_E_eV = 7.5;
    
    double derived_eV0 = alpha_norm / 0.5; 
    
    for(int i=0; i<=MTAB; i++) {
        dos[bandof[PELEC]][i] = 0.0;
        sumdos[i][PELEC] = 0.0;
    }

    for (int itab = 0; itab <= MTAB; itab++) {
        double E_norm = energy[itab]; 
        
        double E_eV = E_norm * eV0; 
        
        if (E_eV > max_E_eV) break;
        
        double dos_real = CalculateAnalyticDOS_Real(E_eV, alpha_real, ml_rel, mt_rel, 6);
        
        double dos_code = dos_real * eV0 * std::pow(spr0, 3.0);
        
        dos[bandof[PELEC]][itab] = dos_code;
        sumdos[itab][PELEC] = dos_code;
        
        out_dos << E_eV << " " << dos_real << " " << dos_code << endl;
    }
    out_dos.close();
    
    DOSMAX[PELEC] = 0.0;
    for(int itab=0; itab<=MTAB; itab++) {
        if(sumdos[itab][PELEC] > DOSMAX[PELEC]) 
            DOSMAX[PELEC] = sumdos[itab][PELEC];
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
    string filename = input_path + "/phonon_dispersion.txt";
    cout << "Reading Phonon Spectrum from: " << filename << endl;

    ifstream infile(filename.c_str());
    if (!infile) {
        cerr << "Error: Cannot open phonon dispersion file!" << endl;
        exit(1);
    }

    string line;
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

        if (ss.fail()) break;

        for(int i=0; i<4; ++i) {
            phonon.omega_table[i].push_back(w[i]);
            phonon.vg_table[i].push_back(v[i]);
        }

    } while (std::getline(infile, line));

    infile.close();

    phonon.nq_tab = phonon.omega_table[0].size();
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

void Band::BuildAnalyticScatteringTable() {
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

            for (int iq = 0; iq < nq_int; ++iq) {
                double q = iq * dq_int;
                if (q < 1e-12) continue;
                if (q > 2.0 * ks) continue;

                double w_LA = GetPhononOmega(PH_LA, q);
                double w_TA = GetPhononOmega(PH_TA, q);

                if (w_LA <= 0 || w_TA <= 0) continue;

                double N_LA = 1.0 / (std::exp(HBAR_SI * w_LA / (KB_SI * T_lattice)) - 1.0);
                double N_TA = 1.0 / (std::exp(HBAR_SI * w_TA / (KB_SI * T_lattice)) - 1.0);

                double Iq = (Rs > 0) ? GetOverlapFactor(q, Rs) : 1.0;
                double q3_I2 = q * q * q * Iq * Iq;

                double hw_LA_eV = HBAR_SI * w_LA / Q_SI;
                double hw_TA_eV = HBAR_SI * w_TA / Q_SI;

                integ_LA += (1.0/w_LA) * N_LA * q3_I2;
                if (E_eV > hw_LA_eV) integ_LA += (1.0/w_LA) * (N_LA + 1.0) * q3_I2;

                integ_TA += (1.0/w_TA) * N_TA * q3_I2;
                if (E_eV > hw_TA_eV) integ_TA += (1.0/w_TA) * (N_TA + 1.0) * q3_I2;
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
    static const double a_lattice = 5.43e-10;
    double conversion = (PI / a_lattice) * spr0;
    k_ticks_code.resize(num_ticks_axis);
    for (int i = 0; i < num_ticks_axis; ++i) {
        k_ticks_code[i] = ticks[i] * conversion;
    }

    k_map_min = -2.1;
    k_map_max = 2.1;
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
    static const double a_lattice = 5.43e-10;
    double to_pi = 1.0 / ((PI / a_lattice) * spr0);

    double min_dt = 1.0e99;

    auto update_dt = [&](double k_val, double F, double &min_val) {
        if (std::fabs(F) < 1.0e-20) return;
        int idx = GetAxisIndex_O1((k_val + (F > 0 ? 1e-9 : -1e-9)) * to_pi);
        double wall;
        if (F > 0) {
            if (idx >= num_ticks_axis - 1) wall = 1.0e99;
            else wall = k_ticks_code[idx + 1];
        } else {
            if (idx <= 0) wall = -1.0e99;
            else wall = k_ticks_code[idx];
        }
        double dt = (wall - k_val) / F;
        if (dt > 0 && dt < min_val) min_val = dt;
    };

    update_dt(p->kx, Fx, min_dt);
    update_dt(p->ky, Fy, min_dt);
    update_dt(p->kz, Fz, min_dt);

    if (min_dt < 1.0e-8) min_dt = 1.0e-8;
    return min_dt;
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

    double a_lattice = 5.43e-10;
    double K_valley_mag = 0.85 * (2.0 * PI / a_lattice);
    double K0_code = K_valley_mag * spr0;

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

    GetAnalyticV_FromTable(p);
}

// -----------------------------------------------------------------------------
// 解析能带：声子散射与能谷操作模块
// -----------------------------------------------------------------------------

void Band::InitValleyConfiguration() {
    static const double a_lattice = 5.43e-10;
    double K_real = 0.85 * (2.0 * PI / a_lattice);
    valley_k0_norm = K_real * spr0;

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

    int current_valley = GetValleyID(p->kx, p->ky, p->kz);
    int target_valley = current_valley;
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
    GetAnalyticV_FromTable(p);
}
void Band::GetAnalyticV_FromTable(Particle* p) {
    static const double a_lattice = 5.43e-10;
    static const double conversion_factor = 1.0 / ((PI / a_lattice) * spr0);

    double kx_pi = p->kx * conversion_factor;
    double ky_pi = p->ky * conversion_factor;
    double kz_pi = p->kz * conversion_factor;

    int ix = GetAxisIndex_O1(kx_pi);
    int iy = GetAxisIndex_O1(ky_pi);
    int iz = GetAxisIndex_O1(kz_pi);

    int N = num_ticks_axis;
    int flat_idx = ix * (N * N) + iy * N + iz;
    if (flat_idx < 0 || flat_idx >= static_cast<int>(analytic_k_grid.size())) {
        analytic_vx = analytic_vy = analytic_vz = 0.0;
        return;
    }

    const AnalyticKPoint& pt = analytic_k_grid[flat_idx];
    analytic_vx = pt.vx;
    analytic_vy = pt.vy;
    analytic_vz = pt.vz;
    p->energy = pt.energy; // 保持 E-k 一致
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
    analytic_vx = pt.vx;
    analytic_vy = pt.vy;
    analytic_vz = pt.vz;
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
    string dos_file = input_path + "/analytic_dos.txt";
    ifstream in_dos(dos_file.c_str());
    if (!in_dos) {
        cerr << "Error: Cannot open " << dos_file << endl;
        exit(1);
    }
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
    string ek_file = input_path + "/analytic_ek.txt";
    ifstream in_ek(ek_file.c_str());
    if (!in_ek) {
        cerr << "Error: Cannot open " << ek_file << endl;
        exit(1);
    }
    in_ek.getline(buffer, 256); // 跳过表头

    analytic_k_grid.clear();
    analytic_k_grid.reserve(8000000);

    // 归一化转换
    double a_lattice = 5.43e-10;
    double k_pi_to_internal = (PI / a_lattice) * spr0; // k_internal = k_pi * (pi/a*spr0)
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
