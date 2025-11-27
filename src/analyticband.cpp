#include "Band.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cstdlib>

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
    
    double step_coarse = 0.1; // 粗网格步长
    double step_fine = 0.01;  // 细网格步长

    while (current <= k_end + 0) {
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
    std::cout << "Initializing Analytic Band (Kane's Model - Gamma-X Only)..." << std::endl;

    // 1. 物理常数准备 (用于生成 Output 表格，SI 单位)
    double m0 = 9.10938356e-31; // kg
    double hbar = 1.0545718e-34; // J*s
    double q = 1.60217662e-19;   // C
    double a_lattice = 5.43e-10; // Si 晶格常数 (m)

    // 恢复 alpha 到真实单位 (1/eV)
    double alpha_real = 0.5; 

    // 质量 (kg)
    double ml_kg = ml_rel * m0;
    double mt_kg = mt_rel * m0;

    // 波谷配置
    double K_valley_norm = 1.7; 
    double valley_dirs[6][3] = {
        {1,0,0}, {-1,0,0}, 
        {0,1,0}, {0,-1,0}, 
        {0,0,1}, {0,0,-1}  
    };
    int valley_l_axis[6] = {0, 0, 1, 1, 2, 2};

    // ----------------------------------------------------------
    // 2. 生成 E-k-v 表 (analytic_ek.txt)
    // ----------------------------------------------------------
    string ek_file = input_path + "/analytic_ek.txt";
    ofstream out_ek(ek_file.c_str());
    
    // 表头: 坐标(pi/a) 能量(eV) 速度模长(m/s)
    out_ek << "kx(pi/a) ky(pi/a) kz(pi/a) Energy(eV) Velocity(m/s)" << endl;
    
    // 生成非均匀网格
    std::vector<double> ticks = GenerateNonUniformTicks();
    int num_ticks = ticks.size();
    
    cout << "  Generating E-k table with grid size: " << num_ticks << "^3 points..." << endl;

    // 单位转换系数: k(pi/a) -> k(1/m)
    double k_conversion = PI / a_lattice;
    double J_to_eV = 1.0 / q;

    for (int i = 0; i < num_ticks; i++) {
        for (int j = 0; j < num_ticks; j++) {
            for (int k = 0; k < num_ticks; k++) {
                // 当前全局坐标 (pi/a)
                double gx = ticks[i];
                double gy = ticks[j];
                double gz = ticks[k];
                
                // --- 第一步：找到最近的波谷 ---
                int best_valley = 0;
                double min_dist_sq = 1.0e99;
                double dx_best=0, dy_best=0, dz_best=0; // 相对坐标 (pi/a)
                
                for (int v = 0; v < 6; v++) {
                    double cx = valley_dirs[v][0] * K_valley_norm;
                    double cy = valley_dirs[v][1] * K_valley_norm;
                    double cz = valley_dirs[v][2] * K_valley_norm;
                    
                    double dx = gx - cx;
                    double dy = gy - cy;
                    double dz = gz - cz;
                    
                    double dist_sq = dx*dx + dy*dy + dz*dz;
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        best_valley = v;
                        dx_best = dx; dy_best = dy; dz_best = dz;
                    }
                }
                
                // --- 第二步：转换到局部坐标系 (Longitudinal / Transversal) ---
                double dkx_real = dx_best * k_conversion;
                double dky_real = dy_best * k_conversion;
                double dkz_real = dz_best * k_conversion;
                
                double kl_real = 0.0; 
                double kt2_real = 0.0; 
                
                int l_axis = valley_l_axis[best_valley];
                
                if (l_axis == 0) { 
                    kl_real = dkx_real;
                    kt2_real = dky_real*dky_real + dkz_real*dkz_real;
                } else if (l_axis == 1) { 
                    kl_real = dky_real;
                    kt2_real = dkx_real*dkx_real + dkz_real*dkz_real;
                } else { 
                    kl_real = dkz_real;
                    kt2_real = dkx_real*dkx_real + dky_real*dky_real;
                }
                
                // --- 第三步：计算能量 E ---
                double Gamma_J = (hbar * hbar / 2.0) * ( (kl_real*kl_real)/ml_kg + kt2_real/mt_kg );
                double Gamma_eV = Gamma_J * J_to_eV;
                
                double E_val = (-1.0 + std::sqrt(1.0 + 4.0 * alpha_real * Gamma_eV)) / (2.0 * alpha_real);
                
                // --- 第四步：计算速度 v ---
                double term_l = (kl_real * kl_real) / (ml_kg * ml_kg);
                double term_t = kt2_real / (mt_kg * mt_kg);
                
                double v_prefactor = hbar / (1.0 + 2.0 * alpha_real * E_val); 
                double v_mag = v_prefactor * std::sqrt(term_l + term_t);
                
                // 输出数据
                out_ek << gx << " " << gy << " " << gz << " " << E_val << " " << v_mag << endl;
            }
        }
    }
    out_ek.close();
    cout << "  E-k table generated: " << ek_file << endl;

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

void Band::BuildAnalyticLists() {
    cout << "Building Analytic Lists (Indexing)..." << endl;
    if (dlist <= 0) dlist = 1.0 / (0.01 * eV0);

    analytic_ntlist.assign(MWLE, 0);
    analytic_ptlist.assign(MWLE, 0);
    analytic_tlist.resize(analytic_k_grid.size());

    for (size_t i = 0; i < analytic_k_grid.size(); ++i) {
        double E = analytic_k_grid[i].energy;
        int itab = (int)((E - emin) * dlist);
        if (itab >= 0 && itab < MWLE) analytic_ntlist[itab]++;
    }

    int current_offset = 0;
    for (int itab = 0; itab < MWLE; ++itab) {
        analytic_ptlist[itab] = current_offset;
        current_offset += analytic_ntlist[itab];
    }

    std::vector<int> temp_counters(MWLE, 0);
    for (size_t i = 0; i < analytic_k_grid.size(); ++i) {
        double E = analytic_k_grid[i].energy;
        int itab = (int)((E - emin) * dlist);
        if (itab >= 0 && itab < MWLE) {
            int pos = analytic_ptlist[itab] + temp_counters[itab];
            analytic_tlist[pos] = (int)i;
            temp_counters[itab]++;
        }
    }
    cout << "  Analytic lists built. Indexed " << analytic_k_grid.size() << " states into " << MWLE << " energy bins." << endl;
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
    for(int it = 0; it < nt; it++) {
        gamtet[it] = max_gamma;
    }
    gamma[PELEC] = max_gamma;

    cout << "  Analytic scattering table built. Max Rate (norm) = " << max_gamma << endl;
}

void Band::BuildAnalyticInjectionTable() {
    cout << "Building Analytic Injection Table (Output to file)..." << endl;

    double Ef_min_eV = -0.4;
    double Ef_max_eV = 0.4;
    int n_Ef_steps = 10000;

    double T_lattice = T0; 
    double kb_T_eV = T_lattice * 8.617333262e-5; 
    double q_SI = 1.602176634e-19;

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

        double g_norm = sumdos[itab][bandof[PELEC]];
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
         << "n= " << Ef_density[mid]*conc0 << " m^-3" << endl;
}

void Band::ReadAnalyticData(string input_path) {
    cout << "Reading Analytic Band Data..." << endl;

    // 1) 读 DOS
    string dos_file = input_path + "/analytic_dos.txt";
    ifstream in_dos(dos_file.c_str());
    if (!in_dos) {
        cerr << "Error: Cannot open " << dos_file << endl;
        exit(1);
    }
    char buffer[256];
    in_dos.getline(buffer, 256); // skip header

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

    // 2) 读 E-k-v
    string ek_file = input_path + "/analytic_ek.txt";
    ifstream in_ek(ek_file.c_str());
    if (!in_ek) {
        cerr << "Error: Cannot open " << ek_file << endl;
        exit(1);
    }
    in_ek.getline(buffer, 256); // skip header

    analytic_k_grid.clear();
    analytic_k_grid.reserve(2000000);

    double kx_pi, ky_pi, kz_pi, v_ms;
    double a_lattice = 5.43e-10;
    double k_to_internal = (PI / a_lattice) * spr0;

    while(in_ek >> kx_pi >> ky_pi >> kz_pi >> E_eV >> v_ms) {
        AnalyticKPoint pt;
        pt.kx = kx_pi * k_to_internal;
        pt.ky = ky_pi * k_to_internal;
        pt.kz = kz_pi * k_to_internal;
        pt.energy = E_eV / eV0;
        pt.velocity = v_ms / velo0;

        double step_x = GetGridStep(kx_pi);
        double step_y = GetGridStep(ky_pi);
        double step_z = GetGridStep(kz_pi);
        pt.weight = step_x * step_y * step_z;
        pt.valley_index = 0;

        analytic_k_grid.push_back(pt);
    }
    in_ek.close();
    cout << "  E-k table loaded. Total points: " << analytic_k_grid.size() << endl;
}
