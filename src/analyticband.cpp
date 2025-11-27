#include "Band.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

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
