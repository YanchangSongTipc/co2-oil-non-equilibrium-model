#include "saft_mie.h" // 自己的头文件放最上面

// C++ 标准库
#include <cmath>      // C++推荐用 <cmath> 替代 "math.h"
#include <iostream>
#include <chrono>
#include <ctime>
#include <algorithm>

// 引入 autodiff 魔法
#include <autodiff/forward/dual.hpp>

// 在 .cpp 文件中放心使用 using，省去敲 std:: 的麻烦
using namespace autodiff;
using std::vector;
using std::string;
using std::vector;


void init_saft_cross_params(add_args& cppargs){
    // Determine number of components from the size of the parameter arrays
    int ncomp = cppargs.m.size();
    if (ncomp == 0) return; // Prevent errors if struct is empty
    int n_pairs = ncomp * ncomp;
    // Allocate memory and initialize to zero
    cppargs.e_ij.assign(n_pairs, 0.0);
    cppargs.s_ij.assign(n_pairs, 0.0);
    cppargs.c_ij.assign(n_pairs, 0.0);
    cppargs.lambda_aij.assign(n_pairs, 0.0);
    cppargs.lambda_rij.assign(n_pairs, 0.0);
    int idx = -1;
    for (int i = 0; i < ncomp; i++) {
        for (int j = 0; j < ncomp; j++) {
            idx++;
            
            // 1. Calculate Segment Diameter cross-term (s_ij)
            if (cppargs.l_ij.empty()) {
                cppargs.s_ij[idx] = (cppargs.s[i] + cppargs.s[j]) * 0.5;
            } else {
                cppargs.s_ij[idx] = (cppargs.s[i] + cppargs.s[j]) * 0.5 * (1.0 - cppargs.l_ij[idx]);
            } 
            
            // 2. Calculate Dispersion Energy cross-term (e_ij)
            // Optimize cubic powers with direct multiplication
            double s_i3 = cppargs.s[i] * cppargs.s[i] * cppargs.s[i];
            double s_j3 = cppargs.s[j] * cppargs.s[j] * cppargs.s[j];
            double s_ij3 = cppargs.s_ij[idx] * cppargs.s_ij[idx] * cppargs.s_ij[idx];
            
            bool calculate_e_ij = true;
            // For ions: e_ij is kept at zero to avoid dispersion between like ions (cations-cations or anions-anions)
            if (!cppargs.z.empty() && (cppargs.z[i] * cppargs.z[j] > 0)) {
                calculate_e_ij = false; 
            }
            
            if (calculate_e_ij) {
                double kij = cppargs.k_ij.empty() ? 0.0 : cppargs.k_ij[idx];
                cppargs.e_ij[idx] = (1.0 - kij) * std::sqrt(s_i3 * s_j3) / s_ij3 * std::sqrt(cppargs.e[i] * cppargs.e[j]);
            }
            
            // 3. Calculate Repulsive/Attractive Exponent cross-terms (lambda_aij, lambda_rij)
            cppargs.lambda_aij[idx] = 3.0 + std::sqrt((cppargs.lambda_a[i] - 3.0) * (cppargs.lambda_a[j] - 3.0));
            cppargs.lambda_rij[idx] = 3.0 + std::sqrt((cppargs.lambda_r[i] - 3.0) * (cppargs.lambda_r[j] - 3.0));
            
            // 4. Calculate pre-factor c_ij
            double l_r = cppargs.lambda_rij[idx];
            double l_a = cppargs.lambda_aij[idx];
            cppargs.c_ij[idx] = (l_r / (l_r - l_a)) * std::pow(l_r / l_a, l_a / (l_r - l_a));
        }
    }
}

template <typename T>
T saft_ares_cpp(T t, T rho, const vector<double>& x, const add_args& cppargs) {
    int ncomp = x.size(); 
    vector<T> d(ncomp); // 直径现在是温度的函数，必须为 T

    for (int i = 0; i < ncomp; i++) {
        double c = cppargs.lambda_r[i] / (cppargs.lambda_r[i] - cppargs.lambda_a[i]) * pow(cppargs.lambda_r[i] / cppargs.lambda_a[i], cppargs.lambda_a[i] / (cppargs.lambda_r[i] - cppargs.lambda_a[i]));
        T theta = c * cppargs.e[i] / t;
        d[i] = saft_gl10_int(theta, 1., cppargs.e[i], cppargs.s[i], t, cppargs.lambda_r[i], cppargs.lambda_a[i]);
    }

    if (!cppargs.z.empty()) {
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z[i] != 0) {
                d[i] = cppargs.s[i] * 0.88; 
            }
        }
    }

    double summ = 0.;
    for (int i = 0; i < ncomp; i++) summ += cppargs.m[i] * x[i];

    vector<double> x_s(ncomp);
    for (int i = 0; i < ncomp; i++) x_s[i] = cppargs.m[i] * x[i] / summ;

    T den = rho * N_AV / 1.0e30;
    T den_s = den * summ;
    T pi_6_dens = PI / 6.0 * den_s;

    vector<T> zeta(4, 0.0);
    for (int j = 0; j < ncomp; j++) {
        T d1 = d[j], d2 = d1 * d1, d3 = d2 * d1;
        zeta[0] += x_s[j];
        zeta[1] += x_s[j] * d1;
        zeta[2] += x_s[j] * d2;
        zeta[3] += x_s[j] * d3;
    }
   
    for(int i=0; i<4; i++) zeta[i] *= pi_6_dens;
    
    T z3_comp = 1.0 - zeta[3];
    T zeta2_3 = zeta[2] * zeta[2] * zeta[2];
    T a_HS = 6. / PI / den_s * (
        ((zeta2_3 / (zeta[3] * zeta[3]) - zeta[0]) * log(z3_comp)) + 
        (3. * zeta[1] * zeta[2] / z3_comp) + 
        (zeta2_3 / (zeta[3] * z3_comp * z3_comp))
    );

    T zeta_x = 0.0, zx_ave = 0.0;
    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            idx++;
            T d_ij = (d[i] + d[j]) * 0.5;
            double s_ij = cppargs.s_ij[idx];
            double xs_ij = x_s[i] * x_s[j];
            zeta_x += pi_6_dens * xs_ij * (d_ij * d_ij * d_ij);
            zx_ave += pi_6_dens * xs_ij * (s_ij * s_ij * s_ij);
        }
    }

    T zx2 = zeta_x * zeta_x, zx3 = zx2 * zeta_x, zx4 = zx2 * zx2;
    T KHS_comp = 1.0 - zeta_x;
    T KHS = (KHS_comp * KHS_comp * KHS_comp * KHS_comp) / (1. + 4. * zeta_x + 4. * zx2 - 4. * zx3 + zx4);

    T zx_ave2 = zx_ave * zx_ave, zx_ave5 = zx_ave2 * zx_ave2 * zx_ave, zx_ave8 = zx_ave2 * zx_ave2 * zx_ave2 * zx_ave2;

    T a1 = 0., a2 = 0., a3 = 0.;
    vector<T> a1_ij(ncomp * ncomp, 0.0);
    idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            idx++;
            T d_ij = (d[i] + d[j]) * 0.5;
            double s_ij = cppargs.s_ij[idx], e_ij = cppargs.e_ij[idx], c_ij = cppargs.c_ij[idx];
            double lambda_aij = cppargs.lambda_aij[idx], lambda_rij = cppargs.lambda_rij[idx];

            T x0 = s_ij / d_ij;
            double xs_ij = x_s[i] * x_s[j];

            T as1 = saft_as_lambda(den_s, d_ij, e_ij, zeta_x, lambda_aij);
            T B1 = saft_B_lambda(den_s, d_ij, e_ij, zeta_x, lambda_aij, x0);
            T as2 = saft_as_lambda(den_s, d_ij, e_ij, zeta_x, lambda_rij);
            T B2 = saft_B_lambda(den_s, d_ij, e_ij, zeta_x, lambda_rij, x0);

            a1_ij[idx] = c_ij * (pow(x0, lambda_aij) * (as1 + B1) - pow(x0, lambda_rij) * (as2 + B2));
            a1 += xs_ij * a1_ij[idx];

            double alpha = c_ij * (1. / (lambda_aij - 3.) - 1. / (lambda_rij - 3.));
            T chi_ij = saft_fk_alpha(1, phi_k, alpha) * zx_ave + saft_fk_alpha(2, phi_k, alpha) * zx_ave5 + saft_fk_alpha(3, phi_k, alpha) * zx_ave8;

            T t_a = pow(x0, 2. * lambda_aij) * (saft_as_lambda(den_s, d_ij, e_ij, zeta_x, 2. * lambda_aij) + saft_B_lambda(den_s, d_ij, e_ij, zeta_x, 2. * lambda_aij, x0));
            T t_b = 2. * pow(x0, lambda_aij + lambda_rij) * (saft_as_lambda(den_s, d_ij, e_ij, zeta_x, lambda_aij + lambda_rij) + saft_B_lambda(den_s, d_ij, e_ij, zeta_x, lambda_aij + lambda_rij, x0));
            T t_c = pow(x0, 2. * lambda_rij) * (saft_as_lambda(den_s, d_ij, e_ij, zeta_x, 2. * lambda_rij) + saft_B_lambda(den_s, d_ij, e_ij, zeta_x, 2. * lambda_rij, x0));

            a2 += xs_ij * 0.5 * KHS * (1. + chi_ij) * e_ij * c_ij * c_ij * (t_a - t_b + t_c);
            a3 += -xs_ij * pow(e_ij, 3.) * saft_fk_alpha(4, phi_k, alpha) * zx_ave * exp(saft_fk_alpha(5, phi_k, alpha) * zx_ave + saft_fk_alpha(6, phi_k, alpha) * zx_ave2);
        }
    }
    T a_MONO = summ * (a_HS + a1 / t + a2 / (t * t) + a3 / (t * t * t));
    
// ================== CHAIN contribution ==================
    T a_CHAIN = 0.;
    
    // k0, k1, k2, k3 都依赖于 zeta_x 和 KHS_comp，所以必须是 T
    T k0 = -log(KHS_comp) + (42. * zeta_x - 39. * zx2 + 9. * zx3 - 2. * zx4) / (6. * KHS_comp * KHS_comp * KHS_comp);
    T k1 = (zx4 + 6. * zx2 - 12. * zeta_x) / (2. * KHS_comp * KHS_comp * KHS_comp);
    T k2 = -3. * zx2 / (8. * KHS_comp * KHS_comp);
    T k3 = (-zx4 + 3. * zx2 + 3. * zeta_x) / (6. * KHS_comp * KHS_comp * KHS_comp);

    for (int i = 0; i < ncomp; ++i) {
        int idx_chain = i * ncomp + i; 
        T d_ii = d[i];
        T x0 = cppargs.s_ij[idx_chain] / d_ii;
        T g_HS = exp(k0 + k1 * x0 + k2 * (x0 * x0) + k3 * (x0 * x0 * x0));

        // ⚠️ 注意：这里调用的 da1ddens_cpp 和 da2ddens_cpp 也必须被模板化！
        T da1dr = da1ddens_cpp(t, den_s, x_s, cppargs, idx_chain, d);
        T da2dr = da2ddens_cpp(t, den_s, x_s, cppargs, idx_chain, d);

        // 纯参数保持 double
        double e_ii = cppargs.e_ij[idx_chain];
        double c_ii = cppargs.c_ij[idx_chain];
        double lambda_aii = cppargs.lambda_aij[idx_chain];
        double lambda_rii = cppargs.lambda_rij[idx_chain];

        // 所有的 std::pow 和 std::exp 改为 pow 和 exp
        T g1 = 1. / (2. * PI * e_ii * (d_ii * d_ii * d_ii)) * (3. * da1dr
            - c_ii * lambda_aii * pow(x0, lambda_aii) * (saft_as_lambda(den_s, d_ii, e_ii, zeta_x, lambda_aii) + saft_B_lambda(den_s, d_ii, e_ii, zeta_x, lambda_aii, x0)) / den_s
            + c_ii * lambda_rii * pow(x0, lambda_rii) * (saft_as_lambda(den_s, d_ii, e_ii, zeta_x, lambda_rii) + saft_B_lambda(den_s, d_ii, e_ii, zeta_x, lambda_rii, x0)) / den_s);

        double alpha = c_ii * (1. / (lambda_aii - 3.) - 1. / (lambda_rii - 3.));
        T gamma = 10. * (-tanh(10. * (0.57 - alpha)) + 1.) * zx_ave * (exp(e_ii / t) - 1.) * exp(-6.7 * zx_ave - 8. * zx_ave2);
        T g2 = (1. + gamma) / (2. * PI * (e_ii * e_ii) * (d_ii * d_ii * d_ii)) * (3. * da2dr
            - e_ii * KHS * (c_ii * c_ii) * lambda_rii * pow(x0, 2. * lambda_rii) * (saft_as_lambda(den_s, d_ii, e_ii, zeta_x, 2. * lambda_rii) + saft_B_lambda(den_s, d_ii, e_ii, zeta_x, 2. * lambda_rii, x0)) / den_s
            + e_ii * KHS * (c_ii * c_ii) * (lambda_aii + lambda_rii) * pow(x0, lambda_aii + lambda_rii) * (saft_as_lambda(den_s, d_ii, e_ii, zeta_x, lambda_aii + lambda_rii) + saft_B_lambda(den_s, d_ii, e_ii, zeta_x, lambda_aii + lambda_rii, x0)) / den_s
            - e_ii * KHS * (c_ii * c_ii) * lambda_aii * pow(x0, 2. * lambda_aii) * (saft_as_lambda(den_s, d_ii, e_ii, zeta_x, 2. * lambda_aii) + saft_B_lambda(den_s, d_ii, e_ii, zeta_x, 2. * lambda_aii, x0)) / den_s);
        // printf("g_HS: %f, g1: %f, g2: %f\n",(double)g_HS,(double)g1,(double)g2);
        T g_Mie = g_HS * exp(e_ii / t * g1 / g_HS + (e_ii / t) * (e_ii / t) * g2 / g_HS);
        
        a_CHAIN += x[i] * (1. - cppargs.m[i]) * log(g_Mie);
    }
    // ================== ASSOC contribution ==================
    T a_ASSOC = 0.;
    if (!cppargs.e_assoc.empty()) {
        vector<T> delta_abij(ncomp * ncomp); // 声明为 T
        idx = -1;
        for (int i = 0; i < ncomp; ++i) {
            for (int j = 0; j < ncomp; ++j) {
                idx++;
                T d_ij = (d[i] + d[j]) * 0.5;
                double s_ij = cppargs.s_ij[idx];
                T d_term = d[i] * d[j] / (d[i] + d[j]);
                
                // 去掉所有的 std::
                T Fabij = exp(sqrt(cppargs.e_assoc[i] * cppargs.e_assoc[j]) / t) - 1.;
                T gHsdij = 1. / z3_comp + 3. * d_term * zeta[2] / (z3_comp * z3_comp) + 2. * (d_term * d_term) * (zeta[2] * zeta[2]) / (z3_comp * z3_comp * z3_comp);
                
                double rcij = (cppargs.rc[i] * cppargs.s[i] + cppargs.rc[j] * cppargs.s[j]) * 0.5;
                double rdij = 0.2 * (cppargs.s[i] + cppargs.s[j]);
                
                T Kabij = 4. * PI * (d_ij * d_ij) * (log((rcij + 2. * rdij) / d_ij)
                    * (6. * rcij * rcij * rcij + 18. * rcij * rcij * rdij - 24. * rdij * rdij * rdij)
                    + (rcij + 2. * rdij - d_ij) * (22. * rdij * rdij - 5. * rcij * rdij - 7. * rdij * d_ij - 8. * rcij * rcij + rcij * d_ij + d_ij * d_ij))
                    / (72. * rdij * rdij * s_ij * s_ij * s_ij);
                
                delta_abij[idx] = s_ij * s_ij * s_ij * Fabij * gHsdij * Kabij;
            }
        }
        
        int num_sites = 0;
        vector<int> iA; 
        for (auto it = cppargs.assoc_num.begin(); it != cppargs.assoc_num.end(); ++it) {
            num_sites += *it;
            for (int i = 0; i < *it; i++) {
                iA.push_back(it - cppargs.assoc_num.begin());
            }
        }
        
        // 摩尔分数不随 t 和 rho 变化，保持 double 即可
        vector<double> x_assoc(num_sites); 
        for (int i = 0; i < num_sites; i++) {
            x_assoc[i] = x[iA[i]];
        }
        
        // XA 和 XA_old 必须升级为 T，因为它们受到方程组迭代的直接影响
        vector<T> XA(num_sites, 0.), XA_old(num_sites, 0.);
        T error_norm = 1.;
        
        while (error_norm > 1e-6) {
            for (int i = 0; i < num_sites; ++i) {
                T tmp = 0;
                for (int j = 0; j < num_sites; ++j) {
                    int idx_assoc = iA[i] * ncomp + iA[j];
                    tmp += x_assoc[j] * XA_old[j] * delta_abij[idx_assoc];
                }
                XA[i] = 1. / (1. + den * tmp);
            }
            error_norm = 0;
            for (int i = 0; i < num_sites; ++i) { 
                error_norm += (XA[i] - XA_old[i]) * (XA[i] - XA_old[i]);
            }
            error_norm = sqrt(error_norm); // 去掉 std::
            XA_old = XA;
        }
        for (int i = 0; i < num_sites; ++i) {
            a_ASSOC += x_assoc[i] * (log(XA[i]) - 0.5 * XA[i] + 0.5); // 去掉 std::
        }
    }
    return a_MONO + a_CHAIN + a_ASSOC;
}

double daresdt_cpp(double t,double rho,vector<double> x, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=t*0.01;
    double fr=saft_ares_cpp(t+hh,rho,x,cppargs);
    double fl=saft_ares_cpp(t-hh,rho,x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_ares_cpp(t+hh,rho,x,cppargs);
        fl=saft_ares_cpp(t-hh,rho,x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
vector<double> tmp_daresdt_cpp(double t,double rho,vector<double> x, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=t*0.01;
    double fr=saft_ares_cpp(t+hh,rho,x,cppargs);
    double fl=saft_ares_cpp(t-hh,rho,x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_ares_cpp(t+hh,rho,x,cppargs);
        fl=saft_ares_cpp(t-hh,rho,x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    vector<double>d(3);
    d[0]=ans;
    d[1]=hh;
    d[2]=err;
    return d;
}
double dares2dt2_cpp(double t,double rho,vector<double> x, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=t*0.1;
    double fr=saft_ares_cpp(t+hh,rho,x,cppargs);
    double f=saft_ares_cpp(t,rho,x,cppargs);
    double fl=saft_ares_cpp(t-hh,rho,x,cppargs);
    a.push_back((fr-2*f+fl)/hh/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_ares_cpp(t+hh,rho,x,cppargs);
        f=saft_ares_cpp(t,rho,x,cppargs);
        fl=saft_ares_cpp(t-hh,rho,x,cppargs);
        a.push_back((fr-2*f+fl)/hh/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
double daresdrho_cpp(double t,double rho,vector<double> x, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=rho*0.01;
    double fr=saft_ares_cpp(t,rho+hh,x,cppargs);
    double fl=saft_ares_cpp(t,rho-hh,x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_ares_cpp(t,rho+hh,x,cppargs);
        fl=saft_ares_cpp(t,rho-hh,x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
/*Calculate the Fugacity Coefficient */
vector<double> saft_fugcoe_cpp(double t, double rho, vector<double> x, add_args& cppargs) {
    int ncomp = x.size();
    vector<double> fugcoe(ncomp,0);
    vector<double> partial_ax(ncomp, 0);
    double ares = saft_ares_cpp(t, rho, x, cppargs);
    double Z = saft_Z_cpp(t, rho, x, cppargs);
    double tmp = 0;
    for (int i = 0; i < ncomp; ++i) {
        vector<double> x_tmp = x;
        x_tmp[i] = x[i] + 1e-6;
        partial_ax[i] = (saft_ares_cpp(t, rho, x_tmp, cppargs) - ares) / 1e-6;
        tmp += x[i] * partial_ax[i];
    }
    for (int i = 0; i < ncomp; ++i) {
        fugcoe[i] = exp(ares + Z - 1 + partial_ax[i] - tmp - log(Z));
    }
    return fugcoe;
}
//inner function of the integral of hard-sphere system diameter d
double saft_func_lj(double e, double s, double t, double r,double lambda_r,double lambda_a ) {
    if (r > 0)
        return 1. - exp(-saft_lj_mie(e, s, r,lambda_r,lambda_a) / t);
    else
        return 1.;
    //return exp(r)/r;
}
//Lennard Jones Mie Potential
double saft_lj_mie(double e, double s, double r,double n,double m) {
    if (r > 0) {
        return n / (n - m) * pow(n / m, m / (n - m)) * e * (pow(s / r, n) - pow(s / r, m));
    }
    else
    {
        return 0;
    }
}
//derivative of inner function
double saft_func_dljdt(double e, double s, double t, double r,double lambda_r,double lambda_a) {
    if (r > 0)
        return -exp(-saft_lj_mie(e, s, r,lambda_r,lambda_a) / t) * saft_lj_mie(e, s, r,lambda_r,lambda_a) / t / t;
    else
        return 0.;
    //return exp(r)/r;
}
// 确保传入的 t 和 x 被声明为 T
template <typename T>
T saft_func_glf(double e, T t, T x, double lambda_r, double lambda_a) {
    // 1. 常数 c 纯由常量组成，可以直接用 double
    double c = lambda_r / (lambda_r - lambda_a) * pow(lambda_r / lambda_a, lambda_a / (lambda_r - lambda_a));
    // 2. theta 包含温度 t，必须是 T 类型
    T theta = c * e / t;
    // 3. 严格保留你原始的物理分段逻辑
    if (x > 0.0) {
        return (1.0 / lambda_r) * pow(x, -1.0 / lambda_r - 1.0) * exp(theta * pow(x, lambda_a / lambda_r));
    } else {
        return T(0.0);
    }
}
double saft_rbg_int(double a, double b, double err, double e, double s, double t, double lambda_r, double lambda_a) {
    double F = 0, h = 0;
    h = (b - a) / 2.;
    int n = 1, k = 1;
    vector<vector<double> >Ts;
    vector<double> line;
    Ts.push_back(line);
    Ts[0].push_back(h * (saft_func_lj(e, s, t, a,lambda_r,lambda_a) + saft_func_lj(e, s, t, b, lambda_r, lambda_a)));

    bool flag = false;
    while (!flag)
    {
        F = 0;
        for (int i = 0; i < n; ++i) {
            F += saft_func_lj(e, s, t, a + (2 * i + 1) * h, lambda_r, lambda_a);
        }
        Ts[0].push_back(Ts[0][k - 1] / 2. + h * F);
        Ts.push_back(line);
        for (int m = 1; m < k + 1; ++m) {
            Ts[m].push_back((pow(4, m) * Ts[m - 1][k - m + 1] - Ts[m - 1][k - m]) / (pow(4, m) - 1));
        }
        if (fabs(Ts[k][0] - Ts[k - 1][0]) < err) {
            break;
        }
        else
        {
            h = h / 2.;
            n = n * 2;
            k += 1;
        }
    }
    return Ts[k][0];
}
template <typename T>
T saft_gl5_int(T r, double a, double e, double s, T t, double lambda_r, double lambda_a){
    static constexpr std::array<double, 10> u={0.26356031971814109102031,1.41340305910651679221800,3.59642577104072208122300,7.08581000585883755692200,12.6408008442757826594300};
    static constexpr std::array<double, 10> w={0.5217556105828086524759,0.3986668110831759274500,7.5942449681707595390e-2,3.6117586799220484545e-3,2.3369972385776227891e-5};
    T k = exp(-r * a) / r;  // 注意：去掉了 std::
    T rinv = 1. / r;
    T res = 0;
    for(size_t i = 0; i < u.size(); ++i){
        T xi = u[i] * rinv + a;
        T fi = saft_func_glf(e, t, xi, lambda_r, lambda_a);
        res += w[i] * fi;
    }
    return s * (1. - res * k);
}
template <typename T>
T saft_gl10_int(T r, double a, double e, double s, T t, double lambda_r, double lambda_a) {
    static constexpr std::array<double, 10> u={0.13779347054049243, 0.7294545495031705, 1.808342901740316, 3.4014336978548996, 5.552496140063804, 8.330152746764497, 11.843785837900066, 16.279257831378104, 21.99658581198076, 29.92069701227389};
    static constexpr std::array<double, 10> w={0.30844111576502015, 0.40111992915527356, 0.2180682876118094, 0.062087456098677746, 0.0095015169751811, 0.0007530083885875388, 2.8259233495995656e-5, 4.2493139849626863e-7, 1.8395648239796308e-9, 9.911827219609008e-13};
    T k = exp(-r * a) / r;  // 注意：去掉了 std::
    T rinv = 1. / r;
    T res = 0;
    for(size_t i = 0; i < u.size(); ++i){
        T xi = u[i] * rinv + a;
        T fi = saft_func_glf(e, t, xi, lambda_r, lambda_a);
        res += w[i] * fi;
    }
    return s * (1. - res * k);
}

double saft_rbg_dintdt(double a, double b, double err, double e, double s, double t,double lambda_r,double lambda_a) {
    double F = 0, h = 0;
    h = (b - a) / 2.;
    int n = 1, k = 1;
    vector<vector<double> >Ts;
    vector<double> line;
    Ts.push_back(line);
    Ts[0].push_back(h * (saft_func_dljdt(e, s, t, a,lambda_r,lambda_a) + saft_func_dljdt(e, s, t, b, lambda_r, lambda_a)));

    bool flag = false;
    while (!flag)
    {
        F = 0;
        for (int i = 0; i < n; ++i) {
            F += saft_func_dljdt(e, s, t, a + (2 * i + 1) * h, lambda_r, lambda_a);
        }
        Ts[0].push_back(Ts[0][k - 1] / 2. + h * F);
        Ts.push_back(line);
        for (int m = 1; m < k + 1; ++m) {
            Ts[m].push_back((pow(4, m) * Ts[m - 1][k - m + 1] - Ts[m - 1][k - m]) / (pow(4, m) - 1));
        }
        if (fabs(Ts[k][0] - Ts[k - 1][0]) < err) {
            break;
        }
        else
        {
            h = h / 2.;
            n = n * 2;
            k += 1;
        }
    }
    return Ts[k][0];
}

//A12
template <typename T>
T saft_B_lambda(T den_s, T d, double e, T z, double lambda, T x) {
    T I_lambda = saft_I_lambda(lambda, x);
    T J_lambda = saft_J_lambda(lambda, x);
    // 【性能优化】：避免在 return 中重复计算 (1-z)^3
    T one_minus_z = 1.0 - z;
    T one_minus_z_3 = pow(one_minus_z, 3.0); 
    return 2. * PI * den_s * pow(d, 3.) * e * (I_lambda * (1. - z / 2.) / one_minus_z_3 - J_lambda * 9. * z * (1. + z) / 2. / one_minus_z_3);
}

// caculate the perturbation terms a^s for the corresponding Sutherland potentials. A16
template <typename T>
T saft_as_lambda(T den_s, T d, double e, T z, double lambda) {
    // 【性能优化】：提前计算 zeff，避免函数被调用两次，减少求导图的开销
    T zeff = saft_zeff_lambda(z, lambda);
    return -2. * den_s * (PI * e * pow(d, 3.) / (lambda - 3.)) * (1. - zeff / 2.) / pow(1. - zeff, 3.);
}

// A14
template <typename T>
T saft_I_lambda(double lambda, T x) {
    // 这里的 pow 没有 std::，autodiff 会自动接管 T 的求导
    return -(pow(x, 3. - lambda) - 1.) / (lambda - 3.);
}

// A15
template <typename T>
T saft_J_lambda(double lambda, T x) {
    return -(pow(x, 4. - lambda) * (lambda - 3.) - pow(x, 3. - lambda) * (lambda - 4.) - 1.) / ((lambda - 3.) * (lambda - 4.));
}

// caculate the effective packing fraction zete^eff_x. A17, A18
template <typename T>
T saft_zeff_lambda(T z, double lambda) {
    // 【性能优化 1】：lambda 是常数，提前算出 lambda 的 2 次方和 3 次方
    double lambda2 = lambda * lambda;
    double lambda3 = lambda2 * lambda;
    
    // c1 ~ c4 全程不涉及变量 z，所以它们依然是普通的 double
    double c1 = 0.81096 + 1.7888 / lambda - 37.578 / lambda2 + 92.284 / lambda3;
    double c2 = 1.0205 - 19.341 / lambda + 151.26 / lambda2 - 463.50 / lambda3;
    double c3 = -1.9057 + 22.845 / lambda - 228.14 / lambda2 + 973.92 / lambda3;
    double c4 = 1.0885 - 6.1962 / lambda + 106.98 / lambda2 - 677.64 / lambda3;
    
    // 【性能优化 2】：对变量 z 的多项式进行乘法展开，避免在 AD 求导图中使用缓慢的 pow
    T z2 = z * z;
    T z3 = z2 * z;
    T z4 = z2 * z2;
    
    return c1 * z + c2 * z2 + c3 * z3 + c4 * z4;
}
//caculate the function of fk in A22 and A25 with coefficients defined by phi. A26
double saft_fk_alpha(int k, vector<vector<double>> phi, double aplha) {
    double up = 0.;
    for (int i = 0; i < 4; ++i) {
        up += phi[k - 1] [i] * pow(aplha, i);
    }
    double down = 1.;
    for (int i = 4; i < 7; ++i) {
        down += phi[k - 1][i] * pow(aplha, i - 3.);
    }
    return up / down;
}
template <typename T>
T saft_a1_tmp(T t, T den_s, const vector<double>& x_s, const add_args& cppargs, int num, const vector<T>& d) {
    int ncomp = x_s.size();
    T pi_6_dens = PI / 6.0 * den_s;

    T zeta_x = 0.0;
    T zx_ave = 0.0;
    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            idx++;
            T d_ij = (d[i] + d[j]) * 0.5;
            double s_ij = cppargs.s_ij[idx];
            zeta_x += pi_6_dens * x_s[i] * x_s[j] * (d_ij * d_ij * d_ij);
            zx_ave += pi_6_dens * x_s[i] * x_s[j] * (s_ij * s_ij * s_ij);
        }
    }

    // 从一维索引 num 恢复出二维的 i 和 j
    int row = num / ncomp;
    int col = num % ncomp;
    T d_ij_num = (d[row] + d[col]) * 0.5;
    
    double s_ij = cppargs.s_ij[num];
    double e_ij = cppargs.e_ij[num];
    double c_ij = cppargs.c_ij[num];
    double lambda_aij = cppargs.lambda_aij[num];
    double lambda_rij = cppargs.lambda_rij[num];

    T x0 = s_ij / d_ij_num;

    T as1 = saft_as_lambda(den_s, d_ij_num, e_ij, zeta_x, lambda_aij);
    T B1 = saft_B_lambda(den_s, d_ij_num, e_ij, zeta_x, lambda_aij, x0);
    T as2 = saft_as_lambda(den_s, d_ij_num, e_ij, zeta_x, lambda_rij);
    T B2 = saft_B_lambda(den_s, d_ij_num, e_ij, zeta_x, lambda_rij, x0);

    // 去掉 std::
    return c_ij * (pow(x0, lambda_aij) * (as1 + B1) - pow(x0, lambda_rij) * (as2 + B2));
}

template <typename T>
T saft_a2_tmp(T t, T den_s, const vector<double>& x_s, const add_args& cppargs, const vector<T>& d) {
    int ncomp = x_s.size();
    T pi_6_dens = PI / 6.0 * den_s;

    T zeta_x = 0.0;
    T zx_ave = 0.0;
    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            idx++;
            T d_ij = (d[i] + d[j]) * 0.5;
            double s_ij = cppargs.s_ij[idx];
            zeta_x += pi_6_dens * x_s[i] * x_s[j] * (d_ij * d_ij * d_ij);
            zx_ave += pi_6_dens * x_s[i] * x_s[j] * (s_ij * s_ij * s_ij);
        }
    }

    T zx2 = zeta_x * zeta_x;
    T zx3 = zx2 * zeta_x;
    T zx4 = zx2 * zx2;
    T KHS_comp = 1.0 - zeta_x;
    T KHS = (KHS_comp * KHS_comp * KHS_comp * KHS_comp) / (1. + 4.*zeta_x + 4.*zx2 - 4.*zx3 + zx4);

    T zx_ave2 = zx_ave * zx_ave;
    T zx_ave5 = zx_ave2 * zx_ave2 * zx_ave;
    T zx_ave8 = zx_ave2 * zx_ave2 * zx_ave2 * zx_ave2;

    T a2_tmp = 0;
    idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            idx++;
            double c_ij = cppargs.c_ij[idx];
            double lambda_aij = cppargs.lambda_aij[idx];
            double lambda_rij = cppargs.lambda_rij[idx];
            double e_ij = cppargs.e_ij[idx];
            T d_ij = (d[i] + d[j]) * 0.5;
            double s_ij = cppargs.s_ij[idx];

            double alpha = c_ij * (1. / (lambda_aij - 3.) - 1. / (lambda_rij - 3.));
            
            // saft_fk_alpha 里面全是常数运算，不需要变成 T，直接返回 double 并与 T 相乘
            T chi_ij = saft_fk_alpha(1, phi_k, alpha) * zx_ave + saft_fk_alpha(2, phi_k, alpha) * zx_ave5 + saft_fk_alpha(3, phi_k, alpha) * zx_ave8;

            T x0 = s_ij / d_ij;
            T t_a = pow(x0, 2. * lambda_aij) * (saft_as_lambda(den_s, d_ij, e_ij, zeta_x, 2. * lambda_aij) + saft_B_lambda(den_s, d_ij, e_ij, zeta_x, 2. * lambda_aij, x0));
            T t_b = 2. * pow(x0, lambda_aij + lambda_rij) * (saft_as_lambda(den_s, d_ij, e_ij, zeta_x, lambda_aij + lambda_rij) + saft_B_lambda(den_s, d_ij, e_ij, zeta_x, lambda_aij + lambda_rij, x0));
            T t_c = pow(x0, 2. * lambda_rij) * (saft_as_lambda(den_s, d_ij, e_ij, zeta_x, 2. * lambda_rij) + saft_B_lambda(den_s, d_ij, e_ij, zeta_x, 2. * lambda_rij, x0));

            a2_tmp += x_s[i] * x_s[j] * 0.5 * KHS * (1. + chi_ij) * e_ij * c_ij * c_ij * (t_a - t_b + t_c);
        }
    }
    return a2_tmp;
}

template <typename T>
T saft_chi_tmp(T t, T den_s, const vector<double>& x_s, const add_args& cppargs, int num, const vector<T>& d) {
    int ncomp = x_s.size();
    T pi_6_dens = PI / 6.0 * den_s;

    T zx_ave = 0.0;
    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            idx++;
            double s_ij = cppargs.s_ij[idx];
            zx_ave += pi_6_dens * x_s[i] * x_s[j] * (s_ij * s_ij * s_ij);
        }
    }

    double alpha = cppargs.c_ij[num] * (1. / (cppargs.lambda_aij[num] - 3.) - 1. / (cppargs.lambda_rij[num] - 3.));
    
    T zx_ave2 = zx_ave * zx_ave;
    T zx_ave5 = zx_ave2 * zx_ave2 * zx_ave;
    T zx_ave8 = zx_ave2 * zx_ave2 * zx_ave2 * zx_ave2;

    return saft_fk_alpha(1, phi_k, alpha) * zx_ave + saft_fk_alpha(2, phi_k, alpha) * zx_ave5 + saft_fk_alpha(3, phi_k, alpha) * zx_ave8;
}

double from_eta_to_rho(double t, double rho, vector<double> x_s, double summ, add_args& cppargs,vector<double> d) {
    int ncomp = x_s.size(); // number of components

    //caculate the total number density of spherical segments.
    double den = rho * N_AV / 1.0e30;
    double den_s = den * summ;

    //caculate the moments of the number density zeta. A7
    vector<double> zeta(4, 0);
    for (int i = 0; i < 4; i++) {
        double sumd = 0;
        for (int j = 0; j < ncomp; j++) {
            sumd += x_s[j] * pow(d[j], i);
        }
        zeta[i] = PI / 6. * den_s * sumd;
    }
    double eta_tmp = zeta[3]+differential_coefficient;
    double sumd_3 = 0;
    for (int j = 0; j < ncomp; j++) {
        sumd_3 += x_s[j] * pow(d[j], 3.);
    }
    double dens_tmp = eta_tmp / PI * 6. / sumd_3;
    double rho_tmp = dens_tmp / summ / N_AV * 1.0e30;
    return rho_tmp;
}
double saft_A_cpp(double t, double rho, vector<double> x, add_args& cppargs) {
    double den = rho * N_AV;
    int ncomp = x.size();
    double a_IDEAL = -1.;
    for (int i = 0; i < ncomp; ++i) {
        double TdB_wavelength = h_planck / sqrt(2. * PI * cppargs.mol_mass[i] / 1000 / N_AV * kb * t);
        a_IDEAL += x[i] * log(den * x[i] * pow(TdB_wavelength, 3.));
    }
    double ares = saft_ares_cpp(t, rho, x, cppargs);
    //A=a*NkT
    double mass_ave = 0;
    for (int i = 0; i < ncomp; ++i) {
        mass_ave += x[i] * cppargs.mol_mass[i];
    }
    return (a_IDEAL + ares)* N_AV * kb * t;
}
double saft_Z_cpp(double t, double rho, const std::vector<double>& x, const add_args& cppargs) {
    // 将密度提升为自动微分变量
    autodiff::dual r = rho;
    double da_drho = autodiff::derivative(saft_ares_cpp<autodiff::dual>, wrt(r), at(t, r, x, cppargs));
    return 1.0 + rho * da_drho;
}
double dZdetat_cpp(double t, double eta, vector<double> x, add_args& cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=eta*0.01;
    double rho_r=eta_to_den(t,eta+hh,x,cppargs);
    double rho_l=eta_to_den(t,eta-hh,x,cppargs);
    double fr=saft_Z_cpp(t,rho_r*1e30/N_AV,x,cppargs);
    double fl=saft_Z_cpp(t,rho_l*1e30/N_AV,x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        rho_r=eta_to_den(t,eta+hh,x,cppargs);
        rho_l=eta_to_den(t,eta-hh,x,cppargs);
        fr=saft_Z_cpp(t,rho_r*1e30/N_AV,x,cppargs);
        fl=saft_Z_cpp(t,rho_l*1e30/N_AV,x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
double dZdtrho_cpp(double t, double rho, vector<double> x, add_args& cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=t*0.01;
    double fr=saft_Z_cpp(t+hh,rho,x,cppargs);
    double fl=saft_Z_cpp(t-hh,rho,x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_Z_cpp(t+hh,rho,x,cppargs);
        fl=saft_Z_cpp(t-hh,rho,x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
double saft_t_cpp(double p,double rho,std::vector<double> x,add_args cppargs){
    /* Newton Method, sometimes return NaN for zero point near 0, the gradient grows rapidly.
    double t_prime=1e3;
    double t=275.15;
    double den=rho*N_AV;
    int k=0;
    while(fabs(t-t_prime)>1e-4){
        t_prime=t;
        double f=support_func(t_prime,p,rho,x,cppargs);
        double fr=support_func(t_prime+1e-6,p,rho,x,cppargs);
        double fl=support_func(t_prime-1e-6,p,rho,x,cppargs);
        t=t_prime-f/((fr-fl)/2e-6);
    }
    return t;
    */
    //Bisection method, firstly, gget lower boundary, then search root space from lower boundary
    /*
    double lower=300.;
    bool run=true;
    
    double fl=support_func(lower,p,rho,x,cppargs);
    double tmp=lower,tmpf=0;
    while(tmpf==tmpf){
        lower=tmp;
        fl=tmpf;
        tmp=lower-10;
        tmpf=support_func(tmp,p,rho,x,cppargs);
    }
    double upper=lower+100;
    double fr=support_func(upper,p,rho,x,cppargs);
    while(fr*fl>0){
        upper+=100;
        fr=support_func(upper,p,rho,x,cppargs);
    }
    */
    //Bisection method, firstly, give upper boundary, then search root space from upper boundary
    double upper=800,lower=750.;
    double fr=support_func(upper,p,rho,x,cppargs);
    bool run=true;
    double fl=support_func(lower,p,rho,x,cppargs);
    while(fr*fl>0){
        upper=lower;
        fr=fl;
        lower-=50;
        fl=support_func(lower,p,rho,x,cppargs);
    }
    double f=support_func((upper+lower)/2.,p,rho,x,cppargs);
    int step=0;
    while (run&&step<15){
        ++step;
        if(fr*f<0){
            lower=(upper+lower)/2.;
            fl=f;
            f=support_func((upper+lower)/2.,p,rho,x,cppargs);
        }    
        else{
            upper=(upper+lower)/2.;
            fr=f;
            f=support_func((upper+lower)/2.,p,rho,x,cppargs);
        }
        if(fabs(f)<1){
            run=false;
        }
    }
    return (upper+lower)/2.;
}

double support_func(double t,double p,double rho,std::vector<double> x,add_args cppargs){
    double z=saft_Z_cpp(t,rho,x,cppargs);
    return p-z*t*kb*rho*N_AV;
}
double saft_newp_cpp(double t, double rho, vector<double> x, add_args& cppargs){
    double dadrho=daresdrho_cpp(t,rho,x,cppargs);
    return rho*kb*N_AV*t+dadrho*rho*rho;
}

double saft_p_cpp(double t, double rho, const std::vector<double>& x, const add_args& cppargs) {
    //Calculate pressure
    double Z = saft_Z_cpp(t, rho, x, cppargs);
    return rho * N_AV * kb * t * Z;
}
double saft_pv_cpp(double t,double v,vector<double> x, add_args &cppargs){
    double rho=1./v;
    return saft_p_cpp(t,rho,x,cppargs)/1e6;
}
double dpdv_cpp(double t,double v,vector<double> x, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=0.001*v;
    double fr=saft_pv_cpp(t,v+hh,x,cppargs);
    double fl=saft_pv_cpp(t,v-hh,x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_pv_cpp(t,v+hh,x,cppargs);
        fl=saft_pv_cpp(t,v-hh,x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}

double dp2dv2_cpp(double t,double v,vector<double> x, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=0.001*v;
    double fr=saft_pv_cpp(t,v+hh,x,cppargs);
    double fl=saft_pv_cpp(t,v-hh,x,cppargs);
    double f=saft_pv_cpp(t,v,x,cppargs);
    a.push_back((fr-2*f+fl)/hh/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_pv_cpp(t,v+hh,x,cppargs);
        fl=saft_pv_cpp(t,v-hh,x,cppargs);
        a.push_back((fr-2*f+fl)/hh/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
double dpdeta_cpp(double t,double eta,vector<double>x,add_args& cppargs){ 
    double errt,fac,hh,ans;
    vector<double> a;
    hh=eta*0.01;
    double rho_r=eta_to_den(t,eta+hh,x,cppargs);
    double rho_l=eta_to_den(t,eta-hh,x,cppargs);
    double fr=saft_p_cpp(t,rho_r*1e30/N_AV,x,cppargs);
    double fl=saft_p_cpp(t,rho_l*1e30/N_AV,x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        rho_r=eta_to_den(t,eta+hh,x,cppargs);
        rho_l=eta_to_den(t,eta-hh,x,cppargs);
        fr=saft_p_cpp(t,rho_r*1e30/N_AV,x,cppargs);
        fl=saft_p_cpp(t,rho_l*1e30/N_AV,x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
double dpdrho_cpp(double t,double rho,vector<double> x, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=rho*0.01;
    double fr=saft_p_cpp(t,rho+hh,x,cppargs);
    double fl=saft_p_cpp(t,rho-hh,x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_p_cpp(t,rho+hh,x,cppargs);
        fl=saft_p_cpp(t,rho-hh,x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
double dp2drho2_cpp(double t,double rho,vector<double> x, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=rho*0.01;
    double fr=saft_p_cpp(t,rho+hh,x,cppargs);
    double f=saft_p_cpp(t,rho,x,cppargs);
    double fl=saft_p_cpp(t,rho-hh,x,cppargs);
    a.push_back((fr-2*f+fl)/hh/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_p_cpp(t,rho+hh,x,cppargs);
        f=saft_p_cpp(t,rho,x,cppargs);
        fl=saft_p_cpp(t,rho-hh,x,cppargs);
        a.push_back((fr-2*f+fl)/hh/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
double dZdt_cpp(double t,double rho,vector<double> x, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=t*0.01;
    double fr=saft_Z_cpp(t+hh,rho,x,cppargs);
    double fl=saft_Z_cpp(t-hh,rho,x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_Z_cpp(t+hh,rho,x,cppargs);
        fl=saft_Z_cpp(t-hh,rho,x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
double dpdt_cpp(double t,double rho,vector<double> x, add_args &cppargs){
    double ans=rho*N_AV*kb*t*(dZdt_cpp(t,rho,x,cppargs)+saft_Z_cpp(t,rho,x,cppargs)/t);
    return ans;
}
//x0=t,x1=rho
double dpfuncdxi(double(*function)(double x1,double x2,vector<double> x, add_args &cppargs),
    vector<double> xi,vector<double> x,add_args &cppargs,int n){
    double errt,fac,hh,ans;
    vector<double> a,x_l=xi,x_r=xi;
    hh=xi[n]*0.01;
    if(n==0){
        x_l[0]=xi[0]-hh;
        x_r[0]=xi[0]+hh;
    }
    else{
        x_l[1]=xi[1]-hh;
        x_r[1]=xi[1]+hh;
    }
    double fr=function(x_r[0],x_r[1],x,cppargs);
    double fl=function(x_l[0],x_l[1],x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        if(n==0){
            x_l[0]=xi[0]-hh;
            x_r[0]=xi[0]+hh;
        }
        else{
            x_l[1]=xi[1]-hh;
            x_r[1]=xi[1]+hh;
        }
        fr=function(x_r[0],x_r[1],x,cppargs);
        fl=function(x_l[0],x_l[1],x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
double pfunction(double(*function)(double t,double rho,vector<double> x, add_args &cppargs),
    vector<double> xi,vector<double> x,add_args &cppargs){
    //x0=t,x1=rho
    double t=xi[0];
    double rho=xi[1];
    return function(t,rho,x,cppargs);
}
double saft_gres_cpp(double t, double rho, vector<double> x, add_args& cppargs) {
    /**
    Calculate the residual Gibbs energy for one phase of the system.
    */
    double ares = saft_ares_cpp(t, rho, x, cppargs);
    double Z = saft_Z_cpp(t, rho, x, cppargs);
    double gres = (ares + (Z - 1) - log(Z)) *kb *N_AV * t; // Equation A.50 from Gross and Sadowski 2001
    return gres;
}

double saft_rho_cpp(double t, double p, vector<double> x, int phase, add_args& cppargs) {
    /**
    Solve for the molar density when temperature and pressure are given.

    Parameters
    ----------
    t : double
        Temperature (K)
    p : double
        Pressure (Pa)
    phase : int
        The phase for which the calculation is performed. Options: 0 (liquid),
        1 (vapor).
    cppargs : add_args
        A struct containing additional arguments that can be passed for
        use in PC-SAFT:

        x : vector<double>, shape (n,)
            Mole fractions of each component. It has a length of n, where n is
            the number of components in the system.
        m : vector<double>, shape (n,)
            Segment number for each component.
        s : vector<double>, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : vector<double>, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : vector<double>, shape (n*n,)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : vector<double>, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : vector<double>, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : vector<double>, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : vector<double>, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Some implementations use this as an adjustable parameter
            that is fit to data.
        z : vector<double>, shape (n,)
            Charge number of the ions
        dielc : double
            Dielectric constant of the medium to be used for electrolyte
            calculations.

    Returns
    -------
    rho : double
        Molar density (mol m^-3)
    */
    // split into grid and find bounds for each root
    int ncomp = x.size(); // number of components
    vector<double> d(ncomp);
    for (int i = 0; i < ncomp; i++) {
        double c=cppargs.lambda_r[i]/(cppargs.lambda_r[i]-cppargs.lambda_a[i])*pow(cppargs.lambda_r[i]/cppargs.lambda_a[i],cppargs.lambda_a[i]/(cppargs.lambda_r[i]-cppargs.lambda_a[i]));
        double theta=c*cppargs.e[i]/t;
        d[i]=saft_gl10_int(theta,1.,cppargs.e[i],cppargs.s[i],t,cppargs.lambda_r[i],cppargs.lambda_a[i]);
    }
    if (!cppargs.z.empty()) {
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z[i] != 0) {
                d[i] = cppargs.s[i] * (1 - 0.12); // for ions the diameter is assumed to be temperature independent (see Held et al. 2014)
            }
        }
    }
    vector<double> x_lo, x_hi;
    int num_pts = 25;
    double err;
    double rho_guess = 1e-13;
    double rho_guess_prev = rho_guess;
    double err_prev = resid_rho(reduced_to_molar(rho_guess, t, ncomp, x, cppargs,d), t, p, x, cppargs);
    for (int i = 0; i < num_pts; i++) {
        rho_guess = 0.65 / (double)num_pts * i + 6e-3;
        err = resid_rho(reduced_to_molar(rho_guess, t, ncomp, x, cppargs,d), t, p, x, cppargs);
        if (err * err_prev < 0) {
            x_lo.push_back(rho_guess_prev);
            x_hi.push_back(rho_guess);
        }
        err_prev = err;
        rho_guess_prev = rho_guess;
    }

    // solve for appropriate root(s)
    double rho = _HUGE;
    double x_lo_molar, x_hi_molar;

    if (x_lo.size() == 1) {
        rho_guess = reduced_to_molar((x_lo[0] + x_hi[0]) / 2., t, ncomp, x, cppargs,d);
        x_lo_molar = reduced_to_molar(x_lo[0], t, ncomp, x, cppargs,d);
        x_hi_molar = reduced_to_molar(x_hi[0], t, ncomp, x, cppargs,d);
        rho = BrentRho(t, p, x, phase, cppargs, x_lo_molar, x_hi_molar, DBL_EPSILON, 1e-8, 200);
    }
    
    else if (x_lo.size() <= 3 && !x_lo.empty()) {
        if (phase == 0) {
            rho_guess = reduced_to_molar((x_lo.back() + x_hi.back()) / 2., t, ncomp, x, cppargs,d);
            x_lo_molar = reduced_to_molar(x_lo.back(), t, ncomp, x, cppargs,d);
            x_hi_molar = reduced_to_molar(x_hi.back(), t, ncomp, x, cppargs,d);
            rho = BrentRho(t, p, x, phase, cppargs, x_lo_molar, x_hi_molar, DBL_EPSILON, 1e-8, 200);
        }
        else if (phase == 1) {
            rho_guess = reduced_to_molar((x_lo[0] + x_hi[0]) / 400., t, ncomp, x, cppargs,d); // starting with a lower guess often provides better results
            x_lo_molar = reduced_to_molar(x_lo[0], t, ncomp, x, cppargs,d);
            x_hi_molar = reduced_to_molar(x_hi[0], t, ncomp, x, cppargs,d);
            rho = BrentRho(t, p, x, phase, cppargs, x_lo_molar, x_hi_molar, DBL_EPSILON, 1e-8, 200);
        }
    }
    
    else if (x_lo.size() > 3) {
        // if multiple roots to check, then find the one with the minimum gibbs energy. Reference: Privat R, Gani R, Jaubert JN. Are safe results obtained when the PC-SAFT equation of state is applied to ordinary pure chemicals?. Fluid Phase Equilibria. 2010 Aug 15;295(1):76-92.
        double g_min = 1e60;
        for (unsigned int i = 0; i < x_lo.size(); i++) {
            rho_guess = reduced_to_molar((x_lo[i] + x_hi[i]) / 2., t, ncomp, x, cppargs,d);
            x_lo_molar = reduced_to_molar(x_lo[i], t, ncomp, x, cppargs,d);
            x_hi_molar = reduced_to_molar(x_hi[i], t, ncomp, x, cppargs,d);
            double rho_i = BrentRho(t, p, x, phase, cppargs, x_lo_molar, x_hi_molar, DBL_EPSILON, 1e-8, 200);
            double g_i = saft_gres_cpp(t, rho_i, x, cppargs);
            if (g_i < g_min) {
                g_min = g_i;
                rho = rho_i;
            }
        }
    }
    else {
        int num_pts = 10;
        double err_min = 1e40;
        double rho_min = _HUGE;
        double err, rho_guess;
        for (int i = 0; i < num_pts; i++) {
            rho_guess = 0.65 / (double)num_pts * i + 1e-8;
            err = resid_rho(reduced_to_molar(rho_guess, t, ncomp, x, cppargs,d), t, p, x, cppargs);
            if (abs(err) < err_min) {
                err_min = abs(err);
                rho_min = reduced_to_molar(rho_guess, t, ncomp, x, cppargs,d);
            }
        }
        rho = rho_min;
    }

    return rho;
}

double resid_rho(double rhomolar, double t, double p, vector<double> x, add_args& cppargs) {
    double peos = saft_p_cpp(t, rhomolar, x, cppargs);
    double cost = (peos - p) / p;
    if (std::isfinite(cost)) {
        return cost;
    }
    else {
        return _HUGE;
    }
}

double reduced_to_molar(double nu, double t, int ncomp, vector<double> x, add_args& cppargs,vector<double> d) {
    double summ = 0.;
    for (int i = 0; i < ncomp; i++) {
        summ += x[i] * cppargs.m[i] * pow(d[i], 3.);
    }

    return 6 / PI * nu / summ * 1.0e30 / N_AV;
}

double BrentRho(double t, double p, vector<double> x, int phase, add_args& cppargs, double a, double b,
    double macheps, double tol_abs, int maxiter)
{
    int iter;
    double fa, fb, c, fc, m, tol, d, e, pp, q, s, r;
    fa = resid_rho(a, t, p, x, cppargs);
    fb = resid_rho(b, t, p, x, cppargs);

    // If one of the boundaries is to within tolerance, just stop
    if (std::abs(fb) < tol_abs) { return b; }
    if (std::isnan(fb)) {
        throw ValueError("BrentRho's method f(b) is NAN for b");
    }
    if (std::abs(fa) < tol_abs) { return a; }
    if (std::isnan(fa)) {
        throw ValueError("BrentRho's method f(a) is NAN for a");
    }
    if (fa * fb > 0) {
        throw ValueError("Inputs in BrentRho do not bracket the root");
    }

    c = a;
    fc = fa;
    iter = 1;
    if (std::abs(fc) < std::abs(fb)) {
        // Goto ext: from BrentRho ALGOL code
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
    }
    d = b - a;
    e = b - a;
    m = 0.5 * (c - b);
    tol = 2 * macheps * std::abs(b) + tol_abs;
    while (std::abs(m) > tol && fb != 0) {
        // See if a bisection is forced
        if (std::abs(e) < tol || std::abs(fa) <= std::abs(fb)) {
            m = 0.5 * (c - b);
            d = e = m;
        }
        else {
            s = fb / fa;
            if (a == c) {
                //Linear interpolation
                pp = 2 * m * s;
                q = 1 - s;
            }
            else {
                //Inverse quadratic interpolation
                q = fa / fc;
                r = fb / fc;
                m = 0.5 * (c - b);
                pp = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
                q = (q - 1) * (r - 1) * (s - 1);
            }
            if (pp > 0) {
                q = -q;
            }
            else {
                pp = -pp;
            }
            s = e;
            e = d;
            m = 0.5 * (c - b);
            if (2 * pp < 3 * m * q - std::abs(tol * q) || pp < std::abs(0.5 * s * q)) {
                d = pp / q;
            }
            else {
                m = 0.5 * (c - b);
                d = e = m;
            }
        }
        a = b;
        fa = fb;
        if (std::abs(d) > tol) {
            b += d;
        }
        else if (m > 0) {
            b += tol;
        }
        else {
            b += -tol;
        }
        fb = resid_rho(b, t, p, x, cppargs);
        if (std::isnan(fb)) {
            throw ValueError("BrentRho's method f(t) is NAN for t");
        }
        if (std::abs(fb) < macheps) {
            return b;
        }
        if (fb * fc > 0) {
            // Goto int: from BrentRho ALGOL code
            c = a;
            fc = fa;
            d = e = b - a;
        }
        if (std::abs(fc) < std::abs(fb)) {
            // Goto ext: from BrentRho ALGOL code
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        m = 0.5 * (c - b);
        tol = 2 * macheps * std::abs(b) + tol_abs;
        iter += 1;
        if (std::isnan(a)) {
            throw ValueError("BrentRho's method a is NAN");
        }
        if (std::isnan(b)) {
            throw ValueError("BrentRho's method b is NAN");
        }
        if (std::isnan(c)) {
            throw ValueError("BrentRho's method c is NAN");
        }
        if (iter > maxiter) {
            throw SolutionError("BrentRho's method reached maximum number of steps");
        }
        if (std::abs(fb) < 2 * macheps * std::abs(b)) {
            return b;
        }
    }
    return b;
}

double resid_inner(double R, double kb0, double Q, vector<double> u, vector<double> x, add_args &cppargs) {
    int ncomp = x.size();
    double error = 0;

    vector<double> pp(ncomp);
    double L = 0;
    for (int i = 0; i < ncomp; i++) {
        if (cppargs.z.empty() || cppargs.z[i] == 0) {
            pp[i] = x[i] / (1 - R + kb0 * R * exp(u[i]));
            L += pp[i];
        } else {
            L += x[i];
        }
    }
    L = (1 - R) * L;

    error = pow((L + Q - 1), 2.);
    return error;
}
/**
In the secant function, a 1-D Newton-Raphson solver is implemented.  An initial guess for the solution is provided.

@param x0 The initial guess for the solution
@param xmax The upper bound for the solution
@param xmin The lower bound for the solution
@param dx The initial amount that is added to x in order to build the numerical derivative
@param tol The absolute value of the tolerance accepted for the objective function
@param maxiter Maximum number of iterations
@returns If no errors are found, the solution, otherwise the value _HUGE, the value for infinity
*/
double BoundedSecantInner(double kb0, double Q, vector<double> u, vector<double> x, add_args &cppargs, double x0, double xmin,
    double xmax, double dx, double tol, int maxiter) {
    double x1=0,x2=0,x3=0,y1=0,y2=0,R,fval=999;
    int iter=1;
    if (std::abs(dx)==0){ throw ValueError("dx cannot be zero"); }
    while (std::abs(fval)>tol)
    {
        if (iter==1){
          x1 = x0;
          R = x1;
          x3 = R;
        }
        else if (iter==2){
          x2 = x0+dx;
          R = x2;
          x3 = R;
        }
        else {R=x2;}
        fval=resid_inner(R, kb0, Q, u, x, cppargs);

        if (iter==1){y1=fval;}
        else
        {
            if (std::isfinite(fval)) {
                y2 = fval;
            }
            else {
                y2 = 1e40;
            }
            x3=x2-y2/(y2-y1)*(x2-x1);
            // Check bounds, go half the way to the limit if limit is exceeded
            if (x3 < xmin)
            {
                x3 = (xmin + x2)/2;
            }
            if (x3 > xmax)
            {
                x3 = (xmax + x2)/2;
            }
            y1=y2; x1=x2; x2=x3;

        }
        if (iter>maxiter){
            throw SolutionError("BoundedSecant reached maximum number of iterations");
        }
        iter=iter+1;
    }
    return x3;
}
add_args get_single_component(int i, add_args &cppargs) {
    add_args args_single;
    args_single.m.push_back(cppargs.m[i]);
    args_single.s.push_back(cppargs.s[i]);
    args_single.e.push_back(cppargs.e[i]);
    args_single.lambda_a.push_back(cppargs.lambda_a[i]);
    args_single.lambda_r.push_back(cppargs.lambda_r[i]);
    if (!cppargs.e_assoc.empty()) {
        args_single.e_assoc.push_back(cppargs.e_assoc[i]);
        args_single.rc.push_back(cppargs.rc[i]);
    }
    if (!cppargs.dipm.empty()) {
        args_single.dipm.push_back(cppargs.dipm[i]);
        args_single.dip_num.push_back(cppargs.dip_num[i]);
    }
    if (!cppargs.z.empty()) {
        args_single.z.push_back(cppargs.z[i]);
        args_single.dielc = cppargs.dielc;
    }
    if (!cppargs.assoc_num.empty()) {
        args_single.assoc_num.push_back(cppargs.assoc_num[i]);

        if (args_single.assoc_num[0] > 0) {
            int nassoc = cppargs.assoc_num.size();
            int start = 0;
            for (int l = 0; l < (int)cppargs.assoc_num.size(); l++) {
                if (l < i) {
                    start += 1;
                }
            }
            for (int j = 0; j < nassoc; j++) {
                for (int k = 0; k < args_single.assoc_num[0]; k++) {
                    //args_single.assoc_matrix.push_back(cppargs.assoc_matrix[j*nassoc + start + k]);
                }
            }
        }
    }

    return args_single;
}

double estimate_flash_p(double t, double Q, vector<double> x, add_args &cppargs) {
    /**
    Get a quick estimate of the pressure at which VLE occurs
    */
    double p_guess;
    int ncomp = x.size();

    double x_ions = 0.; // overall mole fraction of ions in the system
    for (int i = 0; i < ncomp; i++) {
        if (!cppargs.z.empty() && cppargs.z[i] != 0) {
            x_ions += x[i];
        }
    }

    bool guess_found = false;
    double p_start = 10000;
    while (!guess_found && p_start < 1.5e7) {
        // initialize variables
        //std::cout<<"start"<<p_start<<std::endl;
        vector<double> fugcoef_l(ncomp), fugcoef_v(ncomp), k(ncomp), u(ncomp), kprime(ncomp);
        double rhol, rhov;
        double Pprime = 0.99 * p_start;
        double p = p_start;

        // calculate initial guess for compositions based on fugacity coefficients and Raoult's Law.
        rhol = saft_rho_cpp(t, p, x, 0, cppargs);
        rhov = saft_rho_cpp(t, p, x, 1, cppargs);
        if ((rhol - rhov) < 1e-4) {
            p_start = p_start + 2e5;
            continue;
        }
        fugcoef_l = saft_fugcoe_cpp(t, rhol, x, cppargs);
        fugcoef_v = saft_fugcoe_cpp(t, rhov, x, cppargs);

        vector<double> xl(ncomp);
        vector<double> xv(ncomp);
        double xv_sum = 0;
        double xl_sum = 0;
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                k[i] = fugcoef_l[i] / fugcoef_v[i];
            } else {
                k[i] = 0; // set k to 0 for ionic components
            }
            xl[i] = x[i] / (1 + Q * (k[i] - 1));
            xl_sum += xl[i];
            xv[i] = k[i] * x[i] / (1 + Q * (k[i] - 1));
            xv_sum += xv[i];
        }

        if (xv_sum != 1) {
            for (int i = 0; i < ncomp; i++) {
                xv[i] = xv[i] / xv_sum;
            }
        }

        if (xl_sum != 1) {
            for (int i = 0; i < ncomp; i++) {
                xl[i] = xl[i] / xl_sum;
            }
        }

        rhol = saft_rho_cpp(t, p, xl, 0, cppargs);
        rhov = saft_rho_cpp(t, p, xv, 1, cppargs);
        if ((rhol - rhov) < 1e-4) {
            p_start = p_start + 2e5;
            continue;
        }
        fugcoef_l = saft_fugcoe_cpp(t, rhol, xl, cppargs);
        fugcoef_v = saft_fugcoe_cpp(t, rhov, xv, cppargs);
        double numer = 0;
        double denom = 0;
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                numer += xl[i] * fugcoef_l[i];
                denom += xv[i] * fugcoef_v[i];
            }
        }
        double ratio = numer / denom;

        rhol = saft_rho_cpp(t, Pprime, xl, 0, cppargs);
        rhov = saft_rho_cpp(t, Pprime, xv, 1, cppargs);
        if ((rhol - rhov) < 1e-4) {
            p_start = p_start + 2e5;
            continue;
        }
        fugcoef_l = saft_fugcoe_cpp(t, rhol, xl, cppargs);
        fugcoef_v = saft_fugcoe_cpp(t, rhov, xv, cppargs);
        numer = 0;
        denom = 0;
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                numer += xl[i] * fugcoef_l[i];
                denom += xv[i] * fugcoef_v[i];
            }
        }
        double ratio_prime = numer / denom;

        double slope = (std::log10(ratio) - std::log10(ratio_prime)) / (std::log10(p) - std::log10(Pprime));
        double intercept = std::log10(ratio) - slope * std::log10(p);
        p_guess = pow(10, -intercept / slope);

        guess_found = true;
    }

    if (!guess_found) {
        throw SolutionError("an estimate for the VLE pressure could not be found");
    }

    return p_guess;
}

vector<double> outerTQ(double p_guess, double t, double Q, vector<double> x, add_args &cppargs) {
    // Based on the algorithm proposed in H. A. J. Watson, M. Vikse, T. Gundersen, and P. I. Barton, “Reliable Flash Calculations: Part 1. Nonsmooth Inside-Out Algorithms,” Ind. Eng. Chem. Res., vol. 56, no. 4, pp. 960–973, Feb. 2017, doi: 10.1021/acs.iecr.6b03956.
    int ncomp = x.size();
    double TOL = 1e-8;
    double MAXITER = 200;
    double x_ions = 0.; // overall mole fraction of ions in the system
    for (int i = 0; i < ncomp; i++) {
        if (!cppargs.z.empty() && cppargs.z[i] != 0) {
            x_ions += x[i];
        }
    }

    // initialize variables
    vector<double> fugcoef_l(ncomp), fugcoef_v(ncomp), k(ncomp), u(ncomp), kprime(ncomp);
    double rhol, rhov;
    double Pref = p_guess - 0.01 * p_guess;
    double Pprime = p_guess + 0.01 * p_guess;
    if (p_guess > 1e6) { // when close to the critical pressure then we need to have Pprime be less than p_guess
        Pprime = p_guess - 0.005 * p_guess;
    }
    double p = p_guess;

    // calculate initial guess for compositions based on fugacity coefficients and Raoult's Law.
    rhol = saft_rho_cpp(t, p, x, 0, cppargs);
    rhov = saft_rho_cpp(t, p, x, 1, cppargs);
    if ((rhol - rhov) > 1e-4) { // first, simply try calculating with the overall system composition, if the two density roots are found
        fugcoef_l = saft_fugcoe_cpp(t, rhol, x, cppargs);
        fugcoef_v = saft_fugcoe_cpp(t, rhov, x, cppargs);

        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                k[i] = fugcoef_l[i] / fugcoef_v[i];
            } else {
                k[i] = 0; // set k to 0 for ionic components
            }
        }
    } else {
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                try { // if using overall system composition doesn't work, try calculating using the vapor pressure
                    if (ncomp == 1) {
                        throw SolutionError("one component does not have a phase split at these conditions.");
                    }
                    vector<double> x_single(1, 1);
                    add_args args_single = get_single_component(i, cppargs);
                    //std::cout<<"single "<<i<<std::endl;
                    double Psat = flashTQ_cpp(t, Q, x_single, args_single)[0];
                    k[i] = Psat / p;
                }
                catch (const SolutionError& ex) { // if vapor pressure cannot be calculated for a component, then the component might be supercritical at these conditions. Instead, try calculating an initial guess at a very dilute concentration.
                    double DILUTE_MOLE_FRAC = 1e-4;
                    vector<double> x_dilute(ncomp, 0); // to make the calculation simpler, we only include the current component and the major component in the system when estimating an initial k value
                    int idx_largest = 0; // index of the major component in the system
                    for (int j = 0; j < ncomp; j++) {
                        if (x[j] > x[idx_largest]) {
                            idx_largest = j;
                        }
                    }
                    if(i==idx_largest){
                        x_dilute[i] = 1.;
                    }
                    else{
                        x_dilute[i] = DILUTE_MOLE_FRAC;
                        x_dilute[idx_largest] = 1 - DILUTE_MOLE_FRAC;
                    }
                    rhol = saft_rho_cpp(t, p, x_dilute, 0, cppargs);
                    //std::cout<<rhol<<std::endl;
                    if (rhol < 200) { // the density is probably not for the liquid phase
                        throw SolutionError("initial k could not be estimated in outerTQ for one or more components.");
                    }

                    fugcoef_l[i] = saft_fugcoe_cpp(t, rhol, x_dilute, cppargs)[i];
                    //std::cout<<"k="<<fugcoef_l[i]<<std::endl;
                    k[i] = fugcoef_l[i] / 1; // here we assume that the vapor phase behaves like an ideal gas (fugacity coefficient close to 1)
                }
            } else {
                k[i] = 0; // set k to 0 for ionic components
            }
        }
    }

    vector<double> xl(ncomp);
    vector<double> xv(ncomp);
    double xv_sum = 0;
    double xl_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        xl[i] = x[i] / (1 + Q * (k[i] - 1));
        xl_sum += xl[i];
        xv[i] = k[i] * x[i] / (1 + Q * (k[i] - 1));
        xv_sum += xv[i];
    }

    if (xv_sum != 1) {
        for (int i = 0; i < ncomp; i++) {
            xv[i] = xv[i] / xv_sum;
        }
    }

    if (xl_sum != 1) {
        for (int i = 0; i < ncomp; i++) {
            xl[i] = xl[i] / xl_sum;
        }
    }

    rhol = saft_rho_cpp(t, p, xl, 0, cppargs);
    fugcoef_l = saft_fugcoe_cpp(t, rhol, xl, cppargs);
    rhov = saft_rho_cpp(t, p, xv, 1, cppargs);
    fugcoef_v = saft_fugcoe_cpp(t, rhov, xv, cppargs);
    for (int i = 0; i < ncomp; i++) {
        k[i] = fugcoef_l[i] / fugcoef_v[i];
        u[i] = std::log(k[i] / kb);
    }

    rhol = saft_rho_cpp(t, Pprime, xl, 0, cppargs);
    fugcoef_l = saft_fugcoe_cpp(t, rhol, xl, cppargs);
    rhov = saft_rho_cpp(t, Pprime, xv, 1, cppargs);
    fugcoef_v = saft_fugcoe_cpp(t, rhov, xv, cppargs);
    for (int i = 0; i < ncomp; i++) {
        kprime[i] = fugcoef_l[i] / fugcoef_v[i];
    }

    vector<double> t_weight(ncomp);
    double t_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        double dlnk_dt = (kprime[i] - k[i]) / (Pprime - p);
        t_weight[i] = xv[i] * dlnk_dt / (1 + Q * (k[i] - 1));
        t_sum += t_weight[i];
    }

    double kb = 0;
    for (int i = 0; i < ncomp; i++) {
        double wi = t_weight[i] / t_sum;
        if (cppargs.z.empty() || cppargs.z[i] == 0) {
            kb += wi * std::log(k[i]);
        }
    }
    kb = std::exp(kb);

    t_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        double dlnk_dt = (kprime[i] - k[i]) / (Pprime - p);
        t_weight[i] = xv[i] * dlnk_dt / (1 + Q * (kprime[i] - 1));
        t_sum += t_weight[i];
    }

    double kbprime = 0;
    for (int i = 0; i < ncomp; i++) {
        double wi = t_weight[i] / t_sum;
        if (cppargs.z.empty() || cppargs.z[i] == 0) {
            kbprime += wi * std::log(kprime[i]);
        }
    }
    kbprime = std::exp(kbprime);
    double kb0 = kbprime;

    for (int i = 0; i < ncomp; i++) {
        u[i] = std::log(k[i] / kb);
    }

    double B = std::log(kbprime / kb) / (1/Pprime - 1/p);
    double A = std::log(kb) - B * (1/p - 1/Pref);
    if (B < 0) {
        //std::cout<<1<<std::endl;
        throw SolutionError("B < 0 in outerTQ");
    }

    // solve
    vector<double> u_old(ncomp);
    double A_old;
    vector<double> pp(ncomp);
    double maxdif = 1e10 * TOL;
    int itr = 0;
    double Rmin = 0, Rmax = 1;
    while (maxdif > TOL && itr < MAXITER) {
        // save previous values for calculating the difference at the end of the iteration
        u_old = u;
        std::copy(u.begin(), u.end(), u_old.begin());
        A_old = A;
        double R0 = kb * Q / (kb * Q + kb0 * (1 - Q));
        double R = BoundedSecantInner(kb0, Q, u, x, cppargs, R0, Rmin, Rmax, DBL_EPSILON, 1e-8, 200);

        double pp_sum = 0;
        double eupp_sum = 0;
        for (int i = 0; i < ncomp; i++) {
            pp[i] = x[i] / (1 - R + kb0 * R * std::exp(u[i]));
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                pp_sum += pp[i];
                eupp_sum += std::exp(u[i]) * pp[i];
            }
        }
        kb = pp_sum / eupp_sum;

        p = 1 / (1 / Pref + (std::log(kb) - A) / B);
        for (int i = 0; i < ncomp; i++) {
            if (x_ions == 0) {
                xl[i] = pp[i] / pp_sum;
                xv[i] = std::exp(u[i]) * pp[i] / eupp_sum;
            }
            else if (cppargs.z.empty() || cppargs.z[i] == 0) {
                xl[i] = pp[i] / pp_sum * (1 - x_ions/(1-Q));
                xv[i] = std::exp(u[i]) * pp[i] / eupp_sum;
            }
            else {
                xl[i] = x[i] / (1 - Q);
                xv[i] = 0;
            }
        }

        rhol = saft_rho_cpp(t, p, xl, 0, cppargs);
        fugcoef_l = saft_fugcoe_cpp(t, rhol, xl, cppargs);
        rhov = saft_rho_cpp(t, p, xv, 1, cppargs);
        fugcoef_v = saft_fugcoe_cpp(t, rhov, xv, cppargs);
        for (int i = 0; i < ncomp; i++) {
            k[i] = fugcoef_l[i] / fugcoef_v[i];
            u[i] = std::log(k[i] / kb);
        }

        if (itr == 0) {
            B = std::log(kbprime / kb) / (1/Pprime - 1/p);
            if (B < 0) {
                //std::cout<<1<<std::endl;
                throw SolutionError("B < 0 in outerTQ");
            }
        }
        A = std::log(kb) - B * (1./p - 1./Pref);

        maxdif = std::abs(A - A_old);
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                double dif = std::abs(u[i] - u_old[i]);
                if (dif > maxdif) {
                    maxdif = dif;
                } else if (!std::isfinite(dif)) {
                    maxdif = dif;
                }
            }
        }
        itr += 1;
    }

    if (!std::isfinite(p) || !std::isfinite(maxdif) || maxdif > 0.1 || p < 0) {
        //std::cout<<2<<std::endl;
        throw SolutionError("outerTQ did not converge to a solution");
    }

    //std::cout<<p_guess<<std::endl;
    vector<double> result;
    result.push_back(p);
    result.insert(result.end(), xl.begin(), xl.end());
    result.insert(result.end(), xv.begin(), xv.end());
    return result;
}

vector<double> flashTQ_cpp(double t, double Q, vector<double> x, add_args &cppargs) {
    bool solution_found = false;
    double p_guess;
    vector<double> output;
    try {
        p_guess = estimate_flash_p(t, Q, x, cppargs);
        output = outerTQ(p_guess, t, Q, x, cppargs);
        solution_found = true;
    }
    catch (const SolutionError& ex) {}
    catch (const ValueError& ex) {}
    int count=0;
    // if solution hasn't been found, try cycling through a range of pressures
    if (!solution_found) {
        double p_lbound = 5; // here we're using log10 of the pressure
        double p_ubound = 9;
        double p_step = 0.1;
        p_guess = p_lbound;
        
        while (p_guess < p_ubound && !solution_found) {
            count+=1;
            try {
                output = outerTQ(pow(10, p_guess), t, Q, x, cppargs);
                solution_found = true;
            } catch (const SolutionError& ex) {
                p_guess += p_step;
            } catch (const ValueError& ex) {
                p_guess += p_step;
            }
        }
    }

    if (!solution_found) {
        throw SolutionError("solution could not be found for TQ flash");
    }
    return output;
}


vector<double> flashTQ_cpp(double t, double Q, vector<double> x, add_args &cppargs, double p_guess) {
    vector<double> output;
    try {
        output = outerTQ(p_guess, t, Q, x, cppargs);
    } catch (const SolutionError& ex) {
        throw ex;
        output = flashTQ_cpp(t, Q, x, cppargs); // call function without an initial guess
    }

    return output;
}
double saft_hres_cpp(double t, double rho, vector<double> x, add_args &cppargs){
    double Z = saft_Z_cpp(t, rho, x, cppargs);
    double dares_dt=daresdt_cpp(t,rho,x,cppargs);
    double hres = (-t*dares_dt + (Z-1))*kb*N_AV*t; // Equation A.46 from Gross and Sadowski 2001
    return hres;
}

double dhresdt_cpp(double t,double rho,vector<double> x, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=t*0.01;
    double fr=saft_hres_cpp(t+hh,rho,x,cppargs);
    double fl=saft_hres_cpp(t-hh,rho,x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_hres_cpp(t+hh,rho,x,cppargs);
        fl=saft_hres_cpp(t-hh,rho,x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
//residual entropy J/(mol*k)
double saft_sres_cpp(double t, double rho, vector<double> x, add_args &cppargs) {
    /**
    Calculate the residual entropy (constant volume) for one phase of the system.
    */
    double gres = saft_gres_cpp(t, rho, x, cppargs);
    double hres = saft_hres_cpp(t, rho, x, cppargs);

    double sres = (hres - gres)/t;
    return sres;
}
//calculate density mol/m^3
double saft_density(double t,double p,vector<double> x,add_args &cppargs){
    double rhol=saft_rho_cpp(t,p,x,LIQUID,cppargs);
    double rhov=saft_rho_cpp(t,p,x,VAPER,cppargs);
    double gresl=saft_gres_cpp(t,rhol,x,cppargs);
    double gresv=saft_gres_cpp(t,rhov,x,cppargs);
    double rho=0;
    if(gresl<gresv){
        rho=rhol;
    }
    else{
        rho=rhov;
    }
    return rho;
}
double saft_drhodp(double t,double p,vector<double> x,add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=p*0.01;
    double fr=saft_density(t,p+hh,x,cppargs);
    double fl=saft_density(t,p-hh,x,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_density(t,p+hh,x,cppargs);
        fl=saft_density(t,p-hh,x,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
//from packing fraction calculate molecular number density
double eta_to_den(double t, double eta, vector<double> x, add_args &cppargs){
    int ncomp = x.size(); // number of components
    vector<double> d(ncomp);
    for (int i = 0; i < ncomp; i++) {
        double c=cppargs.lambda_r[i]/(cppargs.lambda_r[i]-cppargs.lambda_a[i])*pow(cppargs.lambda_r[i]/cppargs.lambda_a[i],cppargs.lambda_a[i]/(cppargs.lambda_r[i]-cppargs.lambda_a[i]));
        double theta=c*cppargs.e[i]/t;
        d[i]=saft_gl10_int(theta,1.,cppargs.e[i],cppargs.s[i],t,cppargs.lambda_r[i],cppargs.lambda_a[i]);
    }
    if (!cppargs.z.empty()) {
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z[i] != 0) {
                d[i] = cppargs.s[i] * (1 - 0.12); // for ions the diameter is assumed to be temperature independent (see Held et al. 2014)
            }
        }
    }
    //caculate the mole fraction of segments of component i. A8
    vector<double> x_s(ncomp);
    double summ = 0.;
    for (int i = 0; i < ncomp; i++) {
        summ += cppargs.m[i] * x[i];
    }
    for (int i = 0; i < ncomp; i++) {
        x_s[i] = cppargs.m[i] * x[i] / summ;
    }
    //caculate the moments of the number density zeta. A7
    vector<double> zeta(4, 0);
    for (int i = 0; i < 4; i++) {
        double sumd = 0;
        for (int j = 0; j < ncomp; j++) {
            sumd += x_s[j] * pow(d[j], i);
        }
        zeta[i] = PI / 6.* sumd;
    }
    double den=eta/zeta[3]/summ;
    return den;
}
/*from mol density calculate packing fraction*/
double rho_to_eta(double t, double rho, vector<double> x, add_args &cppargs){

    int ncomp = x.size(); // number of components
    vector<double> d(ncomp);
    for (int i = 0; i < ncomp; i++) {
        double c=cppargs.lambda_r[i]/(cppargs.lambda_r[i]-cppargs.lambda_a[i])*pow(cppargs.lambda_r[i]/cppargs.lambda_a[i],cppargs.lambda_a[i]/(cppargs.lambda_r[i]-cppargs.lambda_a[i]));
        double theta=c*cppargs.e[i]/t;
        d[i]=saft_gl10_int(theta,1.,cppargs.e[i],cppargs.s[i],t,cppargs.lambda_r[i],cppargs.lambda_a[i]);
    }
    if (!cppargs.z.empty()) {
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z[i] != 0) {
                d[i] = cppargs.s[i] * (1 - 0.12); // for ions the diameter is assumed to be temperature independent (see Held et al. 2014)
            }
        }
    }
    //caculate the mole fraction of segments of component i. A8
    vector<double> x_s(ncomp);
    double summ = 0.;
    for (int i = 0; i < ncomp; i++) {
        summ += cppargs.m[i] * x[i];
    }
    for (int i = 0; i < ncomp; i++) {
        x_s[i] = cppargs.m[i] * x[i] / summ;
    }

    //caculate the total number density of spherical segments.
    double den = rho * N_AV / 1.0e30;
    double den_s = den * summ;
    //caculate the moments of the number density zeta. A7
    vector<double> zeta(4, 0);
    for (int i = 0; i < 4; i++) {
        double sumd = 0;
        for (int j = 0; j < ncomp; j++) {
            sumd += x_s[j] * pow(d[j], i);
        }
        zeta[i] = PI / 6. * den_s * sumd;
    }
    double eta = zeta[3];
    return eta;
}
//Ideal part: Coniglio[1997]. A group contribution method for estimating ideal gas heat capacitie of hydrocarbons
// doi:10.1016/S0378-3812(97)00143-X
double Plank_Einstein(double x){
    return x*x*exp(x)/pow(exp(x)-1.,2.);
}
double ideal_capacity(double t,vector<double> para_A,int is_line){
    double result=0;
    for(int i=1;i<6;++i){
        result+=N_AV*kb*para_A[i]*Plank_Einstein(i*para_A[0]/t);
    }
    if(is_line==LINEAR){
        result+=3.5*N_AV*kb;
    }
    else if (is_line==NON_LINEAR)
    {
        result+=4*N_AV*kb;
    }
    else{
        result=-1;
    }
    return result;
}

double ideal_capacity(double t,vector<vector<double>> para_A,vector<int> group_num,int is_line){
    double result=0;
    for(int j=0;j<para_A.size();++j){
        for(int i=1;i<6;++i){
            result+=group_num[j]*N_AV*kb*para_A[j][i]*Plank_Einstein(i*para_A[j][0]/t);
        }
    }
    if(is_line==LINEAR){
        result+=3.5*N_AV*kb;
    }
    else if (is_line==NON_LINEAR)
    {
        result+=4*N_AV*kb;
    }
    else{
        result=-1;
    }
    return result;
}

double ideal_enthalpy(double t,vector<double> para_A,int is_line){
    double result=0;
    for(int i=1;i<6;++i){
        result+=N_AV*kb*para_A[i]*i*para_A[0]/(exp(i*para_A[0]/t)-1);
    }
    if(is_line==LINEAR){
        result+=3.5*N_AV*kb*t;
    }
    else if (is_line==NON_LINEAR)
    {
        result+=4*N_AV*kb*t;
    }
    else{
        result=-1;
    }
    return result;
}

double ideal_enthalpy(double t,vector<vector<double>> para_A,vector<int> group_num,int is_line){
    double result=0;
    for(int j=0;j<para_A.size();++j){
        for(int i=1;i<6;++i){
            result+=group_num[j]*N_AV*kb*para_A[j][i]*i*para_A[j][0]/(exp(i*para_A[j][0]/t)-1);
        }
    }
    if(is_line==LINEAR){
        result+=3.5*N_AV*kb*t;
    }
    else if (is_line==NON_LINEAR)
    {
        result+=4*N_AV*kb*t;
    }
    else{
        result=-1;
    }
    return result;
}

/*Calculate by Plank-Einstein function*/
double ideal_entropy(double t,double p,vector<double> para_A,int is_line){
    double result=0;
    for(int i=1;i<6;++i){
        double temp=0;
        temp=i*para_A[0]/t*(1./(exp(i*para_A[0]/t)-1.)+1.)-log(exp(i*para_A[0]/t)-1.);
        result+=N_AV*kb*para_A[i]*temp;
    }
    if(is_line==LINEAR){
        result+=3.5*N_AV*kb*log(t)-N_AV*kb*log(p);
    }
    else if (is_line==NON_LINEAR)
    {
        result+=4*N_AV*kb*log(t)-N_AV*kb*log(p);
    }
    else{
        result=-1;
    }
    return result;
}

double ideal_entropy(double t,double p,vector<vector<double>> para_A,vector<int> group_num,int is_line){
    double result=0;
    for(int j=0;j<para_A.size();++j){
        for(int i=1;i<6;++i){
            double temp=0;
            temp=i*para_A[j][0]/t*(1./(exp(i*para_A[j][0]/t)-1.)+1.)-log(exp(i*para_A[j][0]/t)-1.);
            result+=group_num[j]*N_AV*kb*para_A[j][i]*temp;
        }
    }
    if(is_line==LINEAR){
        result+=3.5*N_AV*kb*log(t)-N_AV*kb*log(p);
    }
    else if (is_line==NON_LINEAR)
    {
        result+=4*N_AV*kb*log(t)-N_AV*kb*log(p);
    }
    else{
        result=-1;
    }
    return result;
}
/*
this section used residual properties and ideal gas properties, the ideal part was calculated 
    by Planck-Einstein function fitted to experimental data. 
Residual part: Lafitte[2006]. Simultaneous estimation of phase behavior and second-derivative properties using the 
statistical associating fluid theory with variable range approach http://dx.doi.org/10.1063/1.2140276
Ideal part: Coniglio[1997]. A group contribution method for estimating ideal gas heat capacitie of hydrocarbons doi:10.1016/S0378-3812(97)00143-X
*/
vector<double> saft_heatcapacity_cpp(double t, double p, vector<double> x, vector<double>paras_A,int is_line, add_args &cppargs){

    double rho=saft_density(t,p,x,cppargs);
    double eta=rho_to_eta(t,rho,x,cppargs);
    
    double dpdeta=dpdeta_cpp(t,eta,x,cppargs);
    double dpdt=dpdt_cpp(t,rho,x,cppargs);
    double beta=1./eta/dpdeta;
    double alpha=beta*dpdt;

    double dadt=daresdt_cpp(t,rho,x,cppargs);
    double da2dt2=dares2dt2_cpp(t,rho,x,cppargs);
    double Cv_res=-2*N_AV*kb*t*dadt-N_AV*kb*t*t*da2dt2;
    double Cp_res=Cv_res-N_AV*kb+alpha*alpha/beta*t/rho;

    double Cp_ideal=ideal_capacity(t,paras_A,is_line);
    double Cv_ideal=Cp_ideal-N_AV*kb;
    vector<double> result;
    result.push_back(Cp_res+Cp_ideal);
    result.push_back(Cv_res+Cv_ideal);
    return result;
}
vector<double> saft_heatcapacity_cpp(double t, double p, vector<double> x,int is_line, add_args &cppargs){

    double rho=saft_density(t,p,x,cppargs);
    double eta=rho_to_eta(t,rho,x,cppargs);
    
    double dpdeta=dpdeta_cpp(t,eta,x,cppargs);
    double dpdt=dpdt_cpp(t,rho,x,cppargs);
    double beta=1./eta/dpdeta;
    double alpha=beta*dpdt;

    double dadt=daresdt_cpp(t,rho,x,cppargs);
    double da2dt2=dares2dt2_cpp(t,rho,x,cppargs);
    double Cv_res=-2*N_AV*kb*t*dadt-N_AV*kb*t*t*da2dt2;
    double Cp_res=Cv_res-N_AV*kb+alpha*alpha/beta*t/rho;

    double Cp_ideal=ideal_capacity(t,cppargs.group_A,cppargs.group_num[0],is_line);
    double Cv_ideal=Cp_ideal-N_AV*kb;
    vector<double> result;
    result.push_back(Cp_res+Cp_ideal);
    result.push_back(Cv_res+Cv_ideal);
    return result;
}


//alpha: Isobaric thermal expansivity (K^-1) &&beta: Isothermal compressibility (Pa^-1)
vector<double> saft_SecondIso_cpp(double t, double p, vector<double> x,add_args &cppargs){
    double rho=saft_density(t,p,x,cppargs);
    double eta=rho_to_eta(t,rho,x,cppargs);    
    double dpdeta=dpdeta_cpp(t,eta,x,cppargs);
    double dpdt=dpdt_cpp(t,rho,x,cppargs);
    double beta=1./eta/dpdeta;
    double alpha=beta*dpdt;
    vector<double> result;
    result.push_back(alpha);
    result.push_back(beta);
    return result;
}
//thermodynamic relationship formulas Cp=(dh/dt)_p and so on.
vector<double> saft_heatcapacity_dir_cpp(double t, double p, vector<double> x, vector<double>paras_A,int mole, add_args &cppargs){
    double Cp=saft_dhdtp_cpp(t,p,x,paras_A,mole,cppargs);
    double rho=saft_density(t,p,x,cppargs);
    double Cv=saft_dhdtrho_cpp(t,rho,x,paras_A,mole,cppargs)-dpdt_cpp(t,rho,x,cppargs)/rho;
    vector<double> result;
    result.push_back(Cp);
    result.push_back(Cv);
    return result;
}
//calculate enthalpy J/mol
double saft_enthalpy_cpp(double t, double p, vector<double> x, vector<double>paras_A,int mole, add_args &cppargs){
    double rho=saft_density(t,p,x,cppargs);
    double h_res=saft_hres_cpp(t,rho,x,cppargs);
    double h_ideal=ideal_enthalpy(t,paras_A,mole);
    return h_ideal+h_res;
}

double saft_enthalpy_cpp(double t, double p, vector<double> x,int is_line, add_args &cppargs){
    double rho=saft_density(t,p,x,cppargs);
    double h_res=saft_hres_cpp(t,rho,x,cppargs);
    double h_ideal=ideal_enthalpy(t,cppargs.group_A,cppargs.group_num[0],is_line);
    return h_ideal+h_res;
}

//calculate derivative (dh/dt)_p
double saft_dhdtp_cpp(double t, double p, vector<double> x, vector<double>paras_A,int mole, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=t*0.01;
    double fr=saft_enthalpy_cpp(t+hh,p,x,paras_A,mole,cppargs);
    double fl=saft_enthalpy_cpp(t-hh,p,x,paras_A,mole,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_enthalpy_cpp(t+hh,p,x,paras_A,mole,cppargs);
        fl=saft_enthalpy_cpp(t-hh,p,x,paras_A,mole,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
//calculate derivative (dh/dt)_p
double saft_enthalpybyrho_cpp(double t, double rho, vector<double> x, vector<double>paras_A,int mole, add_args &cppargs){
    double p=saft_p_cpp(t,rho,x,cppargs);
    return saft_enthalpy_cpp(t,p,x,paras_A,mole,cppargs);
}
//calculate derivative (dh/dt)_rho
double saft_dhdtrho_cpp(double t, double rho, vector<double> x, vector<double>paras_A,int mole, add_args &cppargs){
    double errt,fac,hh,ans;
    vector<double> a;
    hh=t*0.01;
    double fr=saft_enthalpybyrho_cpp(t+hh,rho,x,paras_A,mole,cppargs);
    double fl=saft_enthalpybyrho_cpp(t-hh,rho,x,paras_A,mole,cppargs);
    a.push_back((fr-fl)/2/hh);
    double err=big;
    for(int i=1;i<10;++i){
        hh/=con;
        fr=saft_enthalpybyrho_cpp(t+hh,rho,x,paras_A,mole,cppargs);
        fl=saft_enthalpybyrho_cpp(t-hh,rho,x,paras_A,mole,cppargs);
        a.push_back((fr-fl)/2/hh);
        int k=i*(i+1)/2;
        fac=con2;
        for(int j=0;j<i;++j){
            a.push_back((fac*a[k+j]-a[k+j-i])/(fac-1));
            fac=con2*fac;
            errt=std::max(fabs(a[k+j+1]-a[k+j]),fabs(a[k+j+1]-a[k+j-i]));
            if(errt<=err){
                err=errt;
                ans=a[k+j+1];
            }
        }
        if(i>1){
            if(fabs(a[i*(i+1)/2-1]-a[(i-2)*(i+1)/2])>=safe*err){
                break;
            }
        }
    }
    return ans;
}
//calculate enthalpy J/(mol*k)
double saft_entropy_cpp(double t, double p, vector<double> x, vector<double>paras_A,int mole, add_args &cppargs){
    double rho=saft_density(t,p,x,cppargs);
    double s_res=saft_sres_cpp(t,rho,x,cppargs);
    double s_ideal=ideal_entropy(t,p,paras_A,mole);
    return s_ideal+s_res;
}

template <typename T>
T da1ddens_cpp(T t, T den_s, const vector<double>& x_s, const add_args& cppargs, int num, const vector<T>& d) {
    double den_val = static_cast<double>(den_s);
    double hh = den_val * 0.01;
    // ⚠️ 注意：请确保下方调用的 saft_a1_tmp 参数顺序与你原代码完全一致
    T fr = saft_a1_tmp<T>(t, den_s + hh, x_s, cppargs, num, d);
    T fl = saft_a1_tmp<T>(t, den_s - hh, x_s, cppargs, num, d);

    vector<T> a;
    a.push_back((fr - fl) / 2.0 / hh);
    
    T err = 1e30; // 替换原来的 big 常量
    T ans = a[0];
    
    double con = 1.4; 
    double con2 = con * con;
    double safe = 2.0;

    for (int i = 1; i < 10; ++i) {
        hh /= con;
        fr = saft_a1_tmp<T>(t, den_s + hh, x_s, cppargs, num, d);
        fl = saft_a1_tmp<T>(t, den_s - hh, x_s, cppargs, num, d);
        a.push_back((fr - fl) / 2.0 / hh);

        int k = i * (i + 1) / 2;
        double fac = con2;
        for (int j = 0; j < i; ++j) {
            a.push_back((fac * a[k + j] - a[k + j - i]) / (fac - 1.0));
            fac = con2 * fac;
            
            // 【核心修正】：用 autodiff 兼容的 abs 替换 std::fabs，用三元运算符替换 std::max
            T err1 = abs(a[k + j + 1] - a[k + j]);
            T err2 = abs(a[k + j + 1] - a[k + j - i]);
            T errt = (err1 > err2) ? err1 : err2;
            
            if (errt <= err) {
                err = errt;
                ans = a[k + j + 1];
            }
        }
        T diff = abs(a[i * (i + 1) / 2 - 1] - a[(i - 2) * (i + 1) / 2]);
        if (i > 1 && diff >= safe * err) {
            break;
        }
    }
    return ans;
}

template <typename T>
T da2ddens_cpp(T t, T den_s, const vector<double>& x_s, const add_args& cppargs, int num, const vector<T>& d) {
    double den_val = static_cast<double>(den_s);
    double hh = den_val * 0.01;
    T fr = saft_a2_tmp<T>(t, den_s + hh, x_s, cppargs, d); 
    T fl = saft_a2_tmp<T>(t, den_s - hh, x_s, cppargs, d);

    // chi 保持中心点求值
    T chi_r = saft_chi_tmp<T>(t, den_s + hh, x_s, cppargs, num, d);
    T chi_l = saft_chi_tmp<T>(t, den_s - hh, x_s, cppargs, num, d);

    fr = fr / (1.0 + chi_r);
    fl = fl / (1.0 + chi_l);

    vector<T> a;
    a.push_back((fr - fl) / 2.0 / hh);
    
    T err = 1e30; 
    T ans = a[0];

    double con = 1.4; 
    double con2 = con * con;
    double safe = 2.0;

    for (int i = 1; i < 10; ++i) {
        hh /= con;
        fr = saft_a2_tmp<T>(t, den_s + hh, x_s, cppargs, d);
        fl = saft_a2_tmp<T>(t, den_s - hh, x_s, cppargs, d);
        
        chi_r = saft_chi_tmp<T>(t, den_s + hh, x_s, cppargs, num, d);
        chi_l = saft_chi_tmp<T>(t, den_s - hh, x_s, cppargs, num, d);
        
        fr = fr / (1.0 + chi_r);
        fl = fl / (1.0 + chi_l);
        
        a.push_back((fr - fl) / 2.0 / hh);

        int k = i * (i + 1) / 2;
        double fac = con2;
        for (int j = 0; j < i; ++j) {
            a.push_back((fac * a[k + j] - a[k + j - i]) / (fac - 1.0));
            fac = con2 * fac;
            
            T err1 = abs(a[k + j + 1] - a[k + j]);
            T err2 = abs(a[k + j + 1] - a[k + j - i]);
            T errt = (err1 > err2) ? err1 : err2;
            
            if (errt <= err) {
                err = errt;
                ans = a[k + j + 1];
            }
        }
        T diff = abs(a[i * (i + 1) / 2 - 1] - a[(i - 2) * (i + 1) / 2]);
        if (i > 1 && diff >= safe * err) {
            break;
        }
    }
    return ans;
}


// 通用中心差分理查森外推求导器
double calc_derivative(std::function<double(double)> target_func, double x, double hh) {
    std::vector<double> a; // 存储不同阶数的导数估计
    double fr = target_func(x + hh);
    double fl = target_func(x - hh);
    a.push_back((fr - fl) / (2.0 * hh));
    
    double err = big;
    double ans = a[0];
    
    for (int i = 1; i < 10; ++i) {
        hh /= con;
        fr = target_func(x + hh);
        fl = target_func(x - hh);
        a.push_back((fr - fl) / (2.0 * hh));
        
        int k = i * (i + 1) / 2;
        double fac = con2;
        for (int j = 0; j < i; ++j) {
            a.push_back((fac * a[k + j] - a[k + j - i]) / (fac - 1.0));
            fac *= con2;
            double errt = std::max(std::fabs(a[k + j + 1] - a[k + j]), std::fabs(a[k + j + 1] - a[k + j - i]));
            if (errt <= err) {
                err = errt;
                ans = a[k + j + 1];
            }
        }
        if (i > 1 && std::fabs(a[i * (i + 1) / 2 - 1] - a[(i - 2) * (i + 1) / 2]) >= safe * err) break;
    }
    return ans;
}

// 实例化 saft_ares_cpp
template double saft_ares_cpp<double>(double, double, const vector<double>&, const add_args&);
template autodiff::dual saft_ares_cpp<autodiff::dual>(autodiff::dual, autodiff::dual, const vector<double>&, const add_args&);
template autodiff::dual2nd saft_ares_cpp<autodiff::dual2nd>(autodiff::dual2nd, autodiff::dual2nd, const vector<double>&, const add_args&);

// 实例化 saft_B_lambda
template double saft_B_lambda<double>(double, double, double, double, double, double);
template autodiff::dual saft_B_lambda<autodiff::dual>(autodiff::dual, autodiff::dual, double, autodiff::dual, double, autodiff::dual);
template autodiff::dual2nd saft_B_lambda<autodiff::dual2nd>(autodiff::dual2nd, autodiff::dual2nd, double, autodiff::dual2nd, double, autodiff::dual2nd);

// 实例化 saft_as_lambda
template double saft_as_lambda<double>(double, double, double, double, double);
template autodiff::dual saft_as_lambda<autodiff::dual>(autodiff::dual, autodiff::dual, double, autodiff::dual, double);
template autodiff::dual2nd saft_as_lambda<autodiff::dual2nd>(autodiff::dual2nd, autodiff::dual2nd, double, autodiff::dual2nd, double);

// 实例化 saft_I_lambda
template double saft_I_lambda<double>(double, double);
template autodiff::dual saft_I_lambda<autodiff::dual>(double, autodiff::dual);
template autodiff::dual2nd saft_I_lambda<autodiff::dual2nd>(double, autodiff::dual2nd);

// 实例化 saft_J_lambda
template double saft_J_lambda<double>(double, double);
template autodiff::dual saft_J_lambda<autodiff::dual>(double, autodiff::dual);
template autodiff::dual2nd saft_J_lambda<autodiff::dual2nd>(double, autodiff::dual2nd);

// 实例化 saft_zeff_lambda
template double saft_zeff_lambda<double>(double, double);
template autodiff::dual saft_zeff_lambda<autodiff::dual>(autodiff::dual, double);
template autodiff::dual2nd saft_zeff_lambda<autodiff::dual2nd>(autodiff::dual2nd, double);

// 实例化 saft_a1_tmp
template double saft_a1_tmp<double>(double, double, const vector<double>&, const add_args&, int, const vector<double>&);
template autodiff::dual saft_a1_tmp<autodiff::dual>(autodiff::dual, autodiff::dual, const vector<double>&, const add_args&, int, const vector<autodiff::dual>&);
template autodiff::dual2nd saft_a1_tmp<autodiff::dual2nd>(autodiff::dual2nd, autodiff::dual2nd, const vector<double>&, const add_args&, int, const vector<autodiff::dual2nd>&);

// 实例化 saft_a2_tmp
template double saft_a2_tmp<double>(double, double, const vector<double>&, const add_args&, const vector<double>&);
template autodiff::dual saft_a2_tmp<autodiff::dual>(autodiff::dual, autodiff::dual, const vector<double>&, const add_args&, const vector<autodiff::dual>&);
template autodiff::dual2nd saft_a2_tmp<autodiff::dual2nd>(autodiff::dual2nd, autodiff::dual2nd, const vector<double>&, const add_args&, const vector<autodiff::dual2nd>&);

// 实例化 saft_chi_tmp
template double saft_chi_tmp<double>(double, double, const vector<double>&, const add_args&, int, const vector<double>&);
template autodiff::dual saft_chi_tmp<autodiff::dual>(autodiff::dual, autodiff::dual, const vector<double>&, const add_args&, int, const vector<autodiff::dual>&);
template autodiff::dual2nd saft_chi_tmp<autodiff::dual2nd>(autodiff::dual2nd, autodiff::dual2nd, const vector<double>&, const add_args&, int, const vector<autodiff::dual2nd>&);

// 实例化 saft_func_glf
template double saft_func_glf<double>(double, double, double, double, double);
template autodiff::dual saft_func_glf<autodiff::dual>(double, autodiff::dual, autodiff::dual, double, double);
template autodiff::dual2nd saft_func_glf<autodiff::dual2nd>(double, autodiff::dual2nd, autodiff::dual2nd, double, double);