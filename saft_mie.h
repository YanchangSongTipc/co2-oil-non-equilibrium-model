#pragma once // 防止头文件被重复包含 (如果你用的是 #ifndef ... 也可以)

#include <vector>
#include <array>
#include <string>

#include<limits>
#include <functional>
// 宏定义可以放在这里
#if defined(HUGE_VAL) && !defined(_HUGE)
# define _HUGE HUGE_VAL
#else
// GCC Version of huge value macro
#if defined(HUGE) && !defined(_HUGE)
#  define _HUGE HUGE
#endif
#endif
//Constants

const static double kb = 1.380648465952442093e-23; // Boltzmann constant, J K^-1
const static double PI = 3.141592653589793;
const static double N_AV = 6.022140857e23; // Avagadro's number
const static double E_CHRG = 1.6021766208e-19; // elementary charge, units of coulomb
const static double perm_vac = 8.854187817e-22; //permittivity in vacuum, C V^-1 Angstrom^-1#pragma once
const static double differential_coefficient = 1e-5;
const static double h_planck = 6.62607015e-34;//Planck constant, J*s
const static int LIQUID = 0;
const static int VAPER = 1;
const static int LINEAR = 1;
const static int NON_LINEAR = 0;
const static double big = std::numeric_limits<double>::max();
const static int ntab=10;
const static double con=1.4,con2=(con*con);
const static double safe=2.0;

#ifndef DBL_EPSILON
    #define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif

struct add_args {
    std::vector<double> m;
    std::vector<double> s;
    std::vector<double> e;
    std::vector<double> k_ij;
    std::vector<double> e_assoc;
    std::vector<double> rc;
    std::vector<double> dipm;
    std::vector<double> dip_num;
    std::vector<double> z;
    double dielc;
    std::vector<int> assoc_num;
    std::vector<int> assoc_scheme;
    std::vector<double> k_hb;
    std::vector<double> l_ij;
    std::vector<double> mol_mass;
    std::vector<double> lambda_a;
    std::vector<double> lambda_r;
    std::vector<std::vector<double>> group_A;
    std::vector<std::vector<int>> group_num;
    std::vector<double> e_ij;
    std::vector<double> s_ij;
    std::vector<double> c_ij;
    std::vector<double> lambda_aij;
    std::vector<double> lambda_rij;
};


// 最核心的 ares 函数声明
template <typename T>
T saft_ares_cpp(T t, T rho, const std::vector<double>& x, const add_args& cppargs);
double saft_A_cpp(double t, double rho, std::vector<double> x, add_args& cppargs);
double saft_Z_cpp(double t, double rho, const std::vector<double>& x, const add_args& cppargs);
double saft_p_cpp(double t, double rho, const std::vector<double>& x, const add_args& cppargs);
double saft_pv_cpp(double t,double v,std::vector<double> x, add_args &cppargs);
double saft_newp_cpp(double t, double rho, std::vector<double> x, add_args& cppargs);
double saft_t_cpp(double p,double rho,std::vector<double> x,add_args cppargs);
double saft_gres_cpp(double t, double rho, std::vector<double> x, add_args& cppargs);
double saft_hres_cpp(double t, double rho, std::vector<double> x, add_args &cppargs);
double saft_sres_cpp(double t, double rho, std::vector<double> x, add_args &cppargs);
double saft_rho_cpp(double t, double p, std::vector<double> x, int phase, add_args& cppargs);
double dZdetat_cpp(double t, double eta, std::vector<double> x, add_args& cppargs);
double dZdtrho_cpp(double t, double rho, std::vector<double> x, add_args& cppargs);

void init_saft_cross_params(add_args& cppargs);


double daresdt_cpp(double t,double rho,std::vector<double> x, add_args &cppargs);
std::vector<double> tmp_daresdt_cpp(double t,double rho,std::vector<double> x, add_args &cppargs);
double dares2dt2_cpp(double t,double rho,std::vector<double> x, add_args &cppargs);
double daresdrho_cpp(double t,double rho,std::vector<double> x, add_args &cppargs);
double dZdt_cpp(double t,double rho,std::vector<double> x, add_args &cppargs);
double dpdt_cpp(double t,double rho,std::vector<double> x, add_args &cppargs);
double dpdeta_cpp(double t,double eta,std::vector<double> x,add_args& cppargs);
double dpdv_cpp(double t,double v,std::vector<double> x, add_args &cppargs);
double dp2dv2_cpp(double t,double v,std::vector<double> x, add_args &cppargs);
double dpdrho_cpp(double t,double rho,std::vector<double> x, add_args &cppargs);
double dp2drho2_cpp(double t,double rho,std::vector<double> x, add_args &cppargs);
double dpfuncdxi(double(*function)(double x1,double x2,std::vector<double> x, add_args &cppargs),
    std::vector<double> xi,std::vector<double> x,add_args &cppargs,int n);
double pfunction(double(*function)(double t,double rho,std::vector<double> x, add_args &cppargs),
    std::vector<double> xi,std::vector<double> x,add_args &cppargs);
double saft_density(double t,double p,std::vector<double> x,add_args &cppargs);
double saft_drhodp(double t,double p,std::vector<double> x,add_args &cppargs);
std::vector<double> saft_SecondIso_cpp(double t, double p, std::vector<double> x,add_args &cppargs);
std::vector<double> saft_heatcapacity_cpp(double t, double p, std::vector<double> x, std::vector<double>paras_A,int mole, add_args &cppargs);
std::vector<double> saft_heatcapacity_cpp(double t, double p, std::vector<double> x,int is_line, add_args &cppargs);
std::vector<double> saft_heatcapacity_dir_cpp(double t, double p, std::vector<double> x, std::vector<double>paras_A,int mole, add_args &cppargs);
double saft_enthalpy_cpp(double t, double p, std::vector<double> x, std::vector<double>paras_A,int mole, add_args &cppargs);
double saft_enthalpy_cpp(double t, double p, std::vector<double> x,int is_line, add_args &cppargs);
double saft_entropy_cpp(double t, double p, std::vector<double> x, std::vector<double>paras_A,int mole, add_args &cppargs);

double saft_dhdtp_cpp(double t, double p, std::vector<double> x, std::vector<double>paras_A,int mole, add_args &cppargs);

double saft_enthalpybyrho_cpp(double t, double rho, std::vector<double> x, std::vector<double>paras_A,int mole, add_args &cppargs);
double saft_dhdtrho_cpp(double t, double rho, std::vector<double> x, std::vector<double>paras_A,int mole, add_args &cppargs);

std::vector<double> saft_fugcoe_cpp(double t, double rho, std::vector<double> x, add_args& cppargs);
double saft_func_lj(double e, double s, double t, double r, double lambda_r, double lambda_a);
double saft_lj_mie(double e, double s, double r, double n, double m);
double saft_func_dljdt(double e, double s, double t, double r, double lambda_r, double lambda_a);
template <typename T>
T saft_func_glf(double e, T t, T x, double lambda_r, double lambda_a);

//numerical integration for the integral of hard-sphere system diameter d
double saft_rbg_int(double a, double b, double err, double e, double s, double t, double lambda_r, double lambda_a);
double saft_rbg_dintdt(double a, double b, double err, double e, double s, double t, double lambda_r, double lambda_a);
template <typename T>
T saft_gl5_int(T r, double a, double e, double s, T t, double lambda_r, double lambda_a);
template <typename T>
T saft_gl10_int(T r, double a, double e, double s, T t, double lambda_r, double lambda_a);



//assistant functions for the integral of hard-sphere system diameter d
template <typename T>
T saft_I_lambda(double lambda, T x);
template <typename T>
T saft_J_lambda(double lambda, T x);
template <typename T>
T saft_as_lambda(T den_s, T d_ij, double e_ij, T zeta_x, double lambda);
template <typename T>
T saft_B_lambda(T den_s, T d_ij, double e_ij, T zeta_x, double lambda, T x0);
template <typename T>
T saft_zeff_lambda(T z, double lambda);



template <typename T>
T saft_a1_tmp(T t, T den_s, const std::vector<double>& x_s, const add_args& cppargs, int num, const std::vector<T>& d);
template <typename T>
T saft_a2_tmp(T t, T den_s, const std::vector<double>& x_s, const add_args& cppargs, const std::vector<T>& d);
template <typename T>
T saft_chi_tmp(T t, T den_s, const std::vector<double>& x_s, const add_args& cppargs, int num, const std::vector<T>& d);

template <typename T>
T da1ddens_cpp(T t, T den_s, const std::vector<double>& x_s, const add_args& cppargs, int num, const std::vector<T>& d);
template <typename T>
T da2ddens_cpp(T t, T den_s, const std::vector<double>& x_s, const add_args& cppargs, int num, const std::vector<T>& d);


double saft_fk_alpha(int k, std::vector<std::vector<double>> phi, double aplha);

double resid_rho(double rhomolar, double t, double p, std::vector<double> x, add_args& cppargs);
double reduced_to_molar(double nu, double t, int ncomp, std::vector<double> x, add_args& cppargs,std::vector<double> d);
double rho_to_eta(double t, double rho, std::vector<double> x, add_args &cppargs);
double eta_to_den(double t, double eta, std::vector<double> x, add_args &cppargs);
double from_eta_to_rho(double t, double rho, std::vector<double> x_s, double summ, add_args& cppargs,std::vector<double> d);

std::vector<double> flashTQ_cpp(double t, double Q, std::vector<double> x, add_args &cppargs);
std::vector<double> outerTQ(double p_guess, double t, double Q, std::vector<double> x, add_args &cppargs);
double estimate_flash_p(double t, double Q, std::vector<double> x, add_args &cppargs);
double resid_inner(double R, double kb0, double Q, std::vector<double> u, std::vector<double> x, add_args &cppargs);
double BoundedSecantInner(double kb0, double Q, std::vector<double> u, std::vector<double> x, add_args &cppargs, double x0, double xmin,
    double xmax, double dx, double tol, int maxiter);
add_args get_single_component(int i, add_args &cppargs);

double BrentRho(double t, double p, std::vector<double> x, int phase, add_args& cppargs, double a, double b,
    double macheps, double tol_abs, int maxiter);
double support_func(double t,double p,double rho,std::vector<double> x,add_args cppargs);
double calc_derivative(std::function<double(double)> target_func, double x, double hh);
class ValueError : public std::exception
{
public:
    ValueError(const std::string& err) throw() : m_err(err) {}
    ~ValueError() throw() {};
    virtual const char* what() const throw() { return m_err.c_str(); }
private:
    std::string m_err;
};

class SolutionError : public std::exception
{
public:
    SolutionError(const std::string& err) throw() : m_err(err) {}
    ~SolutionError() throw() {};
    virtual const char* what() const throw() { return m_err.c_str(); }
private:
    std::string m_err;
};

//functions used for heat capacity calculations
double Plank_Einstein(double x);
double ideal_capacity(double t,std::vector<double> para_A,int is_line);
double ideal_capacity(double t,std::vector<std::vector<double>> para_A,std::vector<int> group_num,int is_line);
double ideal_enthalpy(double t,std::vector<double> para_A,int is_line);
double ideal_enthalpy(double t,std::vector<std::vector<double>> para_A,std::vector<int> group_num,int is_line);
double ideal_entropy(double t,double p,std::vector<double> para_A,int is_line);
double ideal_entropy(double t,double p,std::vector<std::vector<double>> para_A,std::vector<int> group_num,int is_line);

const std::vector<double> phi1{ 7.5365557,-37.60463,71.745953,-46.83552,-2.467982,-0.50272,8.0956883 };
const std::vector<double> phi2{ -359.44,1825.6,-3168.0,1884.2,-0.82376,-3.1935,3.7090 };
const std::vector<double> phi3{ 1550.9,-5070.1,6534.6,-3288.7,-2.7171,2.0883,0. };
const std::vector<double> phi4{ -1.19932,9.063632,-17.9482,11.34027,20.52142,-56.6377,40.53683 };
const std::vector<double> phi5{ -1911.28,21390.175,-51320.7,37064.54,1103.742,-3264.61,2556.181 };
const std::vector<double> phi6{ 9236.9,-129430,357230,-315530,1390.2,-4518.2,4241.6 };
const std::vector<std::vector<double>> phi_k = { phi1,phi2,phi3,phi4,phi5,phi6 };