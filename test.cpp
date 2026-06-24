#pragma GCC optimize(3)
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include "saft_mie.h"

using namespace std;

constexpr double mole_mass = 44.0095;
using namespace std;

int main() {
    add_args CO2;
    CO2.m.push_back(1.41222);
    CO2.s.push_back(3.22324);
    CO2.e.push_back(192.341);
    CO2.lambda_a.push_back(5.99299);
    CO2.lambda_r.push_back(11.0216);
    
    // 如果是优化版，确保调用了交叉项初始化！
    init_saft_cross_params(CO2); 

    double para_A[6] = {966.45, 2.0502, 0.72817, 1.1321, -0.29169, 0.53293};
    vector<double> A(para_A, para_A + 6);
    int mole = 1;
    // ==========================================
    // 2. 设定固定工况 (已知温度和真实密度)
    // ==========================================
    double t = 273.15; // 0 °C
    double rho = 21358.0; // 约 940 kg/m^3 的液态 CO2 密度
    vector<double> x = {1.0}; // 纯组分

    cout << "==========================================" << endl;
    cout << "    SAFT-VR Mie a_res 纯净前向计算测试" << endl;
    cout << "==========================================" << endl;
    cout << "设定工况:" << endl;
    cout << "温度 T   = " << t << " K" << endl;
    cout << "密度 rho = " << rho << " mol/m^3" << endl;

    // ==========================================
    // 3. 计算单点 a_res
    // ==========================================
    // 注意：只算值的话，直接传入 <double> 模板参数即可
    double ares_val = saft_ares_cpp<double>(t, rho, x, CO2);

    cout << fixed << setprecision(6);
    cout << "------------------------------------------" << endl;
    cout << "计算结果:" << endl;
    cout << "a_res = " << ares_val << endl;
    cout << "------------------------------------------" << endl;

    return 0;
}