
#include <vector>
#include <cmath>
#include "rates.h"


enum StateIndices {
    
    Su_h, Eu_h, Iu_h, Ru_h, Nu_h, Su_m, Eu_m, Iu_m, Nu_m,
    Sp_h, Ep_h, Ip_h, Rp_h, Np_h, Sp_p, Ep_p, Ip_p, Rp_p, Np_p, Sp_m, Ep_m, Ip_m, Np_m,
    Ss_p, Es_p, Is_p, Rs_p, Ns_p, Ss_m, Es_m, Is_m, Ns_m,
    U_infection_QOI, P_human_infection_QOI, P_NHP_infection_QOI, S_NHP_infection_QOI
};

enum ParamIndices {
    mu_p, vu_p,
    Theta_u, Theta_p, Theta_s,
    b_u, b_p, b_s,
    beta_hm_u, beta_hm_p,
    beta_pm_p, beta_pm_s,
    beta_mh_u, beta_mh_p,
    beta_mp_p, beta_mp_s,
    epsilon_h, epsilon_p,
    epsilon_u_m, epsilon_p_m, epsilon_s_m,
    rho_h, rho_p,
    eta_p,
    tau_u, tau_p,
    varphi_p, varphi_s,
};

//used for carrying capacity
const int total_urban_host_pop = 7800;
const int total_peri_urban_host_pop = 300 + 202;
const int total_sylvatic_host_pop = 358;

const int K_u = 2*total_urban_host_pop;
const int K_p = 2*total_peri_urban_host_pop;
const int K_s = 2*total_sylvatic_host_pop;


double rate_0(const std::vector<int>& state, const std::vector<double>& p) {
    return p[b_u]*p[beta_mh_u]*state[Su_h]*state[Iu_m]/state[Nu_h];
}


double rate_1(const std::vector<int>& state, const std::vector<double>& p) {
    return p[epsilon_h]*state[Eu_h];
}


double rate_2(const std::vector<int>& state, const std::vector<double>& p) {
    return p[rho_h]*state[Iu_h];
}

double rate_3(const std::vector<int>& state, const std::vector<double>& p) {
    
    return p[Theta_u]*state[Nu_m];
}

double rate_4(const std::vector<int>& state, const std::vector<double>& p) {
    return p[b_u]*p[beta_hm_u]*state[Su_m]*state[Iu_h]/state[Nu_h];
}

double rate_5(const std::vector<int>& state, const std::vector<double>& p) {
    return p[epsilon_u_m]*state[Eu_m];
}

double rate_6(const std::vector<int>& state, const std::vector<double>& p) {
    
    return (p[Theta_u]/K_u)*state[Nu_m]*state[Su_m];
}

double rate_7(const std::vector<int>& state, const std::vector<double>& p) {
    return (p[Theta_u]/K_u)*state[Nu_m]*state[Eu_m];
}

double rate_8(const std::vector<int>& state, const std::vector<double>& p) {
    return  (p[Theta_u]/K_u)*state[Nu_m]*state[Iu_m];    
}

double rate_9(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[b_p]*p[beta_mh_p]*state[Sp_h]*state[Ip_m]/(state[Np_h] + state[Np_p]);
}

double rate_10(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[epsilon_h]*state[Ep_h];
}

double rate_11(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[rho_h]*state[Ip_h];
}

double rate_12(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[Theta_p]*state[Np_m];
}

double rate_13(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[b_p]*p[beta_hm_p]*state[Sp_m]*state[Ip_h]/state[Np_h + Np_p] + p[b_p]*p[beta_pm_p]*state[Sp_m]*state[Ip_p]/(state[Np_p] + state[Np_h]);
}

double rate_14(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[epsilon_p_m]*state[Ep_m];
}


double rate_15(const std::vector<int>& state, const std::vector<double>& p)
{
    return (p[Theta_p]/K_p)*state[Np_m]*state[Sp_m];
}


double rate_16(const std::vector<int>& state, const std::vector<double>& p)
{
    return (p[Theta_p]/K_p)*state[Np_m]*state[Ep_m];
}


double rate_17(const std::vector<int>& state, const std::vector<double>& p)
{
    return (p[Theta_p]/K_p)*state[Np_m]*state[Ip_m];
}


double rate_18(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[vu_p]*state[Np_p];
}


double rate_19(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[b_p]*p[beta_mp_p]*state[Sp_p]*state[Ip_m]/(state[Np_p] + state[Np_h]);
}

double rate_20(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[epsilon_p]*state[Ep_p];
}


double rate_21(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[rho_p]*state[Ip_p];
}


double rate_22(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[eta_p]*state[Rp_p];
}


double rate_23(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[mu_p]*state[Sp_p];
}


double rate_24(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[mu_p]*state[Ep_p];
}


double rate_25(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[mu_p]*state[Ip_p];
}


double rate_26(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[mu_p]*state[Rp_p];
}


double rate_27(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[vu_p]*state[Ns_p];
}


double rate_28(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[b_s]*p[beta_mp_s]*state[Ss_p]*state[Is_m]/state[Ns_p];
}


double rate_29(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[epsilon_p]*state[Es_p];
}


double rate_30(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[rho_p]*state[Is_p];
}


double rate_31(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[eta_p]*state[Rs_p];
}


double rate_32(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[mu_p]*state[Ss_p];
}


double rate_33(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[mu_p]*state[Es_p];
}


double rate_34(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[mu_p]*state[Is_p];
}


double rate_35(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[mu_p]*state[Rs_p];
}


double rate_36(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[Theta_s]*state[Ns_m];
}

double rate_37(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[b_p]*p[beta_pm_p]*state[Ss_m]*state[Is_p]/state[Ns_p];
}

double rate_38(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[epsilon_s_m]*state[Es_m];
}

double rate_39(const std::vector<int>& state, const std::vector<double>& p)
{
    return (p[Theta_s]/K_s)*state[Ns_m]*state[Ss_m];
}

double rate_40(const std::vector<int>& state, const std::vector<double>& p)
{
    return (p[Theta_s]/K_s)*state[Ns_m]*state[Es_m];
}


double rate_41(const std::vector<int>& state, const std::vector<double>& p)
{
    return (p[Theta_s]/K_s)*state[Ns_m]*state[Is_m];
}


//travel rates and events
double rate_42(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[tau_u]*state[Su_h];
}

double rate_43(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[tau_u]*state[Eu_h];
}


double rate_44(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[tau_u]*state[Ru_h];
}


double rate_45(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[tau_p]*state[Sp_h];
}

double rate_46(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[tau_p]*state[Ep_h];
}


double rate_47(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[tau_p]*state[Rp_h];
}


double rate_48(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[varphi_p]*state[Sp_p];
}


double rate_49(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[varphi_p]*state[Ep_p];
}


double rate_50(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[varphi_p]*state[Rp_p];
}


double rate_51(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[varphi_s]*state[Ss_p];
}


double rate_52(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[varphi_s]*state[Es_p];
}


double rate_53(const std::vector<int>& state, const std::vector<double>& p)
{
    return p[varphi_s]*state[Rs_p];
}

std::vector<std::function<double(const std::vector<int>&, const std::vector<double>&)>> rates;



void initialize_rates()
{
    rates.push_back(rate_0);
    rates.push_back(rate_1);
    rates.push_back(rate_2);
    rates.push_back(rate_3);
    rates.push_back(rate_4);
    rates.push_back(rate_5);
    rates.push_back(rate_6);
    rates.push_back(rate_7);
    rates.push_back(rate_8);
    rates.push_back(rate_9);
    rates.push_back(rate_10);
    rates.push_back(rate_11);
    rates.push_back(rate_12);
    rates.push_back(rate_13);
    rates.push_back(rate_14);
    rates.push_back(rate_15);
    rates.push_back(rate_16);
    rates.push_back(rate_17);
    rates.push_back(rate_18);
    rates.push_back(rate_19);
    rates.push_back(rate_20);
    rates.push_back(rate_21);
    rates.push_back(rate_22);
    rates.push_back(rate_23);
    rates.push_back(rate_24);
    rates.push_back(rate_25);
    rates.push_back(rate_26);
    rates.push_back(rate_27);
    rates.push_back(rate_28);
    rates.push_back(rate_29);
    rates.push_back(rate_30);
    rates.push_back(rate_31);
    rates.push_back(rate_32);
    rates.push_back(rate_33);
    rates.push_back(rate_34);
    rates.push_back(rate_35);
    rates.push_back(rate_36);
    rates.push_back(rate_37);
    rates.push_back(rate_38);
    rates.push_back(rate_39);
    rates.push_back(rate_40);
    rates.push_back(rate_41);
    rates.push_back(rate_42);
    rates.push_back(rate_43);
    rates.push_back(rate_44);
    rates.push_back(rate_45);
    rates.push_back(rate_46);
    rates.push_back(rate_47);
    rates.push_back(rate_48);
    rates.push_back(rate_49);
    rates.push_back(rate_50);
    rates.push_back(rate_51);
    rates.push_back(rate_52);
    rates.push_back(rate_53);
}


