const real mScaleFactorGks=1.0;
const real mScaleFactorIto=1.0;
const real mScaleFactorGkr=1.0;
const real mScaleFactorGna=1.0;
const real mScaleFactorAch=1e-24;
const real mScaleFactorGNaK=1.0;
const real mScaleFactorGNaCa=1.0;
const real mScaleFactorGCaL=1.0;
const real mScaleFactorGKur=1.0;
const real mScaleFactorGK1=1.0;
const real mScaleFactorAZD=0.0;

const real var_membrane__R = 8314.0;
const real var_membrane__T = 306.15;
const real var_membrane__F = 96487.0;
const real var_membrane__Cm = 0.05;
const real var_sodium_current__F = var_membrane__F;
const real var_sodium_current__R = var_membrane__R;
const real var_sodium_current__T = var_membrane__T;
const real var_sodium_current__Na_c = var_cleft_space_ion_concentrations__Na_c;
const real var_sodium_current__Na_i = var_intracellular_ion_concentrations__Na_i;
const real var_sodium_current__E_Na = ((var_sodium_current__R * var_sodium_current__T) / var_sodium_current__F) * log(var_sodium_current__Na_c / var_sodium_current__Na_i);
const real var_sodium_current__m = var_sodium_current_m_gate__m;
const real var_sodium_current__h2 = var_sodium_current_h2_gate__h2;
const real var_sodium_current__h1 = var_sodium_current_h1_gate__h1;
const real var_sodium_current__V = var_membrane__V;
const real var_sodium_current__P_Na = 0.0018;
const real var_sodium_current__i_Na = mScaleFactorGna*(((var_sodium_current__P_Na * var_sodium_current__m * var_sodium_current__m * var_sodium_current__m * ((0.9 * var_sodium_current__h1) + (0.1 * var_sodium_current__h2)) * var_sodium_current__Na_c * var_sodium_current__V * var_sodium_current__F * var_sodium_current__F) / (var_sodium_current__R * var_sodium_current__T)) * (exp(((var_sodium_current__V - var_sodium_current__E_Na) * var_sodium_current__F) / (var_sodium_current__R * var_sodium_current__T)) - 1.0)) / (exp((var_sodium_current__V * var_sodium_current__F) / (var_sodium_current__R * var_sodium_current__T)) - 1.0);
const real var_membrane__i_Na = var_sodium_current__i_Na;
const real var_L_type_Ca_channel__k_Ca = 0.025;
const real var_L_type_Ca_channel__Ca_d = var_intracellular_ion_concentrations__Ca_d;
const real var_L_type_Ca_channel__f_Ca = var_L_type_Ca_channel__Ca_d / (var_L_type_Ca_channel__Ca_d + var_L_type_Ca_channel__k_Ca);
const real var_L_type_Ca_channel__E_Ca_app = 60.0;
const real var_L_type_Ca_channel__g_Ca_L = 6.75*mScaleFactorGCaL;
const real var_L_type_Ca_channel__d_L = var_L_type_Ca_channel_d_L_gate__d_L;
const real var_L_type_Ca_channel__f_L2 = var_L_type_Ca_channel_f_L2_gate__f_L2;
const real var_L_type_Ca_channel__V = var_membrane__V;
const real var_L_type_Ca_channel__f_L1 = var_L_type_Ca_channel_f_L1_gate__f_L1;
const real var_L_type_Ca_channel__i_Ca_L = var_L_type_Ca_channel__g_Ca_L * var_L_type_Ca_channel__d_L * ((var_L_type_Ca_channel__f_Ca * var_L_type_Ca_channel__f_L1) + ((1.0 - var_L_type_Ca_channel__f_Ca) * var_L_type_Ca_channel__f_L2)) * (var_L_type_Ca_channel__V - var_L_type_Ca_channel__E_Ca_app);
const real var_membrane__i_Ca_L = var_L_type_Ca_channel__i_Ca_L;
const real var_Ca_independent_transient_outward_K_current__V = var_membrane__V;
const real var_Ca_independent_transient_outward_K_current__g_t = 8.25*mScaleFactorIto;
const real var_Ca_independent_transient_outward_K_current__F = var_membrane__F;
const real var_Ca_independent_transient_outward_K_current__K_c = var_cleft_space_ion_concentrations__K_c;
const real var_Ca_independent_transient_outward_K_current__T = var_membrane__T;
const real var_Ca_independent_transient_outward_K_current__K_i = var_intracellular_ion_concentrations__K_i;
const real var_Ca_independent_transient_outward_K_current__R = var_membrane__R;
const real var_Ca_independent_transient_outward_K_current__E_K = ((var_Ca_independent_transient_outward_K_current__R * var_Ca_independent_transient_outward_K_current__T) / var_Ca_independent_transient_outward_K_current__F) * log(var_Ca_independent_transient_outward_K_current__K_c / var_Ca_independent_transient_outward_K_current__K_i);
const real var_Ca_independent_transient_outward_K_current__r = var_Ca_independent_transient_outward_K_current_r_gate__r;
const real var_Ca_independent_transient_outward_K_current__s = var_Ca_independent_transient_outward_K_current_s_gate__s;
const real var_Ca_independent_transient_outward_K_current__i_t = var_Ca_independent_transient_outward_K_current__g_t * var_Ca_independent_transient_outward_K_current__r * var_Ca_independent_transient_outward_K_current__s * (var_Ca_independent_transient_outward_K_current__V - var_Ca_independent_transient_outward_K_current__E_K);
const real var_membrane__i_t = var_Ca_independent_transient_outward_K_current__i_t;
const real var_ultra_rapid_K_current__g_kur = 2.25*mScaleFactorGKur;
const real var_ultra_rapid_K_current__V = var_membrane__V;
const real var_ultra_rapid_K_current__E_K = var_Ca_independent_transient_outward_K_current__E_K;
const real var_ultra_rapid_K_current__a_ur = var_ultra_rapid_K_current_aur_gate__a_ur;
const real var_ultra_rapid_K_current__i_ur = var_ultra_rapid_K_current_iur_gate__i_ur;
const real var_ultra_rapid_K_current__i_Kur = var_ultra_rapid_K_current__g_kur * var_ultra_rapid_K_current__a_ur * var_ultra_rapid_K_current__i_ur * (var_ultra_rapid_K_current__V - var_ultra_rapid_K_current__E_K);
const real var_membrane__i_Kur = var_ultra_rapid_K_current__i_Kur;
const real var_inward_rectifier__T = var_membrane__T;
const real var_inward_rectifier__V = var_membrane__V;
const real var_inward_rectifier__g_K1 = 3.1*mScaleFactorGK1;
const real var_inward_rectifier__R = var_membrane__R;
const real var_inward_rectifier__K_c = var_cleft_space_ion_concentrations__K_c;
const real var_inward_rectifier__F = var_membrane__F;
const real var_inward_rectifier__E_K = var_Ca_independent_transient_outward_K_current__E_K;
const real var_inward_rectifier__i_K1 = (var_inward_rectifier__g_K1 * pow(var_inward_rectifier__K_c / 1.0, 0.4457) * (var_inward_rectifier__V - var_inward_rectifier__E_K)) / (1.0 + exp((1.5 * ((var_inward_rectifier__V - var_inward_rectifier__E_K) + 3.6) * var_inward_rectifier__F) / (var_inward_rectifier__R * var_inward_rectifier__T)));
const real var_membrane__i_K1 = var_inward_rectifier__i_K1;
const real var_delayed_rectifier_K_currents__V = var_membrane__V;
const real var_delayed_rectifier_K_currents_pi_gate__V = var_delayed_rectifier_K_currents__V;
const real var_delayed_rectifier_K_currents_pi_gate__pip = 1.0 / (1.0 + exp((var_delayed_rectifier_K_currents_pi_gate__V + 55.0) / 24.0));
const real var_delayed_rectifier_K_currents__pip = var_delayed_rectifier_K_currents_pi_gate__pip;
const real var_delayed_rectifier_K_currents_pa_gate__AZD = mScaleFactorAZD;
const real var_delayed_rectifier_K_currents__AZD = var_delayed_rectifier_K_currents_pa_gate__AZD;
const real var_delayed_rectifier_K_currents__factorAZD = (1.0 - (1.0 / (1.0 + exp((-(var_delayed_rectifier_K_currents__AZD - 0.6)) / 0.12)))) / 0.9933;
const real var_delayed_rectifier_K_currents__E_K = var_Ca_independent_transient_outward_K_current__E_K;
const real var_delayed_rectifier_K_currents__g_Kr = 0.5*mScaleFactorGkr;
const real var_delayed_rectifier_K_currents__pa = var_delayed_rectifier_K_currents_pa_gate__pa;
const real var_delayed_rectifier_K_currents__i_Kr = var_delayed_rectifier_K_currents__g_Kr * var_delayed_rectifier_K_currents__factorAZD * var_delayed_rectifier_K_currents__pa * var_delayed_rectifier_K_currents__pip * (var_delayed_rectifier_K_currents__V - var_delayed_rectifier_K_currents__E_K);
const real var_membrane__i_Kr = var_delayed_rectifier_K_currents__i_Kr;
const real var_delayed_rectifier_K_currents__n = var_delayed_rectifier_K_currents_n_gate__n;
const real var_delayed_rectifier_K_currents__g_Ks = 1.0*mScaleFactorGks;
const real var_delayed_rectifier_K_currents__i_Ks = var_delayed_rectifier_K_currents__g_Ks * var_delayed_rectifier_K_currents__n * (var_delayed_rectifier_K_currents__V - var_delayed_rectifier_K_currents__E_K);
const real var_membrane__i_Ks = var_delayed_rectifier_K_currents__i_Ks;
const real var_background_currents__E_Na = var_sodium_current__E_Na;
const real var_background_currents__g_B_Na = 0.060599;
const real var_background_currents__V = var_membrane__V;
const real var_background_currents__i_B_Na = var_background_currents__g_B_Na * (var_background_currents__V - var_background_currents__E_Na);
const real var_membrane__i_B_Na = var_background_currents__i_B_Na;
const real var_background_currents__Ca_c = var_cleft_space_ion_concentrations__Ca_c;
const real var_background_currents__R = var_membrane__R;
const real var_background_currents__Ca_i = var_intracellular_ion_concentrations__Ca_i;
const real var_background_currents__F = var_membrane__F;
const real var_background_currents__T = var_membrane__T;
const real var_background_currents__E_Ca = ((var_background_currents__R * var_background_currents__T) / (2.0 * var_background_currents__F)) * log(var_background_currents__Ca_c / var_background_currents__Ca_i);
const real var_background_currents__g_B_Ca = 0.078681;
const real var_background_currents__i_B_Ca = var_background_currents__g_B_Ca * (var_background_currents__V - var_background_currents__E_Ca);
const real var_membrane__i_B_Ca = var_background_currents__i_B_Ca;
const real var_sodium_potassium_pump__K_c = var_cleft_space_ion_concentrations__K_c;
const real var_sodium_potassium_pump__pow_K_NaK_Na_15 = 36.4829;
const real var_sodium_potassium_pump__Na_i = var_intracellular_ion_concentrations__Na_i;
const real var_sodium_potassium_pump__pow_Na_i_15 = pow(var_sodium_potassium_pump__Na_i, 1.5);
const real var_sodium_potassium_pump__V = var_membrane__V;
const real var_sodium_potassium_pump__i_NaK_max = 68.55*mScaleFactorGNaK;
const real var_sodium_potassium_pump__K_NaK_K = 1.0;
const real var_sodium_potassium_pump__i_NaK = (((((var_sodium_potassium_pump__i_NaK_max * var_sodium_potassium_pump__K_c) / (var_sodium_potassium_pump__K_c + var_sodium_potassium_pump__K_NaK_K)) * var_sodium_potassium_pump__pow_Na_i_15) / (var_sodium_potassium_pump__pow_Na_i_15 + var_sodium_potassium_pump__pow_K_NaK_Na_15)) * (var_sodium_potassium_pump__V + 150.0)) / (var_sodium_potassium_pump__V + 200.0);
const real var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
const real var_sarcolemmal_calcium_pump_current__i_CaP_max = 4.0;
const real var_sarcolemmal_calcium_pump_current__k_CaP = 0.0002;
const real var_sarcolemmal_calcium_pump_current__Ca_i = var_intracellular_ion_concentrations__Ca_i;
const real var_sarcolemmal_calcium_pump_current__i_CaP = (var_sarcolemmal_calcium_pump_current__i_CaP_max * var_sarcolemmal_calcium_pump_current__Ca_i) / (var_sarcolemmal_calcium_pump_current__Ca_i + var_sarcolemmal_calcium_pump_current__k_CaP);
const real var_membrane__i_CaP = var_sarcolemmal_calcium_pump_current__i_CaP;
const real var_Na_Ca_ion_exchanger_current__Ca_c = var_cleft_space_ion_concentrations__Ca_c;
const real var_Na_Ca_ion_exchanger_current__Na_i = var_intracellular_ion_concentrations__Na_i;
const real var_Na_Ca_ion_exchanger_current__gamma_Na = 0.45;
const real var_Na_Ca_ion_exchanger_current__Ca_i = var_intracellular_ion_concentrations__Ca_i;
const real var_Na_Ca_ion_exchanger_current__K_NaCa = 0.0374842*mScaleFactorGNaCa;
const real var_Na_Ca_ion_exchanger_current__F = var_membrane__F;
const real var_Na_Ca_ion_exchanger_current__d_NaCa = 0.0003;
const real var_Na_Ca_ion_exchanger_current__R = var_membrane__R;
const real var_Na_Ca_ion_exchanger_current__Na_c = var_cleft_space_ion_concentrations__Na_c;
const real var_Na_Ca_ion_exchanger_current__V = var_membrane__V;
const real var_Na_Ca_ion_exchanger_current__T = var_membrane__T;
const real var_Na_Ca_ion_exchanger_current__i_NaCa = (var_Na_Ca_ion_exchanger_current__K_NaCa * ((var_Na_Ca_ion_exchanger_current__Na_i * var_Na_Ca_ion_exchanger_current__Na_i * var_Na_Ca_ion_exchanger_current__Na_i * var_Na_Ca_ion_exchanger_current__Ca_c * exp((var_Na_Ca_ion_exchanger_current__F * var_Na_Ca_ion_exchanger_current__V * var_Na_Ca_ion_exchanger_current__gamma_Na) / (var_Na_Ca_ion_exchanger_current__R * var_Na_Ca_ion_exchanger_current__T))) - (var_Na_Ca_ion_exchanger_current__Na_c * var_Na_Ca_ion_exchanger_current__Na_c * var_Na_Ca_ion_exchanger_current__Na_c * var_Na_Ca_ion_exchanger_current__Ca_i * exp(((var_Na_Ca_ion_exchanger_current__gamma_Na - 1.0) * var_Na_Ca_ion_exchanger_current__V * var_Na_Ca_ion_exchanger_current__F) / (var_Na_Ca_ion_exchanger_current__R * var_Na_Ca_ion_exchanger_current__T))))) / (1.0 + (var_Na_Ca_ion_exchanger_current__d_NaCa * ((var_Na_Ca_ion_exchanger_current__Na_c * var_Na_Ca_ion_exchanger_current__Na_c * var_Na_Ca_ion_exchanger_current__Na_c * var_Na_Ca_ion_exchanger_current__Ca_i) + (var_Na_Ca_ion_exchanger_current__Na_i * var_Na_Ca_ion_exchanger_current__Na_i * var_Na_Ca_ion_exchanger_current__Na_i * var_Na_Ca_ion_exchanger_current__Ca_c))));
const real var_membrane__i_NaCa = var_Na_Ca_ion_exchanger_current__i_NaCa;
const real var_ACh_dependent_K_current__Cm = var_membrane__Cm;
const real var_ACh_dependent_K_current__V = var_membrane__V;
const real var_ACh_dependent_K_current__E_K = var_Ca_independent_transient_outward_K_current__E_K;
const real var_ACh_dependent_K_current__ACh = mScaleFactorAch;
const real var_ACh_dependent_K_current__i_KACh = (10000.0 / (1.0 + ((9.13652 * pow(1.0, 0.477811)) / pow(var_ACh_dependent_K_current__ACh, 0.477811)))) * (0.0517 + (0.4516 / (1.0 + exp((var_ACh_dependent_K_current__V + 59.53) / 17.18)))) * (var_ACh_dependent_K_current__V - var_ACh_dependent_K_current__E_K) * var_ACh_dependent_K_current__Cm;
const real var_membrane__i_KACh = var_ACh_dependent_K_current__i_KACh;
real var_membrane__i_Stim = stim_current;
const real var_membrane__I = var_membrane__i_Na + var_membrane__i_Ca_L + var_membrane__i_t + var_membrane__i_Kur + var_membrane__i_K1 + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_B_Na + var_membrane__i_B_Ca + var_membrane__i_NaK + var_membrane__i_CaP + var_membrane__i_NaCa + var_membrane__i_KACh + var_membrane__i_Stim;
const real var_sodium_current_m_gate__V = var_sodium_current__V;
const real var_sodium_current_m_gate__m_infinity = 1.0 / (1.0 + exp((var_sodium_current_m_gate__V + 27.12) / (-8.21)));
const real var_sodium_current_m_gate__m_factor = (var_sodium_current_m_gate__V + 25.57) / 28.8;
const real var_sodium_current_m_gate__tau_m = (4.2e-05 * exp((-var_sodium_current_m_gate__m_factor) * var_sodium_current_m_gate__m_factor)) + 2.4e-05;
const real var_sodium_current_h1_gate__V = var_sodium_current__V;
const real var_sodium_current_h1_gate__h_infinity = 1.0 / (1.0 + exp((var_sodium_current_h1_gate__V + 63.6) / 5.3));
const real var_sodium_current_h1_gate__h_factor = 1.0 / (1.0 + exp((var_sodium_current_h1_gate__V + 35.1) / 3.2));
const real var_sodium_current_h1_gate__tau_h1 = (0.03 * var_sodium_current_h1_gate__h_factor) + 0.0003;
const real var_sodium_current_h2_gate__h_infinity = var_sodium_current_h1_gate__h_infinity;
const real var_sodium_current_h2_gate__h_factor = var_sodium_current_h1_gate__h_factor;
const real var_sodium_current_h2_gate__tau_h2 = (0.12 * var_sodium_current_h2_gate__h_factor) + 0.003;
const real var_L_type_Ca_channel_d_L_gate__V = var_L_type_Ca_channel__V;
const real var_L_type_Ca_channel_d_L_gate__d_L_infinity = 1.0 / (1.0 + exp((var_L_type_Ca_channel_d_L_gate__V + 9.0) / (-5.8)));
const real var_L_type_Ca_channel_d_L_gate__d_L_factor = (var_L_type_Ca_channel_d_L_gate__V + 35.0) / 30.0;
const real var_L_type_Ca_channel_d_L_gate__tau_d_L = (0.0027 * exp((-var_L_type_Ca_channel_d_L_gate__d_L_factor) * var_L_type_Ca_channel_d_L_gate__d_L_factor)) + 0.002;
const real var_L_type_Ca_channel_f_L1_gate__V = var_L_type_Ca_channel__V;
const real var_L_type_Ca_channel_f_L1_gate__f_L_infinity = 1.0 / (1.0 + exp((var_L_type_Ca_channel_f_L1_gate__V + 27.4) / 7.1));
const real var_L_type_Ca_channel_f_L1_gate__f_L_factor = var_L_type_Ca_channel_f_L1_gate__V + 40.0;
const real var_L_type_Ca_channel_f_L1_gate__tau_f_L1 = (0.161 * exp((((-var_L_type_Ca_channel_f_L1_gate__f_L_factor) * var_L_type_Ca_channel_f_L1_gate__f_L_factor) / 14.4) / 14.4)) + 0.01;
const real var_L_type_Ca_channel_f_L2_gate__f_L_factor = var_L_type_Ca_channel_f_L1_gate__f_L_factor;
const real var_L_type_Ca_channel_f_L2_gate__tau_f_L2 = (1.3323 * exp((((-var_L_type_Ca_channel_f_L2_gate__f_L_factor) * var_L_type_Ca_channel_f_L2_gate__f_L_factor) / 14.2) / 14.2)) + 0.0626;
const real var_L_type_Ca_channel_f_L2_gate__f_L_infinity = var_L_type_Ca_channel_f_L1_gate__f_L_infinity;
const real var_Ca_independent_transient_outward_K_current_r_gate__V = var_Ca_independent_transient_outward_K_current__V;
const real var_Ca_independent_transient_outward_K_current_r_gate__tau_r = (0.0035 * exp((((-var_Ca_independent_transient_outward_K_current_r_gate__V) * var_Ca_independent_transient_outward_K_current_r_gate__V) / 30.0) / 30.0)) + 0.0015;
const real var_Ca_independent_transient_outward_K_current_r_gate__r_infinity = 1.0 / (1.0 + exp((var_Ca_independent_transient_outward_K_current_r_gate__V - 1.0) / (-11.0)));
const real var_Ca_independent_transient_outward_K_current_s_gate__V = var_Ca_independent_transient_outward_K_current__V;
const real var_Ca_independent_transient_outward_K_current_s_gate__s_factor = (var_Ca_independent_transient_outward_K_current_s_gate__V + 52.45) / 15.8827;
const real var_Ca_independent_transient_outward_K_current_s_gate__tau_s = (0.025635 * exp((-var_Ca_independent_transient_outward_K_current_s_gate__s_factor) * var_Ca_independent_transient_outward_K_current_s_gate__s_factor)) + 0.01414;
const real var_Ca_independent_transient_outward_K_current_s_gate__s_infinity = 1.0 / (1.0 + exp((var_Ca_independent_transient_outward_K_current_s_gate__V + 40.5) / 11.5));
const real var_ultra_rapid_K_current_aur_gate__V = var_ultra_rapid_K_current__V;
const real var_ultra_rapid_K_current_aur_gate__a_ur_infinity = 1.0 / (1.0 + exp((-(var_ultra_rapid_K_current_aur_gate__V + 6.0)) / 8.6));
const real var_ultra_rapid_K_current_aur_gate__tau_a_ur = (0.009 / (1.0 + exp((var_ultra_rapid_K_current_aur_gate__V + 5.0) / 12.0))) + 0.0005;
const real var_ultra_rapid_K_current_iur_gate__V = var_ultra_rapid_K_current__V;
const real var_ultra_rapid_K_current_iur_gate__i_ur_infinity = 1.0 / (1.0 + exp((var_ultra_rapid_K_current_iur_gate__V + 7.5) / 10.0));
const real var_ultra_rapid_K_current_iur_gate__tau_i_ur = (0.59 / (1.0 + exp((var_ultra_rapid_K_current_iur_gate__V + 60.0) / 10.0))) + 3.05;
const real var_delayed_rectifier_K_currents_n_gate__V = var_delayed_rectifier_K_currents__V;
const real var_delayed_rectifier_K_currents_n_gate__n_factor = (var_delayed_rectifier_K_currents_n_gate__V - 20.0) / 20.0;
const real var_delayed_rectifier_K_currents_n_gate__tau_n = 0.7 + (0.4 * exp((-var_delayed_rectifier_K_currents_n_gate__n_factor) * var_delayed_rectifier_K_currents_n_gate__n_factor));
const real var_delayed_rectifier_K_currents_n_gate__n_infinity = 1.0 / (1.0 + exp((var_delayed_rectifier_K_currents_n_gate__V - 19.9) / (-12.7)));
const real var_delayed_rectifier_K_currents_pa_gate__V = var_delayed_rectifier_K_currents__V;
const real var_delayed_rectifier_K_currents_pa_gate__pa_factor = (var_delayed_rectifier_K_currents_pa_gate__V + 20.1376) / 22.1996;
const real var_delayed_rectifier_K_currents_pa_gate__tau_pa = 0.03118 + (0.21718 * exp((-var_delayed_rectifier_K_currents_pa_gate__pa_factor) * var_delayed_rectifier_K_currents_pa_gate__pa_factor));
const real var_delayed_rectifier_K_currents_pa_gate__p_a_infinity = 1.0 / (1.0 + exp((var_delayed_rectifier_K_currents_pa_gate__V + 15.0 + (10.0 * var_delayed_rectifier_K_currents_pa_gate__AZD)) / (-6.0)));
const real var_intracellular_ion_concentrations__phi_Na_en = 0.0;
const real var_intracellular_ion_concentrations__Vol_i = 0.005884;
const real var_intracellular_ion_concentrations__Vol_d = 0.00011768;
const real var_intracellular_ion_concentrations__tau_di = 0.01;
const real var_intracellular_ion_concentrations__F = var_membrane__F;
const real var_intracellular_ion_concentrations__i_di = ((var_intracellular_ion_concentrations__Ca_d - var_intracellular_ion_concentrations__Ca_i) * 2.0 * var_intracellular_ion_concentrations__Vol_d * var_intracellular_ion_concentrations__F) / var_intracellular_ion_concentrations__tau_di;
const real var_intracellular_ion_concentrations__i_Na = var_sodium_current__i_Na;
const real var_intracellular_ion_concentrations__i_Ca_L = var_L_type_Ca_channel__i_Ca_L;
const real var_intracellular_ion_concentrations__i_t = var_Ca_independent_transient_outward_K_current__i_t;
const real var_intracellular_ion_concentrations__i_Kur = var_ultra_rapid_K_current__i_Kur;
const real var_intracellular_ion_concentrations__i_K1 = var_inward_rectifier__i_K1;
const real var_intracellular_ion_concentrations__i_Kr = var_delayed_rectifier_K_currents__i_Kr;
const real var_intracellular_ion_concentrations__i_Ks = var_delayed_rectifier_K_currents__i_Ks;
const real var_intracellular_ion_concentrations__i_B_Na = var_background_currents__i_B_Na;
const real var_intracellular_ion_concentrations__i_B_Ca = var_background_currents__i_B_Ca;
const real var_intracellular_ion_concentrations__i_NaK = var_sodium_potassium_pump__i_NaK;
const real var_intracellular_ion_concentrations__i_CaP = var_sarcolemmal_calcium_pump_current__i_CaP;
const real var_intracellular_ion_concentrations__i_NaCa = var_Na_Ca_ion_exchanger_current__i_NaCa;
const real var_intracellular_ion_concentrations__i_KACh = var_ACh_dependent_K_current__i_KACh;
const real var_Ca_handling_by_the_SR__k_cyca = 0.0003;
const real var_Ca_handling_by_the_SR__k_xcs = 0.4;
const real var_Ca_handling_by_the_SR__k_srca = 0.5;
const real var_Ca_handling_by_the_SR__Ca_i = var_intracellular_ion_concentrations__Ca_i;
const real var_Ca_handling_by_the_SR__I_up_max = 2800.0;
const real var_Ca_handling_by_the_SR__i_up = (var_Ca_handling_by_the_SR__I_up_max * ((var_Ca_handling_by_the_SR__Ca_i / var_Ca_handling_by_the_SR__k_cyca) - ((var_Ca_handling_by_the_SR__k_xcs * var_Ca_handling_by_the_SR__k_xcs * var_Ca_handling_by_the_SR__Ca_up) / var_Ca_handling_by_the_SR__k_srca))) / (((var_Ca_handling_by_the_SR__Ca_i + var_Ca_handling_by_the_SR__k_cyca) / var_Ca_handling_by_the_SR__k_cyca) + ((var_Ca_handling_by_the_SR__k_xcs * (var_Ca_handling_by_the_SR__Ca_up + var_Ca_handling_by_the_SR__k_srca)) / var_Ca_handling_by_the_SR__k_srca));
const real var_intracellular_ion_concentrations__i_up = var_Ca_handling_by_the_SR__i_up;
const real var_Ca_handling_by_the_SR__alpha_rel = 200000.0;
const real var_Ca_handling_by_the_SR__i_rel_f2 = var_Ca_handling_by_the_SR__F2 / (var_Ca_handling_by_the_SR__F2 + 0.25);
const real var_Ca_handling_by_the_SR__i_rel_factor = var_Ca_handling_by_the_SR__i_rel_f2 * var_Ca_handling_by_the_SR__i_rel_f2;
const real var_Ca_handling_by_the_SR__i_rel = var_Ca_handling_by_the_SR__alpha_rel * var_Ca_handling_by_the_SR__i_rel_factor * (var_Ca_handling_by_the_SR__Ca_rel - var_Ca_handling_by_the_SR__Ca_i);
const real var_intracellular_ion_concentrations__i_rel = var_Ca_handling_by_the_SR__i_rel;
const real var_intracellular_ion_concentrations__i_Stim = var_membrane__i_Stim;
const real var_intracellular_Ca_buffering__Ca_i = var_intracellular_ion_concentrations__Ca_i;
const real var_intracellular_Ca_buffering__J_O_TMgC = (200000.0 * var_intracellular_Ca_buffering__Ca_i * ((1.0 - var_intracellular_Ca_buffering__O_TMgC) - var_intracellular_Ca_buffering__O_TMgMg)) - (6.6 * var_intracellular_Ca_buffering__O_TMgC);
const real var_intracellular_Ca_buffering__J_O_TC = (78400.0 * var_intracellular_Ca_buffering__Ca_i * (1.0 - var_intracellular_Ca_buffering__O_TC)) - (392.0 * var_intracellular_Ca_buffering__O_TC);
const real var_intracellular_Ca_buffering__J_O_C = (200000.0 * var_intracellular_Ca_buffering__Ca_i * (1.0 - var_intracellular_Ca_buffering__O_C)) - (476.0 * var_intracellular_Ca_buffering__O_C);
const real var_intracellular_Ca_buffering__J_O = (0.08 * var_intracellular_Ca_buffering__J_O_TC) + (0.16 * var_intracellular_Ca_buffering__J_O_TMgC) + (0.045 * var_intracellular_Ca_buffering__J_O_C);
const real var_intracellular_ion_concentrations__J_O = var_intracellular_Ca_buffering__J_O;
const real var_intracellular_Ca_buffering__Mg_i = 2.5;
const real var_intracellular_Ca_buffering__J_O_TMgMg = (2000.0 * var_intracellular_Ca_buffering__Mg_i * ((1.0 - var_intracellular_Ca_buffering__O_TMgC) - var_intracellular_Ca_buffering__O_TMgMg)) - (666.0 * var_intracellular_Ca_buffering__O_TMgMg);
const real var_cleft_space_ion_concentrations__Vol_c = 0.000800224;
const real var_cleft_space_ion_concentrations__tau_Na = 14.3;
const real var_cleft_space_ion_concentrations__tau_K = 10.0;
const real var_cleft_space_ion_concentrations__tau_Ca = 24.7;
const real var_cleft_space_ion_concentrations__Na_b = 130.0;
const real var_cleft_space_ion_concentrations__Ca_b = 1.8;
const real var_cleft_space_ion_concentrations__K_b = 5.4;
const real var_cleft_space_ion_concentrations__F = var_membrane__F;
const real var_cleft_space_ion_concentrations__i_Na = var_sodium_current__i_Na;
const real var_cleft_space_ion_concentrations__i_Ca_L = var_L_type_Ca_channel__i_Ca_L;
const real var_cleft_space_ion_concentrations__i_t = var_Ca_independent_transient_outward_K_current__i_t;
const real var_cleft_space_ion_concentrations__i_Kur = var_ultra_rapid_K_current__i_Kur;
const real var_cleft_space_ion_concentrations__i_K1 = var_inward_rectifier__i_K1;
const real var_cleft_space_ion_concentrations__i_Kr = var_delayed_rectifier_K_currents__i_Kr;
const real var_cleft_space_ion_concentrations__i_Ks = var_delayed_rectifier_K_currents__i_Ks;
const real var_cleft_space_ion_concentrations__i_B_Na = var_background_currents__i_B_Na;
const real var_cleft_space_ion_concentrations__i_B_Ca = var_background_currents__i_B_Ca;
const real var_cleft_space_ion_concentrations__i_NaK = var_sodium_potassium_pump__i_NaK;
const real var_cleft_space_ion_concentrations__i_CaP = var_sarcolemmal_calcium_pump_current__i_CaP;
const real var_cleft_space_ion_concentrations__i_NaCa = var_Na_Ca_ion_exchanger_current__i_NaCa;
const real var_cleft_space_ion_concentrations__phi_Na_en = var_intracellular_ion_concentrations__phi_Na_en;
const real var_Ca_handling_by_the_SR__tau_tr = 0.01;
const real var_Ca_handling_by_the_SR__Vol_rel = 4.41e-05;
const real var_Ca_handling_by_the_SR__F = var_membrane__F;
const real var_Ca_handling_by_the_SR__i_tr = ((var_Ca_handling_by_the_SR__Ca_up - var_Ca_handling_by_the_SR__Ca_rel) * 2.0 * var_Ca_handling_by_the_SR__Vol_rel * var_Ca_handling_by_the_SR__F) / var_Ca_handling_by_the_SR__tau_tr;
const real var_Ca_handling_by_the_SR__Vol_up = 0.0003969;
const real var_Ca_handling_by_the_SR__k_rel_i = 0.0003;
const real var_Ca_handling_by_the_SR__r_Ca_i_term = var_Ca_handling_by_the_SR__Ca_i / (var_Ca_handling_by_the_SR__Ca_i + var_Ca_handling_by_the_SR__k_rel_i);
const real var_Ca_handling_by_the_SR__r_Ca_i_factor = var_Ca_handling_by_the_SR__r_Ca_i_term * var_Ca_handling_by_the_SR__r_Ca_i_term * var_Ca_handling_by_the_SR__r_Ca_i_term * var_Ca_handling_by_the_SR__r_Ca_i_term;
const real var_Ca_handling_by_the_SR__Ca_d = var_intracellular_ion_concentrations__Ca_d;
const real var_Ca_handling_by_the_SR__k_rel_d = 0.003;
const real var_Ca_handling_by_the_SR__r_Ca_d_term = var_Ca_handling_by_the_SR__Ca_d / (var_Ca_handling_by_the_SR__Ca_d + var_Ca_handling_by_the_SR__k_rel_d);
const real var_Ca_handling_by_the_SR__r_Ca_d_factor = var_Ca_handling_by_the_SR__r_Ca_d_term * var_Ca_handling_by_the_SR__r_Ca_d_term * var_Ca_handling_by_the_SR__r_Ca_d_term * var_Ca_handling_by_the_SR__r_Ca_d_term;
const real var_Ca_handling_by_the_SR__r_act = 203.8 * (var_Ca_handling_by_the_SR__r_Ca_i_factor + var_Ca_handling_by_the_SR__r_Ca_d_factor);
const real var_Ca_handling_by_the_SR__r_inact = 33.96 + (339.6 * var_Ca_handling_by_the_SR__r_Ca_i_factor);
const real var_Ca_handling_by_the_SR__r_recov = 0.815;
const real var_Ca_handling_by_the_SR__J_O_Calse = (480.0 * var_Ca_handling_by_the_SR__Ca_rel * (1.0 - var_Ca_handling_by_the_SR__O_Calse)) - (400.0 * var_Ca_handling_by_the_SR__O_Calse);

real d_dt_membrane__V = (-var_membrane__I) / var_membrane__Cm;

const real d_dt_sodium_current_m_gate__m = (var_sodium_current_m_gate__m_infinity - var_sodium_current_m_gate__m) / var_sodium_current_m_gate__tau_m;
const real d_dt_sodium_current_h1_gate__h1 = (var_sodium_current_h1_gate__h_infinity - var_sodium_current_h1_gate__h1) / var_sodium_current_h1_gate__tau_h1;
const real d_dt_sodium_current_h2_gate__h2 = (var_sodium_current_h2_gate__h_infinity - var_sodium_current_h2_gate__h2) / var_sodium_current_h2_gate__tau_h2;
const real d_dt_L_type_Ca_channel_d_L_gate__d_L = (var_L_type_Ca_channel_d_L_gate__d_L_infinity - var_L_type_Ca_channel_d_L_gate__d_L) / var_L_type_Ca_channel_d_L_gate__tau_d_L;
const real d_dt_L_type_Ca_channel_f_L1_gate__f_L1 = (var_L_type_Ca_channel_f_L1_gate__f_L_infinity - var_L_type_Ca_channel_f_L1_gate__f_L1) / var_L_type_Ca_channel_f_L1_gate__tau_f_L1;
const real d_dt_L_type_Ca_channel_f_L2_gate__f_L2 = (var_L_type_Ca_channel_f_L2_gate__f_L_infinity - var_L_type_Ca_channel_f_L2_gate__f_L2) / var_L_type_Ca_channel_f_L2_gate__tau_f_L2;
const real d_dt_Ca_independent_transient_outward_K_current_r_gate__r = (var_Ca_independent_transient_outward_K_current_r_gate__r_infinity - var_Ca_independent_transient_outward_K_current_r_gate__r) / var_Ca_independent_transient_outward_K_current_r_gate__tau_r;
const real d_dt_Ca_independent_transient_outward_K_current_s_gate__s = (var_Ca_independent_transient_outward_K_current_s_gate__s_infinity - var_Ca_independent_transient_outward_K_current_s_gate__s) / var_Ca_independent_transient_outward_K_current_s_gate__tau_s;
const real d_dt_ultra_rapid_K_current_aur_gate__a_ur = (var_ultra_rapid_K_current_aur_gate__a_ur_infinity - var_ultra_rapid_K_current_aur_gate__a_ur) / var_ultra_rapid_K_current_aur_gate__tau_a_ur;
const real d_dt_ultra_rapid_K_current_iur_gate__i_ur = (var_ultra_rapid_K_current_iur_gate__i_ur_infinity - var_ultra_rapid_K_current_iur_gate__i_ur) / var_ultra_rapid_K_current_iur_gate__tau_i_ur;
const real d_dt_delayed_rectifier_K_currents_n_gate__n = (var_delayed_rectifier_K_currents_n_gate__n_infinity - var_delayed_rectifier_K_currents_n_gate__n) / var_delayed_rectifier_K_currents_n_gate__tau_n;
const real d_dt_delayed_rectifier_K_currents_pa_gate__pa = (var_delayed_rectifier_K_currents_pa_gate__p_a_infinity - var_delayed_rectifier_K_currents_pa_gate__pa) / var_delayed_rectifier_K_currents_pa_gate__tau_pa;
const real d_dt_intracellular_ion_concentrations__K_i = (-(((var_intracellular_ion_concentrations__i_t + var_intracellular_ion_concentrations__i_Kur + var_intracellular_ion_concentrations__i_K1 + var_intracellular_ion_concentrations__i_Ks + var_intracellular_ion_concentrations__i_Kr + var_intracellular_ion_concentrations__i_KACh) - (2.0 * var_intracellular_ion_concentrations__i_NaK)) + var_intracellular_ion_concentrations__i_Stim)) / (var_intracellular_ion_concentrations__Vol_i * var_intracellular_ion_concentrations__F);
const real d_dt_intracellular_ion_concentrations__Na_i = (-(var_intracellular_ion_concentrations__i_Na + var_intracellular_ion_concentrations__i_B_Na + (3.0 * var_intracellular_ion_concentrations__i_NaCa) + (3.0 * var_intracellular_ion_concentrations__i_NaK) + var_intracellular_ion_concentrations__phi_Na_en)) / (var_intracellular_ion_concentrations__Vol_i * var_intracellular_ion_concentrations__F);
const real d_dt_intracellular_ion_concentrations__Ca_i = ((-((var_intracellular_ion_concentrations__i_B_Ca + var_intracellular_ion_concentrations__i_CaP + var_intracellular_ion_concentrations__i_up) - (var_intracellular_ion_concentrations__i_di + var_intracellular_ion_concentrations__i_rel + (2.0 * var_intracellular_ion_concentrations__i_NaCa)))) / (2.0 * var_intracellular_ion_concentrations__Vol_i * var_intracellular_ion_concentrations__F)) - (1.0 * var_intracellular_ion_concentrations__J_O);
const real d_dt_intracellular_ion_concentrations__Ca_d = (-(var_intracellular_ion_concentrations__i_Ca_L + var_intracellular_ion_concentrations__i_di)) / (2.0 * var_intracellular_ion_concentrations__Vol_d * var_intracellular_ion_concentrations__F);
const real d_dt_intracellular_Ca_buffering__O_C = var_intracellular_Ca_buffering__J_O_C;
const real d_dt_intracellular_Ca_buffering__O_TC = var_intracellular_Ca_buffering__J_O_TC;
const real d_dt_intracellular_Ca_buffering__O_TMgC = var_intracellular_Ca_buffering__J_O_TMgC;
const real d_dt_intracellular_Ca_buffering__O_TMgMg = var_intracellular_Ca_buffering__J_O_TMgMg;
const real d_dt_intracellular_Ca_buffering__O = var_intracellular_Ca_buffering__J_O;
const real d_dt_cleft_space_ion_concentrations__Ca_c = ((var_cleft_space_ion_concentrations__Ca_b - var_cleft_space_ion_concentrations__Ca_c) / var_cleft_space_ion_concentrations__tau_Ca) + (((var_cleft_space_ion_concentrations__i_Ca_L + var_cleft_space_ion_concentrations__i_B_Ca + var_cleft_space_ion_concentrations__i_CaP) - (2.0 * var_cleft_space_ion_concentrations__i_NaCa)) / (2.0 * var_cleft_space_ion_concentrations__Vol_c * var_cleft_space_ion_concentrations__F));
const real d_dt_cleft_space_ion_concentrations__K_c = ((var_cleft_space_ion_concentrations__K_b - var_cleft_space_ion_concentrations__K_c) / var_cleft_space_ion_concentrations__tau_K) + (((var_cleft_space_ion_concentrations__i_t + var_cleft_space_ion_concentrations__i_Kur + var_cleft_space_ion_concentrations__i_K1 + var_cleft_space_ion_concentrations__i_Ks + var_cleft_space_ion_concentrations__i_Kr) - (2.0 * var_cleft_space_ion_concentrations__i_NaK)) / (var_cleft_space_ion_concentrations__Vol_c * var_cleft_space_ion_concentrations__F));
const real d_dt_cleft_space_ion_concentrations__Na_c = ((var_cleft_space_ion_concentrations__Na_b - var_cleft_space_ion_concentrations__Na_c) / var_cleft_space_ion_concentrations__tau_Na) + ((var_cleft_space_ion_concentrations__i_Na + var_cleft_space_ion_concentrations__i_B_Na + (3.0 * var_cleft_space_ion_concentrations__i_NaCa) + (3.0 * var_cleft_space_ion_concentrations__i_NaK) + var_cleft_space_ion_concentrations__phi_Na_en) / (var_cleft_space_ion_concentrations__Vol_c * var_cleft_space_ion_concentrations__F));
const real d_dt_Ca_handling_by_the_SR__F1 = (var_Ca_handling_by_the_SR__r_recov * ((1.0 - var_Ca_handling_by_the_SR__F1) - var_Ca_handling_by_the_SR__F2)) - (var_Ca_handling_by_the_SR__r_act * var_Ca_handling_by_the_SR__F1);
const real d_dt_Ca_handling_by_the_SR__F2 = (var_Ca_handling_by_the_SR__r_act * var_Ca_handling_by_the_SR__F1) - (var_Ca_handling_by_the_SR__r_inact * var_Ca_handling_by_the_SR__F2);
const real d_dt_Ca_handling_by_the_SR__O_Calse = var_Ca_handling_by_the_SR__J_O_Calse;
const real d_dt_Ca_handling_by_the_SR__Ca_up = (var_Ca_handling_by_the_SR__i_up - var_Ca_handling_by_the_SR__i_tr) / (2.0 * var_Ca_handling_by_the_SR__Vol_up * var_Ca_handling_by_the_SR__F);
const real d_dt_Ca_handling_by_the_SR__Ca_rel = ((var_Ca_handling_by_the_SR__i_tr - var_Ca_handling_by_the_SR__i_rel) / (2.0 * var_Ca_handling_by_the_SR__Vol_rel * var_Ca_handling_by_the_SR__F)) - (31.0 * var_Ca_handling_by_the_SR__J_O_Calse);

rDY[0]  = 0.001*d_dt_membrane__V;
rDY[1]  = 0.001*d_dt_sodium_current_m_gate__m;
rDY[2]  = 0.001*d_dt_sodium_current_h1_gate__h1;
rDY[3]  = 0.001*d_dt_sodium_current_h2_gate__h2;
rDY[4]  = 0.001*d_dt_L_type_Ca_channel_d_L_gate__d_L;
rDY[5]  = 0.001*d_dt_L_type_Ca_channel_f_L1_gate__f_L1;
rDY[6]  = 0.001*d_dt_L_type_Ca_channel_f_L2_gate__f_L2;
rDY[7]  = 0.001*d_dt_Ca_independent_transient_outward_K_current_r_gate__r;
rDY[8]  = 0.001*d_dt_Ca_independent_transient_outward_K_current_s_gate__s;
rDY[9]  = 0.001*d_dt_ultra_rapid_K_current_aur_gate__a_ur;
rDY[10] = 0.001*d_dt_ultra_rapid_K_current_iur_gate__i_ur;
rDY[11] = 0.001*d_dt_delayed_rectifier_K_currents_n_gate__n;
rDY[12] = 0.001*d_dt_delayed_rectifier_K_currents_pa_gate__pa;
rDY[13] = 0.001*d_dt_intracellular_ion_concentrations__Na_i;
rDY[14] = 0.001*d_dt_intracellular_ion_concentrations__Ca_i;
rDY[15] = 0.001*d_dt_intracellular_ion_concentrations__K_i;
rDY[16] = 0.001*d_dt_intracellular_ion_concentrations__Ca_d;
rDY[17] = 0.001*d_dt_intracellular_Ca_buffering__O_C;
rDY[18] = 0.001*d_dt_intracellular_Ca_buffering__O_TC;
rDY[19] = 0.001*d_dt_intracellular_Ca_buffering__O_TMgC;
rDY[20] = 0.001*d_dt_intracellular_Ca_buffering__O_TMgMg;
rDY[21] = 0.001*d_dt_intracellular_Ca_buffering__O;
rDY[22] = 0.001*d_dt_cleft_space_ion_concentrations__Na_c;
rDY[23] = 0.001*d_dt_cleft_space_ion_concentrations__Ca_c;
rDY[24] = 0.001*d_dt_cleft_space_ion_concentrations__K_c;
rDY[25] = 0.001*d_dt_Ca_handling_by_the_SR__Ca_rel;
rDY[26] = 0.001*d_dt_Ca_handling_by_the_SR__Ca_up;
rDY[27] = 0.001*d_dt_Ca_handling_by_the_SR__O_Calse;
rDY[28] = 0.001*d_dt_Ca_handling_by_the_SR__F1;
rDY[29] = 0.001*d_dt_Ca_handling_by_the_SR__F2;