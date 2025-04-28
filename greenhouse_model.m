%% Greenhouse Thermal Model with Heat and Mass Transfer
% Complete MATLAB script to solve the greenhouse thermal model
% Based on equations 16-51 with heat and mass transfer
 
%% Define parameters
%clear all;
close all;
clc;

% Basic parameters
T0 = 273.15;           % Reference temperature (K)
T_a_init = 25 + T0;    % Initial air temperature (K)
T_as_init = 25 + T0;   % Initial air above screen temperature (K) 
T_s_init = 20 + T0;    % Initial soil temperature (K)
T_c_init = 25 + T0;    % Initial crop temperature (K)
T_ri_init = 20 + T0;   % Initial inner roof temperature (K)
T_sc_init = 22 + T0;   % Initial screen temperature (K)
T_ro_init = 20 + T0;   % Initial outer roof temperature (K)
T_sk_init = 15 + T0;   % Initial sky temperature (K)
C_e_H2O_a_init = 0.010;% Initial water vapor concentration in air (kg/m³)
C_e_H2O_as_init = 0.008;% Initial water vapor concentration above screen (kg/m³)
C_e_H2O_c_init = 0.015;% Initial water vapor concentration at crop (kg/m³)
C_e_H2O_ri_init = 0.007;% Initial water vapor concentration at inner roof (kg/m³)
C_e_H2O_sc_init = 0.009;% Initial water vapor concentration at screen (kg/m³)

% Define the input parameters (from your shared list)
eta_ri_Is = 0.0173;  % η_{ri-ls}
E_s = 0.7;           % E_s
E_sk = 0.8;          % E_{sk}
F_ri_sk = 0.86;      % F_{ri-sk}
A_nw = 11.52;        % A_{nw}
lambda_nw = 0.397;   % λ_{nw}
d_nw = 0.25;         % d_{nw}
cp_r = 840;          % c_{p-r}
E_ri = 0.95;         % E_{ri}
LAI = 1;             % LAI
A_s = 15.36;         % A_s
A_c = 2 * LAI * A_s; % A_c = 2 * LAI * A_s
A_sc = 15.36;        % A_{sc}
A_r = 17.7;          % A_r
E_sc = 0.9;          % E_{sc}
F_ro_sk = A_s / A_r; % F_{ro-sk}
rho_c = 700;         % ρ_c
R_s_min = 82.003;    % R_{min}
d_sc = 0.002;        % d_{sc}
V_sc = 3.07;         % V_{sc}
cp_a = 1000;         % c_{p-a}
f_a = 1;             % f_a
V_r = 0.0708;        % V_r
V_a = 26.4;          % V_a
V_s = 9.984;         % V_s
ro_r = 2500;         % ρ_r
cp_s = 800;          % c_{p-s}
ro_s = 1400;         % ρ_s
F_s_ri = 0.8;        % F_{s-ri}
eta_s_Is = 0.86;     % η_{ls-sk}
ro_H2O = 998;        % ρ_{H2O}
v_a = 0.09;          % v_a
eta_c_Is = 0.5;      
V_as = 12.25;        % V_{as}
E_ro = 0.95;         % E_{ro}
ro_sc = 200;         % ρ_sc
cp_H2O = 4186;       % c_{p-H2O}
cp_sc = 1500;        % c_{p-sc}
sigma = 5.67051e-8;  % σ (Stefan-Boltzmann constant)
F_sc_ri = 0.8;       % F_{sc-ri}
L_e = 2.45e6;        % Latent heat of vaporization (J/kg)
Cl_sc = 0.9;         % Screen closure coefficient
lambda_s = 0.6;      % λ_s
d_s = 0.65;          % d_s

% Weather-related parameters
I_0 = 600;           % Solar radiation intensity (W/m²)
v_o = 1.0;           % Outside wind velocity (m/s)
RH_a = 60;           % Relative humidity (%)
T_o = 30 + T0;       % Outside temperature (K)
I_s = 450;           % Solar radiation on soil (W/m²)
I_r = 550;           % Solar radiation on roof (W/m²)
I_in = 500;          % Indoor solar radiation (W/m²)
J_s = 1.0;           % Factor for Rb_heat calculation

% Screen status (0 = no screen, 1 = screen)
csc = 1;
v_as = 0.1;       % Air velocity above screen (m/s)
cp_c = 3500;      % Specific heat of crop (J/kg·K)
V_c = 0.2 * V_a;  % Crop volume (m³)


%% Define function for equations 16-51
% Create function to calculate all coefficients and heat/mass transfers

%% Set up simulation
tspan = [0 24*3600];  % Simulate for 24 hours
y0 = [T_a_init; T_as_init; T_s_init; T_c_init; T_ri_init; T_sc_init; T_ro_init; 
      C_e_H2O_a_init; C_e_H2O_as_init; C_e_H2O_c_init; C_e_H2O_ri_init; C_e_H2O_sc_init];
options = odeset('RelTol', 1e-1, 'AbsTol', 1e-1);

%% Define the ODE system
function dydt = greenhouse_ode_system(~, y, ~)
    % Extract state variables
   % Basic parameters
T0 = 273.15;           % Reference temperature (K)
T_a_init = 25 + T0;    % Initial air temperature (K)
T_as_init = 25 + T0;   % Initial air above screen temperature (K) 
T_s_init = 20 + T0;    % Initial soil temperature (K)
T_c_init = 25 + T0;    % Initial crop temperature (K)
T_ri_init = 20 + T0;   % Initial inner roof temperature (K)
T_sc_init = 22 + T0;   % Initial screen temperature (K)
T_ro_init = 20 + T0;   % Initial outer roof temperature (K)
T_sk_init = 15 + T0;   % Initial sky temperature (K)
C_e_H2O_a_init = 0.010;% Initial water vapor concentration in air (kg/m³)
C_e_H2O_as_init = 0.008;% Initial water vapor concentration above screen (kg/m³)
C_e_H2O_c_init = 0.015;% Initial water vapor concentration at crop (kg/m³)
C_e_H2O_ri_init = 0.007;% Initial water vapor concentration at inner roof (kg/m³)
C_e_H2O_sc_init = 0.009;% Initial water vapor concentration at screen (kg/m³)

% Define the input parameters (from your shared list)
eta_ri_Is = 0.0173;  % η_{ri-ls}
E_s = 0.7;           % E_s
E_sk = 0.8;          % E_{sk}
F_ri_sk = 0.86;      % F_{ri-sk}
A_nw = 11.52;        % A_{nw}
lambda_nw = 0.397;   % λ_{nw}
d_nw = 0.25;         % d_{nw}
cp_r = 840;          % c_{p-r}
E_ri = 0.95;         % E_{ri}
LAI = 1;             % LAI
A_s = 15.36;         % A_s
A_c = 2 * LAI * A_s; % A_c = 2 * LAI * A_s
A_sc = 15.36;        % A_{sc}
A_r = 17.7;          % A_r
E_sc = 0.9;          % E_{sc}
F_ro_sk = A_s / A_r; % F_{ro-sk}
rho_c = 700;         % ρ_c
R_s_min = 82.003;    % R_{min}
d_sc = 0.002;        % d_{sc}
V_sc = 3.07;         % V_{sc}
cp_a = 1000;         % c_{p-a}
f_a = 1;             % f_a
V_r = 0.0708;        % V_r
V_a = 26.4;          % V_a
V_s = 9.984;         % V_s
ro_r = 2500;         % ρ_r
cp_s = 800;          % c_{p-s}
ro_s = 1400;         % ρ_s
F_s_ri = 0.8;        % F_{s-ri}
eta_s_Is = 0.86;     % η_{ls-sk}
ro_H2O = 998;        % ρ_{H2O}
v_a = 0.09;          % v_a
eta_c_Is = 0.5;      
V_as = 12.25;        % V_{as}
E_ro = 0.95;         % E_{ro}
ro_sc = 200;         % ρ_sc
cp_H2O = 4186;       % c_{p-H2O}
cp_sc = 1500;        % c_{p-sc}
sigma = 5.67051e-8;  % σ (Stefan-Boltzmann constant)
F_sc_ri = 0.8;       % F_{sc-ri}
L_e = 2.45e6;        % Latent heat of vaporization (J/kg)
Cl_sc = 0.9;         % Screen closure coefficient
lambda_s = 0.6;      % λ_s
d_s = 0.65;          % d_s

% Weather-related parameters
I_0 = 600;           % Solar radiation intensity (W/m²)
v_o = 1.0;           % Outside wind velocity (m/s)
RH_a = 60;           % Relative humidity (%)
T_o = 30 + T0;       % Outside temperature (K)
I_s = 450;           % Solar radiation on soil (W/m²)
I_r = 550;           % Solar radiation on roof (W/m²)
I_in = 500;          % Indoor solar radiation (W/m²)
J_s = 1.0;           % Factor for Rb_heat calculation

% Screen status (0 = no screen, 1 = screen)
csc = 1;
v_as = 0.1;       % Air velocity above screen (m/s)
cp_c = 3500;      % Specific heat of crop (J/kg·K)
V_c = 0.2 * V_a;  % Crop volume (m³)

    T_a = y(1);      % Air temperature
    T_as = y(2);     % Air above screen temperature
    T_s = y(3);      % Soil temperature
    T_c = y(4);      % Crop temperature
    T_ri = y(5);     % Inside roof temperature
    T_sc = y(6);     % Screen temperature
    T_ro = y(7);     % Outside roof temperature
    C_e_H2O_a = y(8);   % Water vapor concentration in air
    C_e_H2O_as = y(9);  % Water vapor concentration above screen
    C_e_H2O_c = y(10);  % Water vapor concentration at crop surface
    C_e_H2O_ri = y(11); % Water vapor concentration at inner roof
    C_e_H2O_sc = y(12); % Water vapor concentration at screen

    % Extract parameters from params structure
    % [Code to extract all the parameters]
    
    % Update temperature-dependent parameters
  ro_a = 1.29 * T0 / max(T_a,1e-3); 
ro_as = 1.29 * T0 / max(T_as,1e-3);
    %p_av = 1.29 * (2 * T0) / (T_a + T_as); % Average air density
    T_sk = 0.0552 * (T_s)^1.5; % Sky temperature (Eq. 34)
    
    % Calculate heat transfer coefficients (Eqs. 16-21)
    % Equation (16)
    denominator = max((J_s * abs(T_c - T_a)) + 207 * v_a^2, 1e-3);
Rb_heat = 1174 / sqrt(I_s) / denominator;% Eq. 23
    alpha_a_c = ro_a * cp_a / Rb_heat;

    % Equation (17)
    if T_a < T_s
        alpha_a_s = 1.7 * abs(T_a - T_s)^(1/3);
    else
        alpha_a_s = 1.3 * abs(T_a - T_s)^0.25;
    end

    % Equation (18)
    alpha_as_ri = 3 * abs(T_as - T_ri)^(1/3);

    % Equation (19)
    if v_o < 4
        alpha_ro_o = 2.8 + 1.2 * v_o;
    else
        alpha_ro_o = 2.5 * v_o^0.8;
    end

    % Equation (20)
    alpha_a_sc = Cl_sc * 3 * abs(T_a - T_sc)^(1/3);

    % Equation (21)
    alpha_as_sc = Cl_sc * 3 * abs(T_as - T_sc)^(1/3);
    
    % Calculate mass transfer coefficients (Eqs. 36-43)
    % Equation (36)
    u_a = v_a; % Assuming air velocity equals v_a
    r_b_c = 200 * (u_a)^(-0.5); % Eq. 37
    
    % Equation (38-41)
    f_I = I_s / (I_s + 0.54); % Eq. 39
    
    if T_a <= 3
        f_T = 1 + 0.05 * (T_a - 336);
    else
        f_T = 1 + 2.9593 * 10^(-3) * (T_a - 24.512);
    end
    
    f_H = 1 / sqrt(1 + 255 * exp(-0.136 * RH_a)); % Eq. 41
    r_s_c = R_s_min / (f_I * f_T * f_H); % Eq. 38
    k_c_a_H2O = 1 / (r_b_c + r_s_c); % Eq. 36
    
    % Equation (42-43)
    Le = 0.89; % Lewis number
    k_a_as_H2O = (alpha_a_sc / (ro_a * cp_a)) * Le^(-2/3); % Eq. 47
    k_as_ri_H2O = (alpha_as_ri / (ro_a * cp_a)) * Le^(-2/3); % Eq. 48
    
    % Calculate mass flow rates (Eqs. 45-46, 51)
    Phi_m_as_ri_H2O = max(A_sc * k_as_ri_H2O * (C_e_H2O_as - C_e_H2O_ri), 0); % Eq. 45
    Phi_m_a_as_H2O = max(A_sc * k_a_as_H2O * (C_e_H2O_a - C_e_H2O_as), 0); % Eq. 46
    
    % Calculate f_as_sc (volumetric airflow rate) based on pressure difference
    f_as_sc = 0.1; % Simplified estimation
    Phi_m_as_sc_H2O = f_as_sc * (C_e_H2O_as - C_e_H2O_sc); % Eq. 51
    
    % Calculate latent heat transfers (Eqs. 44, 49, 50)
    %Q_e_as_ri = ro_as * L_e * Phi_m_as_ri_H2O; % Eq. 44
    %Q_e_a_as = ro_a * L_e * Phi_m_a_as_H2O; % Eq. 49
    %Q_e_as_sc = ro_as * L_e * Phi_m_as_sc_H2O; % Eq. 50
    
    % Water flux from crop to air (Eq. 35)
    Phi_m_c_a_H2O = max(A_c * k_c_a_H2O * (C_e_H2O_as - C_e_H2O_c), 0);
    
    % Calculate radiative heat exchanges
    Q_r_ri = A_r * eta_ri_Is * I_r; % Eq. 24
    Q_r_c = A_c * eta_c_Is * I_in; % Eq. 25
    Q_r_s = A_s * eta_s_Is * I_in; % Eq. 26
    
    % Define missing variables and calculate radiative heat exchanges
    tau_cIl = exp(-0.64 * LAI); % Transmission coefficient
    E_c = 1 - tau_cIl; % Crop emissivity
    
    F_s_c = 1 - tau_cIl; % View factor soil to crop
    F_c_sc = Cl_sc * (1 - tau_cIl); % View factor crop to screen
    
    % Calculate all other heat transfers
    Qa_c = A_c * alpha_a_c * (T_a - T_c);
    Qa_s = A_s * alpha_a_s * (T_a - T_s);
    Qas_ri = A_r * alpha_as_ri * (T_as - T_ri);
    Qro_o = A_r * alpha_ro_o * (T_ro - T_o);
    Qa_sc = A_sc * alpha_a_sc * (T_a - T_sc);
    Qas_sc = A_sc * alpha_as_sc * (T_as - T_sc);
    
    % Air exchange between compartments
    Phi_as = v_as * A_sc * (1 - Cl_sc); % Eq. 22
    Qa_as = ro_a * cp_a * Phi_as * (T_a - T_as);
    
    % Conduction through soil
    Qs_ss = A_s * (lambda_s / d_s) * (T_s - T_s); % Here T_s - T_s should be replaced with your actual soil depth temperature difference
    
    % Conduction through north wall (if applicable)
    Qnwi_nwo = A_nw * (lambda_nw / d_nw) * (T_a - T_o); % Simplification
    
    % Radiative heat exchanges (Eqs. 27-33)
    Q_s_c = A_s * E_s * E_c * F_s_c * sigma * (T_s^4 - T_c^4); % Eq. 27
    Q_ro_sk = A_r * E_ro * E_sk * F_ro_sk * sigma * (T_ro^4 - T_sk^4); % Eq. 32
    
    % Additional radiative exchanges needed
    Q_s_ri = A_s * E_s * E_ri * F_s_ri * sigma * (T_s^4 - T_ri^4);
    Q_s_sc = A_s * E_s * E_sc * Cl_sc * sigma * (T_s^4 - T_sc^4);
    Q_sc_ri = A_sc * E_sc * E_ri * F_sc_ri * sigma * (T_sc^4 - T_ri^4); % Eq. 31
    Q_ri_c = A_r * E_ri * E_c * 0.8 * sigma * (T_ri^4 - T_c^4); % Assuming F_ri_c = 0.8
    
    % Latent heat from crop to air
    Q_c_a_H2O = ro_a * L_e * Phi_m_c_a_H2O;
    
    % Heat from crop to screen
    Q_c_sc = A_c * E_c * E_sc * F_c_sc * sigma * (T_c^4 - T_sc^4);
    
    % Heat from screen to air
    Q_sc_a = A_sc * alpha_a_sc * (T_sc - T_a);
    
    % Heat from screen to airspace
    Q_sc_as = A_sc * alpha_as_sc * (T_sc - T_as);
    
    % Radiation from outside to roof
    Q_rd_ri = Q_r_ri;
    
    % Radiation to soil
    Q_rd_s = Q_r_s;
    
    % Radiation to crop
    Q_rd_c = Q_r_c;
    
    % Latent heat transfer from air to screen
    Q_a_sc_H2O = ro_a * L_e * max(A_sc * k_a_as_H2O * (C_e_H2O_a - C_e_H2O_sc), 0);
    
    % Calculate derivatives (temperature change rates)
    if csc == 0
        dTa_dt = (Qa_s - Qa_c - Qas_ri - Qnwi_nwo) / (ro_a * cp_a * V_a);
        dTas_dt = dTa_dt;  % Same as air temperature when no screen
    else
        dTa_dt = (Qa_sc + Qas_sc - Qa_c - Qa_s - Qa_as - Qnwi_nwo) / (ro_a * cp_a * V_a);
        dTas_dt = (Qa_as + Qas_sc - Qas_ri - Qnwi_nwo) / (ro_as * cp_a * V_as);
    end
    
    dTs_dt = (Q_rd_s + Qa_s - Q_s_c - Q_s_ri - Qs_ss - Q_s_sc) / ((0.7 * ro_s * cp_s + 0.2 * ro_H2O * cp_H2O + 0.1 * ro_a * cp_a) * V_s);
    
    % For crop, include cp_c (specific heat) and V_c (volume) - assume these values
    cp_c = 3500; % Specific heat of crop (J/kg·K)
    V_c = 0.2 * V_a; % Volume of crop (m³)
    dTc_dt = (Q_rd_c + Qa_c + Q_ri_c + Q_s_c - Q_c_a_H2O - Q_c_sc) / (rho_c * cp_c * V_c);
    
    dTri_dt = (Q_rd_ri + Qas_ri + Q_s_ri + Q_sc_ri - Q_ri_c - Qro_o - Q_ro_sk) / (ro_r * cp_r * V_r);
    
    dTsc_dt = (Q_c_sc + Q_s_sc + Q_a_sc_H2O - Q_sc_a - Q_sc_as - Q_sc_ri) / (ro_sc * cp_sc * V_sc);
    
    % Outer roof temperature (simplified)
    dTro_dt = (Qro_o + Q_ro_sk) / (ro_r * cp_r * V_r * 0.5);
    
    % Water vapor concentration derivatives (simplified mass balance)
    dC_e_H2O_a_dt = (Phi_m_c_a_H2O - Phi_m_a_as_H2O) / V_a;
    dC_e_H2O_as_dt = (Phi_m_a_as_H2O - Phi_m_as_ri_H2O - Phi_m_as_sc_H2O) / V_as;
    
    % Simplified water vapor concentration at surfaces (assuming equilibrium)
    dC_e_H2O_c_dt = 0;  % Assume crop water content is maintained
    dC_e_H2O_ri_dt = 0; % Assume condensation on roof maintains equilibrium
    dC_e_H2O_sc_dt = 0; % Assume screen water content reaches equilibrium
    
    % Return derivatives
    dydt = [dTa_dt; dTas_dt; dTs_dt; dTc_dt; dTri_dt; dTsc_dt; dTro_dt;
            dC_e_H2O_a_dt; dC_e_H2O_as_dt; dC_e_H2O_c_dt; dC_e_H2O_ri_dt; dC_e_H2O_sc_dt];
end

%% Simulate the model using ode45
tspan = [0 24*3600];  % Time span (24 hours)
y0 = [T_a_init; T_as_init; T_s_init; T_c_init; T_ri_init; T_sc_init; T_ro_init;
      C_e_H2O_a_init; C_e_H2O_as_init; C_e_H2O_c_init; C_e_H2O_ri_init; C_e_H2O_sc_init];

% Solve the system of ODEs
options = odeset('RelTol',1e-2,'AbsTol',1e-2,'MaxStep',3600);  % Adjust tolerances
[t, y] = ode15s(@greenhouse_ode_system, tspan, y0,options);

%% Plot Results
% Plot temperature and humidity profiles over time
figure;
subplot(2,1,1);
plot(t/3600, y(:,1) - T0, 'DisplayName', 'Air Temp (°C)');
hold on;
plot(t/3600, y(:,4) - T0, 'DisplayName', 'Crop Temp (°C)');
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend;

subplot(2,1,2);
plot(t/3600, y(:,8), 'DisplayName', 'Water Vapor (Air) (kg/m³)');
hold on;
plot(t/3600, y(:,10), 'DisplayName', 'Water Vapor (Crop) (kg/m³)');
xlabel('Time (hours)');
ylabel('Water Vapor Concentration (kg/m³)');
legend;


