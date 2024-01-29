%--------------------------------------------------------------------------
% DEBabj model implemented in Ichthyop
% Numerical integration : Euler
% Author                : Laure Pecquerie
% Modified              : J. Flores
% 2023/11/28

%% INITIALIZATION - TIME STEP, SIMULATION DURATION, VECTOR LENGTH
dt     = 0.0833;                        % d, time step of the model = 2h ==> dt = 2/24h
t_0    = 0;                             % d, January 1st
t_end  = 4*365 + t_0 - 1;               % d, last day of the simulation
n_iter = ceil((t_end - t_0 + 1) / dt);  % number of integration loop iterations
T_K    = 273.15;                        % Kelvin

%% FORCING VARIABLES
T = repmat(20 + T_K, n_iter,1); % K, Temperature
% X = repmat(20000, n_iter, 1);   % J/l, Food density (here in Joules per liter of water but can be different depending on the study)
% We assume abundant food / ad libitum for now

%% PARAMETER VALUES

% Primary parameters
T_ref = 20 + T_K;      % K, Reference temperature (not to be changed) [Pethybridge et al 2013]
T_A   = 10000;         % K, Arrhenius temperature [Pethybridge et al 2013]

% % In case you want to use the complex temperature correction equation...
% % Temperature correction - case 1
% T_L  = 6 + T_K;       % K, Lower boundary of the thermal range
% T_H  = 21 + T_K;      % K, Upper boundary of the thermal range
% T_AL = 20000;         % K, Arrhenius temperature at the lower boundary
% T_AH = 95000;         % K, Arrhenius temperature at the upper boundary
% 
% % Temperature correction - case 2
% T_L  = 6 + T_K;       % K, Lower boundary of the thermal range
% T_H  = 24 + T_K;      % K, Upper boundary of the thermal range
% T_AL = 20000;         % K, Arrhenius temperature at the lower boundary
% T_AH = 570000;        % K, Arrhenius temperature at the upper boundary

% if T_L > T_ref || T_H < T_ref
%      fprintf('Warning from temp_corr: invalid parameter combination, T_L > T_ref and/or T_H < T_ref\n')
% end

% kap_X = 0.71;    % -, digestion efficiency of food to reserve
% p_Xm  = 325;     % J.cm-2.d-1 , Surface-area-specific maximum ingestion rate [Pethybridge et al 2013] Cambia con el tipo de alimento
p_Am    = 84.97;   % J.cm-2.d-1 , Surface-area-specific maximum assimilation rate. In DEBstd = kap_X * p_Xm = 230.75
% E_m   = 2700;    % J.cm^(-3), maximum reserve density
V_dot   = 0.04124; % (cm d-1); Energy conductance
E_G     = 5283;    % J/cm^3, spec cost for structure [Pethybridge et al 2013] lo que hay que pagar para generar 1 cm3 de estructura
p_M     = 80.71;   % J/d.cm^3' , vol-spec somatic maintenance rate [Pethybridge et al 2013]
kap     = 0.5512;  % - , allocation fraction to soma [Pethybridge et al 2013] Para mantenimiento y crecimiento
kap_R   = 0.95;    % - , reproduction efficiency [Pethybridge et al 2013]
% L_wb  = 0.25;    % cm, total length at mouth opening --> ?? total or fork length? [add my pet]
% L_wp  = 9.077;   % cm, total length at puberty --> ?? guess [add my pet]
% L_b     = 0.0445;  % cm; Volumetric length at birth
% L_j     = 0.2612;  % cm, Volumetric length at metamorphosis
L_b     = 0.1038;  % cm; Volumetric length at birth (calculo JORGE)
L_j     = 0.6093;  % cm, Volumetric length at metamorphosis (calculo JORGE)
% E_Hb    = 0.3889;  % J, Maturity threshold at birth % ouverture de la bouche
E_Hb    = 0.335;  % J, Maturity threshold at birth % ouverture de la bouche % A 18 ºC se busca una edad de 5d
E_Hj    = 83.22;   % J, Maturity threshold at metamorphosis
E_Hp    = 42160;   % J, Maturity threshold at puberty
k_J     = 0.002;   % d-1, Maturity maintenance rate coefficient

% Auxiliary parameters
% del_M_t = 0.154; % - , shape coefficient (Total Length)[Pethybridge et al 2013]
% del_M_s = 0.166; % - , This value produces a standard length of 1.51 cm less than the total length, which is observed in ichthyometer figures.
del_M   = 0.1889  ;% -

% Compound parameters
% V_b = (L_wb * del_M_t)^3; % cm^3, structural volume at birth (first feeding)
% V_p = (L_wp * del_M_t)^3; % cm^3, structural volume at puberty

%% INITIAL CONDITIONS FOR THE STATE VARIABLES = EGG STAGE
E_0  = 1;         % J, egg content
V_0  = 0.0000001; % revisar valor (0.0025 * del_M_t)^3; % cm, structural volume --> !! try different values
E_H0 = 0;         % J, development
E_R0 = 0;         % J, reproduction buffer

%% NUMERICAL INTEGRATION - EULER METHOD
t      = zeros(n_iter,1);
E      = zeros(n_iter,1);
V      = zeros(n_iter,1);
E_H    = zeros(n_iter,1);
E_R    = zeros(n_iter,1);
F      = zeros(n_iter,1);

t(1)   = t_0;    % d,    Time vector initialization 
E(1)   = E_0;    % J,    Initial reserve
V(1)   = V_0;    % cm^3, Initial structure
E_R(1) = E_R0;   % J,    Reproduction buffer   ------ I put an indice to try to run the script
E_H(1) = E_H0;	 % J, Maturity

for i = 1:n_iter-1 

%     % Temperature correction
%     % In case you want to use the complex temperature correction equation...
%     s_A = exp(T_A/ T_ref - T_A ./ T(i));  % Arrhenius factor
%     s_L_ratio = (1 + exp(T_AL/ T_ref - T_AL/ T_L)) ./ ...
% 	           (1 + exp(T_AL ./ T(i)   - T_AL/ T_L));
%     s_H_ratio = (1 + exp(T_AH/ T_H - T_AH/ T_ref)) ./ ...
% 	           (1 + exp(T_AH/ T_H - T_AH ./ T(i)  ));
%     c_T = s_A .* ((T(i) <= T_ref) .* s_L_ratio + (T(i) > T_ref) .* s_H_ratio); 

    % Temperature correction
    % In case you want to use the simple temperature correction equation...
    c_T    = exp(T_A/ T_ref - T_A ./ T(i));  % simple Arrhenius correction factor
	
	% Correction of physiology parameters for temperature :
    p_AmT  = c_T * p_Am;
    V_dotT = c_T * V_dot;
    p_MT   = c_T * p_M;
	k_JT   = c_T * k_J;
 
    % Scaled functional response
%     f = X(i) / (X(i) + K); % -, scaled functional response
    f = 0.9999; % We assume food to satiation
%     f = 0.5; % We assume food to satiation

    %% Metabolic acceleration – abj model
	if E_H(i) < E_Hb
		s_M = 1;
	elseif (E_Hb <= E_H(i) && (E_H(i) < E_Hj))
%         s_M = (V(i)^(1/3))/L_b;
 		s_M = ((V(i)^(1/3))/del_M)/L_b; % Falta division del_M?
    else
		s_M = L_j / L_b;
	end
	
	% Only two parameters are accelerated by s_M, p_Am and v
	p_AmT  = s_M * p_AmT;
    V_dotT = s_M * V_dotT;
	
	% Fluxes , j/d
    if E_H(i) < E_Hb % antes DEBstd : V(i) < V_b 
        p_A_flux = 0;
    else
        p_A_flux = f * p_AmT * V(i).^(2/3); % J/d assimilation rate 
    end
    
    p_M_flux = p_MT * V(i); % J/d , Volume-related somatic maintenance
	p_C_flux = (E(i) / V(i)) * (E_G * V_dotT * V(i)^(2/3) + p_M_flux) ./ (kap * E(i) ./ V(i) + E_G); % Energy for utilisation
%    p_C_flux = E(i) .* ( E_G * p_AmT / E_m .* V(i)^(-1/3) + p_MT) ./ ( kap .* E(i) ./ V(i) + E_G); % antes DEBstd
    p_G_flux = max(0, kap * p_C_flux - p_M_flux);% J/d growth
  %	p_G_flux = kap * p_C_flux - p_M_flux;% antes DEBstd
    p_J_flux = k_JT * E_H(i)  ; % J/d, Maturity maintenance !!
  % p_J_flux = (1- kap) / kap * p_MT * min(V(i), V_p)  ; % antes DEBstd

  % if or( (p_G_flux < 0), ((1 - kap) * p_C_flux - p_J_flux <0)) % Antes DEBstd
    if or( kap * p_C_flux < p_M_flux, ((1 - kap) * p_C_flux - p_J_flux <0))
        fprintf('starvation = %d\n', t(i))
%         disp('starvation')
        return
    else

        dE = p_A_flux - p_C_flux; % J/d, 
        dV = ((kap * p_C_flux) - p_M_flux) / E_G;      % cm^3/d,

		if E_H(i) < E_Hp % antes DEBstd : V(i) < V_p
			dE_H = (1 - kap) * p_C_flux - p_J_flux; % J, Cumulated energy invested into development
			dE_R = 0; % J/d, No reproduction buffer
		else
			dE_H = 0; % J, No energy invested into development
			dE_R = (1 - kap) * p_C_flux - p_J_flux;  % J/d, Reproduction buffer
		end

		% Antes DEBstd
        %if V(i) < V_p
         %   dE_R = 0;             % J/d, No reproduction buffer
        %else
         %   dE_R = (1 - kap) * p_C_flux - p_J_flux; % J/d, energy allocated to reproduction
        %end 
    end
    
    E(i+1)   = E(i) + dE * dt;     % J/d, Energy available into reserve
    V(i+1)   = V(i) + dV * dt;     % cm^3, Volume of the structure
	E_H(i+1) = E_H(i) + dE_H * dt; % j/d, Energy invested to development
    E_R(i+1) = E_R(i) + dE_R * dt; % J/d, Energy invested to reproduction

%     % Spawning rule 1: Every 30 days
%     if (i/360 - round(i/360) == 0)
%         E_R(i+1) = 0; % E_R regresa a cero
%     end
               
%     % Spawning rule 2: One spawning peak in September
%     if (i == 3240 || i == 7560 || i == 11880 || i == 16200) % iterator indices for the September months
%         E_R(i+1) = E_R(i) * 0.10; % El E_R es el 10% de la cantidad anterior
%     end

    % Spawning rule 3: Two spawning peaks in September & March
    if (i+1 == 3240 || i+1 == 7560 || i+1 == 11880 || i+1 == 16200 || ... % iterator indices for the September months (el desove se hace entre i y i+1)
        i+1 == 1080 || i+1 == 5400 || i+1 == 9720  || i+1 == 14040)       % iterator indices for the March months
        F(i+1)   = E_R(i+1) * kap_R / E_0;
        E_R(i+1) = 0; % E_R regresa a cero
    end

    t(i+1) = t(i) + dt;
end

% Lab data
highT = [ % time since hatch days; standard length cm
0	0.3653
0	0.3675
0	0.3626
0	0.3585
0	0.3612
0	0.3404
0	0.3628
0	0.3336
0	0.3633
0	0.3582
0	0.3706
0	0.3729
0	0.369
0	0.3726
0	0.3901
0	0.3918
0	0.3662
0	0.2871
0	0.4025
0	0.2822
0	0.3786
0	0.3647
0	0.3681
4	0.4263
4	0.3992
4	0.4228
4	0.5208
4	0.5012
4	0.4351
4	0.4608
4	0.4001
4	0.4982
4	0.4852
5	0.5234
5	0.4865
5	0.5026
5	0.4799
5	0.5418
5	0.538
5	0.5585
11	0.8581
11	0.8803
11	0.891
11	0.7403
11	1.0529
11	0.8984
11	1.0262
11	0.84
11	0.9251
11	0.813
11	0.9109
11	0.9527
11	0.8612
12	0.9228
12	0.8507
12	0.6765
12	0.8326
12	0.9708
12	0.8124
12	0.7705
12	0.8151
12	0.7141
12	0.7001
19	1.2754
19	1.1911
19	1.1621
19	1.1283
19	1.246
19	1.3706
19	1.2236
19	1.0893
19	1.2838
19	1.1196
19	1.2051
26	1.3558
26	1.5027
26	1.4334
26	1.4733
26	1.2944
26	1.4202
26	1.6701
26	1.7034
26	1.849
26	1.7012
26	1.5582
];

lowT = [ % time since hatching d, standard length cm
8	0.4045
8	0.5254
8	0.5322
11	0.5343
11	0.5431
11	0.5799
17	0.6631
17	0.8126
17	0.7915
17	0.5038
17	0.5958
17	0.5449
24	1.0108
24	0.9449
24	0.6064
24	0.8379
24	1.0477
24	1.0781
24	1.1711
24	1.0785
24	1.0329
24	0.8123
31	1.1071
31	1.0117
31	1.2131
31	1.1086
31	1.1175
31	1.2745
31	1.2837
31	1.1877
31	1.134
31	1.2973
31	1.0665
31	0.9507
31	1.2379
];

%% OBSERVABLE VARIABLES

% Physical length
L_w = (V.^(1/3))./del_M; % cm, Physical length
% L_w = (V.^(1/3))./del_M_t; % antes DEBstd

% Wet weight
d_V  = 0.23;   % g/cm^3, specific density of structure (dry weight)
mu_V = 500000; % J/mol, specific chemical potential of structure
mu_E = 550000; % J/mol, specific chemical potential of reserve
w_V  = 23.9;   % g/mol, molecular weight of structure
w_E  = 23.9;   % g/mol, molecular weight of reserve
c_w  = 0.756;  % (c_w * W_w = total water weight)

W_V        = d_V.*V;
W_E        = (w_E/mu_E).*E;
W_ER       = (w_E/mu_E).*E_R;
Dry_weight = [W_V W_E W_ER];        % g, Dry weight
Wet_weight = Dry_weight./(1 - c_w); % g, Wet weight

Ww = sum(Wet_weight,2);

%% PLOTS
% State variables plots
% figure(1)
% 
% subplot(221)
% plot(t,E) % time in days VS Energy available into reserve in J 
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Reserve E (J)', 'Fontsize', 15)  
% 
% subplot (222)
% plot(t,V) % time in days VS Volume of the structure in cm^3
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Structure V (cm^3)', 'Fontsize', 15)  
% 
% subplot (223)
% plot(t,E_R) % time in days VS Energy invested to reproduction in J/d
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('E_R (J)', 'Fontsize', 15)  

% Observable variables plots
% figure(2)

% subplot(221)
% plot(t,L_w) % time in days VS Length in cm
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Length (cm)', 'Fontsize', 15)

% subplot(222)
% area(t,Wet_weight)
% legend('DEB Structure', 'DEB Reserve', 'DEB Reproduction')
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Wet Weight (g)', 'Fontsize', 15)
% 
% subplot(223)
% plot(L_w,Ww)
% xlabel('Length (cm)', 'Fontsize', 15)
% ylabel('Wet Weight (g)', 'Fontsize', 15)
% 
% subplot (224)
% plot(t, F) % time in days VS fecundity 
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Number of oocytes (#)', 'Fontsize', 15)

% Plot Fecundidad
% figure(3)
% 
% subplot(122)
% plot(t, F./7.5) % Por qué dividio por 7.5
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Relative Fecundidad (#eggs/batch)', 'Fontsize', 15)
% 
% subplot(121)
% plot(t,Ww)
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Wet Weight (g)', 'Fontsize', 15)
% 
% index = find(F > 0);
% figure(4)
% plot(Ww(index),F(index))
% xlabel('Wet Weight (g)', 'Fontsize', 15)
% ylabel('Relative Fecundidad (#eggs/batch)', 'Fontsize', 15)

% Length at transitions 
i_b = find(E_H >= E_Hb,1); 
fprintf('age at birth = %d\n', t(i_b) - t(1)) % d, age at birth
fprintf('length at birth = %d\n',L_w(i_b))    %cm, length at birth

i_j = find(E_H >= E_Hj,1); 
fprintf('age at metamorphosis = %d\n', t(i_j) - t(1)) % d, age at puberty 
fprintf('length at metamorphosis = %d\n',L_w(i_j))    %cm, length at puberty

i_p = find(E_H >= E_Hp,1); 
fprintf('age at puberty = %d\n', t(i_p) - t(1)) % d, age at puberty 
fprintf('length at puberty = %d\n',L_w(i_p))    %cm, length at puberty

i_l = find(L_w >= 2,1);
fprintf('edad a los 2 cm = %d\n', t(i_l) - t(1)) % d, age at puberty 

% figure(5)
% plot(t(1:i_l),L_w(1:i_l)) % time in days VS Length in cm
% line([t(i_b) t(i_b)],[0 2]);
% line([t(i_j) t(i_j)],[0 2]);
% hold on
% scatter( highT(:, 1)+2,  highT(:, 2))
% scatter( lowT(:, 1)+2,  lowT(:, 2))

% plot(t,L_w) % time in days VS Length in cm
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Length (cm)', 'Fontsize', 15)