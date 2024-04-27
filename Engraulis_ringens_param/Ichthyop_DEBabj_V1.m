%--------------------------------------------------------------------------
% DEBabj model implemented in Ichthyop
% Numerical integration : Euler
% Author                : Laure Pecquerie
% Modified              : J. Flores
% 2024/03/17

%% INITIALIZATION - TIME STEP, SIMULATION DURATION, VECTOR LENGTH
dt     = 0.0833;                        % d, time step of the model = 2h ==> dt = 2/24h
t_0    = 0;                             % d, January 1st
t_end  = 4*365 + t_0 - 1;               % d, last day of the simulation
n_iter = ceil((t_end - t_0 + 1) / dt);  % number of integration loop iterations
T_K    = 273.15;                        % Kelvin degrees

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

p_Am    = 84.97;   % J.cm-2.d-1 , Surface-area-specific maximum assimilation rate. In DEBstd = kap_X * p_Xm = 230.75
V_dot   = 0.04124; % (cm d-1); Energy conductance
E_G     = 5283;    % J/cm^3, spec cost for structure [Pethybridge et al 2013] lo que hay que pagar para generar 1 cm3 de estructura
p_M     = 80.71;   % J/d.cm^3' , vol-spec somatic maintenance rate [Pethybridge et al 2013]
kap     = 0.5512;  % - , allocation fraction to soma [Pethybridge et al 2013] Para mantenimiento y crecimiento
kap_R   = 0.95;    % - , reproduction efficiency [Pethybridge et al 2013]
L_b     = 0.0445;  % cm; Volumetric length at birth
L_j     = 0.2612;  % cm, Volumetric length at metamorphosis
E_Hb    = 0.3350;  % J, Maturity threshold at birth % ouverture de la bouche % A 18.5 ºC, hatch = 5d
% E_Hb    = 0.3889;  % J, Maturity threshold at birth % ouverture de la bouche
E_Hj    = 83.22;   % J, Maturity threshold at metamorphosis
E_Hp    = 42160;   % J, Maturity threshold at puberty
k_J     = 0.002;   % d-1, Maturity maintenance rate coefficient
del_M1  = 0.08095; % -, shape coefficient for standard length of larvae
del_M2  = 0.1889;  % -, shape coefficient for standard length

% % For a dynamic shape factor (Jusup et al 2011)
% E_Hy    = E_Hp;    % J, Maturity at the end of the early juvenile stage
% E_H2    = (E_Hb + E_Hj)/2; % Half-saturation maturity, i.e. the level of maturity at which the shape factor is an arithmetic mean of del_M1 and del_M2

%% Create a directory to store the results
subdir = 'C:/Users/jflores/Documents/JORGE/TESIS/TESIS_PHD/DEB/ichthyop_DEB/Engraulis_ringens_param/DEBoutV1/';
mkdir(subdir);

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
acc    = zeros(n_iter,1);
del    = zeros(n_iter,1);

t(1)   = t_0;    % d,    Time vector initialization 
E(1)   = E_0;    % J,    Initial reserve
V(1)   = V_0;    % cm^3, Initial structure
E_R(1) = E_R0;   % J,    Reproduction buffer   ------ I put an indice to try to run the script
E_H(1) = E_H0;	 % J, Maturity

for i = 1:n_iter-1 

%     %% Temperature correction
%     % In case you want to use the complex temperature correction equation...
%     s_A = exp(T_A/ T_ref - T_A ./ T(i));  % Arrhenius factor
%     s_L_ratio = (1 + exp(T_AL/ T_ref - T_AL/ T_L)) ./ ...
% 	           (1 + exp(T_AL ./ T(i)   - T_AL/ T_L));
%     s_H_ratio = (1 + exp(T_AH/ T_H - T_AH/ T_ref)) ./ ...
% 	           (1 + exp(T_AH/ T_H - T_AH ./ T(i)  ));
%     c_T = s_A .* ((T(i) <= T_ref) .* s_L_ratio + (T(i) > T_ref) .* s_H_ratio); 

    %% Temperature correction
    % In case you want to use the simple temperature correction equation...
    c_T    = exp(T_A/ T_ref - T_A ./ T(i));  % simple Arrhenius correction factor
	
	%% Correction of physiology parameters for temperature :
    p_AmT  = c_T * p_Am;
    V_dotT = c_T * V_dot;
    p_MT   = c_T * p_M;
	k_JT   = c_T * k_J;

    %% Scaled functional response
%     f = X(i) / (X(i) + K); % -, scaled functional response
    f = 0.9999; % We assume food to satiation
   
    %% Metabolic acceleration – abj model
    if E_H(i) < E_Hb
        s_M = 1;
    elseif (E_Hb <= E_H(i) && E_H(i) < E_Hj)
        s_M = (V(i)^(1/3))/L_b;
    else
    	s_M = L_j / L_b;
    end

    %% Shape factor – abj model
    if E_H(i) < E_Hb
       del_M = del_M1; % shape coefficient for standard length of larvae
    elseif (E_Hb <= E_H(i) && E_H(i) < E_Hj)
       del_M = del_M1; % shape coefficient for standard length of larvae
    else
       del_M = del_M2; % shape coefficient for standard length
    end

%     %% Shape factor – abj model (Jusup et al 2011)
%     if E_H(i) < E_Hb
%        del_M = del_M1; % shape coefficient for standard length of larvae
%     elseif (E_Hb <= E_H(i) && E_H(i) < E_Hy)
%        del_M = ( del_M1*(E_H2 - E_Hb) + del_M2*(E_H(i)- E_Hb) ) / (E_H(i) + E_H2 - 2*E_Hb);
%     else
%        del_M = del_M2; % shape coefficient for standard length
%     end

	acc(i) = s_M;
    del(i) = del_M;

	% Only two parameters are accelerated by s_M: p_Am and V_dot
	p_AmT  = s_M * p_AmT;
    V_dotT = s_M * V_dotT;
	
	%% Fluxes , j/d
    if E_H(i) < E_Hb
        p_A_flux = 0;
    else
        p_A_flux = f * p_AmT * V(i).^(2/3); % J/d assimilation rate 
    end
    
    p_M_flux = p_MT * V(i); % J/d , Volume-related somatic maintenance
	p_C_flux = (E(i) / V(i)) * (E_G * V_dotT * V(i)^(2/3) + p_M_flux) ./ (kap * E(i) ./ V(i) + E_G); % Energy for utilisation
    p_J_flux = k_JT * E_H(i)  ; % J/d, Maturity maintenance !!

    if or( kap * p_C_flux < p_M_flux, ((1 - kap) * p_C_flux - p_J_flux <0))
        fprintf('starvation = %d\n', t(i))
%         disp('starvation')
        return
    else
        dE = p_A_flux - p_C_flux; % J/d, 
        dV = ((kap * p_C_flux) - p_M_flux) / E_G; % cm^3/d,

        if E_H(i) < E_Hp
			dE_H = (1 - kap) * p_C_flux - p_J_flux; % J, Cumulated energy invested into development
			dE_R = 0; % J/d, No reproduction buffer
		else
			dE_H = 0; % J, No energy invested into development
			dE_R = (1 - kap) * p_C_flux - p_J_flux;  % J/d, Reproduction buffer
        end
    end
    
    E(i+1)   = E(i) + dE * dt;     % J/d, Energy available into reserve
    V(i+1)   = V(i) + dV * dt;     % cm^3, Volume of the structure
	E_H(i+1) = E_H(i) + dE_H * dt; % j/d, Energy invested to development
    E_R(i+1) = E_R(i) + dE_R * dt; % J/d, Energy invested to reproduction

%     %% Spawning rule 1: Every 30 days
%     if (i/360 - round(i/360) == 0)
%         E_R(i+1) = 0; % E_R resets to zero
%     end
               
%     %% Spawning rule 2: One spawning peak in September
%     if (i == 3240 || i == 7560 || i == 11880 || i == 16200) % Indices for September months' iterators (the spawning is done between i and i+1)
%         E_R(i+1) = E_R(i) * 0.10; % E_R is 10 percent of the previous amount.
%     end

    %% Spawning rule 3: Two spawning peaks in September & March
    if (i+1 == 3240 || i+1 == 7560 || i+1 == 11880 || i+1 == 16200 || ... % Indices for September months' iterators (the spawning is done between i and i+1)
        i+1 == 1080 || i+1 == 5400 || i+1 == 9720  || i+1 == 14040)       % Indices for March months' iterators (the spawning is done between i and i+1)
        F(i+1)   = E_R(i+1) * kap_R / E_0;
        E_R(i+1) = 0; % E_R resets to zero
    end

    t(i+1) = t(i) + dt;
end

out_mat = table(t,E,V,E_H,E_R,F,T-T_K,acc,del,...
                        'VariableNames',...
                        {'t','E','V','E_H','E_R','Fec','temp','acc','delta'});
writetable(out_mat, strcat(subdir, 'DEB_out.txt'))

%% Physical length
Lw = (V.^(1/3))./del; % cm, Physical length

%% Wet weight
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
% 
% %% PLOTS
% 
% % % State variables plots
% % 
% % figure(1)
% % subplot(221)
% % plot(t,E) % time in days VS Energy available into reserve in J 
% % xlabel('time (d)', 'Fontsize', 15)
% % ylabel('Reserve E (J)', 'Fontsize', 15)  
% % 
% % subplot (222)
% % plot(t,V) % time in days VS Volume of the structure in cm^3
% % xlabel('time (d)', 'Fontsize', 15)
% % ylabel('Structure V (cm^3)', 'Fontsize', 15)  
% % 
% % subplot (223)
% % plot(t,E_R) % time in days VS Energy invested to reproduction in J/d
% % xlabel('time (d)', 'Fontsize', 15)
% % ylabel('E_R (J)', 'Fontsize', 15)  
% % 
% % subplot (224)
% % plot(t,E_H) % time in days VS Maturity level J/d
% % xlabel('time (d)', 'Fontsize', 15)
% % ylabel('E_H (J)', 'Fontsize', 15)  
% 
% % % Observable variables plots
% 
% figure(2)
% subplot(221)
% plot(t,Lw) % time in days VS Length in cm
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Length (cm)', 'Fontsize', 15)
% % 
% % subplot(222)
% % area(t,Wet_weight)
% % legend('DEB Structure', 'DEB Reserve', 'DEB Reproduction')
% % xlabel('time (d)', 'Fontsize', 15)
% % ylabel('Wet Weight (g)', 'Fontsize', 15)
% % 
% % subplot(223)
% % plot(Lw,Ww)
% % xlabel('Length (cm)', 'Fontsize', 15)
% % ylabel('Wet Weight (g)', 'Fontsize', 15)
% % 
% % subplot (224)
% % plot(t, F) % time in days VS fecundity 
% % xlabel('time (d)', 'Fontsize', 15)
% % ylabel('Number of oocytes (#)', 'Fontsize', 15)
% % 
% % % Plot Fecundidad
% % 
% % figure(3)
% % subplot(131)
% % plot(t, F./7.5) % Por qué dividio por 7.5
% % xlabel('time (d)', 'Fontsize', 15)
% % ylabel('Relative Fecundidad (#eggs/batch)', 'Fontsize', 15)
% % 
% % subplot(132)
% % plot(t,Ww)
% % xlabel('time (d)', 'Fontsize', 15)
% % ylabel('Wet Weight (g)', 'Fontsize', 15)
% % 
% % index = find(F > 0);
% % subplot(133)
% % plot(Ww(index),F(index))
% % xlabel('Wet Weight (g)', 'Fontsize', 15)
% % ylabel('Relative Fecundidad (#eggs/batch)', 'Fontsize', 15)
% % 
%% Length at transitions 
i_b = find(E_H >= E_Hb,1); 
fprintf('age at birth = %d\n', t(i_b) - t(1)) % d, age at birth
fprintf('length at birth = %d\n',Lw(i_b))    %cm, length at birth

i_j = find(E_H >= E_Hj,1); 
fprintf('age at metamorphosis = %d\n', t(i_j) - t(1)) % d, age at puberty 
fprintf('length at metamorphosis = %d\n',Lw(i_j))    %cm, length at puberty

% i_y = find(E_H >= E_Hy, 1);
% fprintf('age at early juvenile stage = %d\n', t(i_y) - t(1)) % d, age at maturity at the end of the early juvenile stage
% fprintf('length at early juvenile stage = %d\n',Lw(i_y))    %cm, length at maturity at the end of the early juvenile stage

i_p = find(E_H >= E_Hp,1); 
fprintf('age at puberty = %d\n', t(i_p) - t(1)) % d, age at puberty 
fprintf('length at puberty = %d\n',Lw(i_p))    %cm, length at puberty
% 
% 
% 
% i_2cm = find(Lw >= 2,1);
% fprintf('age at 2 cm = %d\n', t(i_2cm) - t(1)) % d, age at puberty 
% 
% figure(4)
% plot(t,acc) % time in days VS acceleration factor
% xlim([0 150])
% ylim([0 6])
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Acceleration Factor', 'Fontsize', 15)
%  
% figure(5)
% plot(t,del) % time in days VS shape coefficient
% xlim([0 200])
% % ylim([0 1])
% xlabel('-', 'Fontsize', 15)
% ylabel('Shape factor', 'Fontsize', 15)
% 
% figure(6)
% plot(t(1:i_2cm),Lw(1:i_2cm), 'Color', 'black') % time in days VS Length in cm
% line([t(i_b) t(i_b)],[0 2],'Color', 'green');
% line([t(i_j) t(i_j)],[0 2],'Color', 'green');
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Length (cm)', 'Fontsize', 15)
% hold on
% scatter(highT(:, 1)+2,  highT(:, 2), 'filled', 'red')
% scatter(lowT(:, 1)+2,  lowT(:, 2), 'filled', 'blue')
% 
% figure(7)
% plot(t(1:i_y),del(1:i_y), 'Color', 'black') % time in days VS Length in cm
% line([t(i_b) t(i_b)],[0.08 0.2],'Color', 'green');
% line([t(i_j) t(i_j)],[0.08 0.2],'Color', 'green');
% xlabel('time (d)', 'Fontsize', 15)
% ylabel('Shape coefficient (-)', 'Fontsize', 15)
% % hold on
% % scatter(highT(:, 1)+2,  highT(:, 2), 'filled', 'red')
% % scatter(lowT(:, 1)+2,  lowT(:, 2), 'filled', 'blue')
% 
% % figure(4)
% % plot(t,Lw) % time in days VS Length in cm
% % xlim([0 70])
% % xlabel('time (d)', 'Fontsize', 15)
% % ylabel('Length (cm)', 'Fontsize', 15)
% 