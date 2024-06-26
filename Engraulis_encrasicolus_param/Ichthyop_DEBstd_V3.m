%--------------------------------------------------------------------------
% DEBstd model implemented in Ichthyop
% Numerical integration : Euler
% Author                : Laure Pecquerie
% Modified              : J. Flores
% 2023/11/28

%% INITIALIZATION - TIME STEP, SIMULATION DURATION, VECTOR LENGTH
dt     = 0.0833;                        % d, time step of the model = 2h ==> dt = 2/24h
t_0    = 0;                             % d, January 1st
t_end  = 4*365 + t_0 - 1;               % d, last day of the simulation
n_iter = ceil((t_end - t_0 + 1) / dt);  % number of integration loop iterations
T_K    = 273.15;                        % Kelvin degrees

%% FORCING VARIABLES
temp  = 10:30;   % Temperaturas en C� para testear
f_res = 0.1 : 0.1 : 1; % Functional response to test
% We assume abundant food / ad libitum for now

%% PARAMETER VALUES

% Primary parameters
T_ref = 16 + T_K;      % K, Reference temperature (not to be changed) [Pethybridge et al 2013]
T_A   = 9800;          % K, Arrhenius temperature [Pethybridge et al 2013]

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

kap_X = 0.71;         % -, digestion efficiency of food to reserve [Pethybridge et al 2013]
p_Xm  = 325;          % J.cm-2.d-1 , Surface-area-specific maximum ingestion rate [Pethybridge et al 2013] Cambia con el tipo de alimento
p_Am  = kap_X * p_Xm; % 325 * 0.71 = 230.75 ;J.cm-2.d-1 , Surface-area-specific maximum assimilation rate [Pethybridge et al 2013]
E_m   = 2700;         % J.cm^(-3), maximum reserve density [Pethybridge et al 2013]
E_G   = 4000;         % J/cm^3, spec cost for structure [Pethybridge et al 2013] lo que hay que pagar para generar 1 cm3 de estructura
p_M   = 48;           % J/d.cm^3' , vol-spec somatic maintenance rate [Pethybridge et al 2013]
kap   = 0.7;          % - , allocation fraction to soma [Pethybridge et al 2013] Para mantenimiento y crecimiento
kap_R = 0.95;         % - , reproduction efficiency [Pethybridge et al 2013]
L_wb  = 0.25;         % cm, total length at mouth opening --> ?? total or fork length? [add my pet]
L_wp  = 9.077;        % cm, total length at puberty --> ?? guess [add my pet]

%% Auxiliary parameters
del_M = 0.154; % - , shape coefficient (Total Length)[Pethybridge et al 2013]
% del_M = 0.166; % - , Este valor produce una talla de 1.51 cm menos que la talla total, que es lo observado en figuras de ictiometro

%% Compound parameters
V_b = (L_wb * del_M)^3; % cm^3, structural volume at birth (first feeding)
V_p = (L_wp * del_M)^3; % cm^3, structural volume at puberty

%% Create a directory to store the results
subdir = 'C:/Users/jflores/Documents/JORGE/TESIS/TESIS_PHD/DEB/ichthyop_DEB/Engraulis_encrasicolus_param/DEBoutV3';
mkdir(subdir);

%% INITIAL CONDITIONS FOR THE STATE VARIABLES = EGG STAGE
E_in = [0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6];
V_0  = (0.0025 * del_M)^3; % cm, structural volume --> !! try different values
E_R0 = 0;                  % J, reproduction buffer

%% NUMERICAL INTEGRATION - EULER METHOD
t      = zeros(n_iter,1);
E      = zeros(n_iter,1);
V      = zeros(n_iter,1);
E_R    = zeros(n_iter,1);
F      = zeros(n_iter,1);
del    = zeros(n_iter,1);

p_M_flux_vec = zeros(n_iter,1); % J/d , Volume-related somatic maintenance
p_C_flux_vec = zeros(n_iter,1); % Energy for utilisation
p_G_flux_vec = zeros(n_iter,1); % J/d growth
p_J_flux_vec = zeros(n_iter,1); % J/d, Maturity maintenance !!

t(1)   = t_0;    % d,    Time vector initialization 
% E(1)   = E_0;    % J,    Initial reserve
V(1)   = V_0;    % cm^3, Initial structure
E_R(1) = E_R0;   % J,    Reproduction buffer   ------ I put an indice to try to run the script

for m = 1:size(E_in, 2)
    
    for j = 1:size(temp,2)
    
    for k = 1:size(f_res,2)
        
        T    = repmat(temp(j) + T_K, n_iter,1); % K, Temperature
        f    = f_res(k);
        E_0  = E_in(m);
        E(1) = E_0;
        
        for i = 1:n_iter-1

    %% Temperature correction
    % In case you want to use the complex temperature correction equation...
    s_A = exp(T_A/ T_ref - T_A ./ T(i));  % Arrhenius factor
    s_L_ratio = (1 + exp(T_AL/ T_ref - T_AL/ T_L)) ./ ...
	           (1 + exp(T_AL ./ T(i)   - T_AL/ T_L));
    s_H_ratio = (1 + exp(T_AH/ T_H - T_AH/ T_ref)) ./ ...
	           (1 + exp(T_AH/ T_H - T_AH ./ T(i)  ));
    c_T = s_A .* ((T(i) <= T_ref) .* s_L_ratio + (T(i) > T_ref) .* s_H_ratio); 

% 		%% Temperature correction
% 		% In case you want to use the simple temperature correction equation...
%         c_T   = exp(T_A/ T_ref - T_A ./ T(i));  % simple Arrhenius correction factor
		
		%% Correction of physiology parameters for temperature :
        p_AmT = c_T * p_Am;
        p_MT  = c_T * p_M;
 
        %% Scaled functional response
        % f = X(i) / (X(i) + K); % -, scaled functional response
        % f = 0.9999;

        %% Shape factor � std model
        del(i) = del_M;

        %% Fluxes , j/d
        if V(i) < V_b 
            p_A_flux = 0;
        else
            p_A_flux = f * p_AmT * V(i).^(2/3); % J/d assimilation rate 
        end

        p_M_flux = p_MT * V(i); % J/d , Volume-related somatic maintenance
        p_C_flux = E(i) .* ( E_G * p_AmT / E_m .* V(i)^(-1/3) + p_MT) ./ ( kap .* E(i) ./ V(i) + E_G); % utilization rate
        p_G_flux = kap * p_C_flux - p_M_flux;% J/d growth
        p_J_flux = (1- kap) / kap * p_MT * min(V(i), V_p)  ; % J/d, Maturity maintenance
        
        p_M_flux_vec(i) = p_M_flux; % J/d , Volume-related somatic maintenance
        p_C_flux_vec(i) = p_C_flux; % Energy for utilisation
        p_G_flux_vec(i) = p_G_flux; % J/d growth
        p_J_flux_vec(i) = p_J_flux; % J/d, Maturity maintenance !!

         if or( (p_G_flux < 0), ((1 - kap) * p_C_flux - p_J_flux <0))
%             disp('starvation')
%             return
            E(i) = 0; % When the value of E = 0, I can identify starvation in the output.
         else
		
            dE = p_A_flux - p_C_flux ; % J/d, 
            dV = p_G_flux / E_G;       % cm^3/d,
			
            if V(i) < V_p
                dE_R = 0;              % J/d, No reproduction buffer
            else  
                dE_R = (1 - kap) * p_C_flux - p_J_flux; % J/d, energy allocated to reproduction
            end 
         end

        E(i+1)   = E(i) + dE * dt ;     % J/d, Energy available into reserve
        V(i+1)   = V(i) + dV * dt ;     % cm^3, Volume of the structure
        E_R(i+1) = E_R(i) + dE_R * dt ; % J/d, Energy invested to reproduction
		
% 		%% Spawning rule 1: Every 30 days
% 		if (i/360 - round(i/360) == 0)
% 			E_R(i+1) = 0; % E_R regresa a cero
% 		end
%                
% 		%% Spawning rule 2: One spawning peak in September
% 		if (i == 3240 || i == 7560 || i == 11880 || i == 16200) % iterator indices for the September months
% 			E_R(i+1) = E_R(i) * 0.10; % El E_R es el 10% de la cantidad anterior
% 		end

		%% Spawning rule 3: Two spawning peaks in September & March
		if (i+1 == 3240 || i+1 == 7560 || i+1 == 11880 || i+1 == 16200 || ... % iterator indices for the September months (el desove se hace entre i y i+1)
			i+1 == 1080 || i+1 == 5400 || i+1 == 9720  || i+1 == 14040)       % iterator indices for the March months
			F(i+1)   = E_R(i+1) * kap_R / E_0;
			E_R(i+1) = 0; % E_R regresa a cero
		end
        
        t(i+1)   = t(i) + dt ;
        
        t_vec = repmat(temp(j), n_iter,1);  % C, Temperature
        f_vec = repmat(f_res(k), n_iter,1); % f, Functional response
        E_vec = repmat(E_in(m), n_iter,1); % E_0, Initial reserve
        
        end
	
        out_mat = table(t,E,V,E_R,F,t_vec,f_vec,E_vec,p_M_flux_vec,p_C_flux_vec,p_G_flux_vec,p_J_flux_vec,del,...
                        'VariableNames',...
                        {'t','E','V','E_R','F','temp','f','E0','p_M_flux','p_C_flux','p_G_flux','p_J_flux','delta'});
        writetable(out_mat, strcat(subdir, '/DEB_out','T',num2str(temp(j)),'f',num2str(f_res(k)),'E0_',num2str(E_in(m)),'.txt'))
    end
    end
end