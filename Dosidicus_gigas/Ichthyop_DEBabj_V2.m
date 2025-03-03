%--------------------------------------------------------------------------
% DEBabj model implemented in Ichthyop
% Numerical integration : Euler
% Author                : J. Flores
% Modified              : J. Flores
% 2025/03/03

%% INITIALIZATION - TIME STEP, SIMULATION DURATION, VECTOR LENGTH
dt     = 0.0833;                        % d, time step of the model = 2h ==> dt = 2/24h
t_0    = 0;                             % d, January 1st
t_end  = 1.5*365 + t_0 - 1;               % d, last day of the simulation
n_iter = ceil((t_end - t_0 + 1) / dt);  % number of integration loop iterations
T_K    = 273.15;                        % Kelvin degrees

%% FORCING VARIABLES
temp  = 10:30;          % Temperaturas en C� para testear
f_res = 0.1 : 0.1 : 1;  % Functional response to test

%% PARAMETER VALUES

% Primary parameters
T_ref = 20 + T_K;       % K, Reference temperature
T_A   = 9000;           % K, Arrhenius temperature

% In case you want to use the complex temperature correction equation...
% Temperature correction - case 1
T_L  = 6 + T_K;       % K, Lower boundary of the thermal range
T_H  = 21 + T_K;      % K, Upper boundary of the thermal range
T_AL = 20000;         % K, Arrhenius temperature at the lower boundary
T_AH = 95000;         % K, Arrhenius temperature at the upper boundary
%  
% % Temperature correction - case 2
% T_L  = 6 + T_K;       % K, Lower boundary of the thermal range
% T_H  = 24 + T_K;      % K, Upper boundary of the thermal range
% T_AL = 20000;         % K, Arrhenius temperature at the lower boundary
% T_AH = 570000;        % K, Arrhenius temperature at the upper boundary
%  
% if T_L > T_ref || T_H < T_ref
%      fprintf('Warning from temp_corr: invalid parameter combination, T_L > T_ref and/or T_H < T_ref\n')
% end

p_Am    = 1950.6;  % J.cm-2.d-1 , Surface-area-specific maximum assimilation rate. Male = 1839.5
V_dot   = 0.02134; % (cm d-1); Energy conductance
E_G     = 5500;    % J/cm^3, spec cost for structure [Pethybridge et al 2013] lo que hay que pagar para generar 1 cm3 de estructura
p_M     = 345.04;  % J/d.cm^3' , vol-spec somatic maintenance rate [Pethybridge et al 2013]
kap     = 0.7354;  % - , allocation fraction to soma [Pethybridge et al 2013] Para mantenimiento y crecimiento
kap_R   = 0.95;    % - , reproduction efficiency [Pethybridge et al 2013]
L_b     = 0.0299;  % cm; Volumetric length at birth
L_j     = 0.2234;  % cm, Volumetric length at metamorphosis
E_Hb    = 0.0566;  % J, Maturity threshold at birth % ouverture de la bouche % A 18.5 �C, hatch = 5d
E_Hj    = 24.09;   % J, Maturity threshold at metamorphosis
E_Hp    = 232700;  % J, Maturity threshold at puberty
k_J     = 0.002;   % d-1, Maturity maintenance rate coefficient
del_M1  = 0.1254;  % -, shape coefficient for standard length of larvae
% del_M2  = 0.1879;  % -, shape coefficient for standard length

% For a dynamic shape factor (Jusup et al 2011)
E_Hy    = E_Hj;    % J, Maturity at the end of the early juvenile stage ('metamorphosis')
E_H2    = (E_Hb + E_Hj)/2; % Half-saturation maturity, i.e. the level of maturity at which the shape factor is an arithmetic mean of del_M1 and del_M2 

%% Create a directory to store the results
subdir = 'C:/Users/jflores/Documents/JORGE/AFIIMM/models/Dosidicus_gigas/DEBoutV2/';
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

for j = 1:size(temp,2)

    for k = 1:size(f_res,2)

        T = repmat(temp(j) + T_K, n_iter,1); % K, Temperature
        f = f_res(k); % Scaled functional response

        for i = 1:n_iter-1
        %% Temperature correction
        % In case you want to use the complex temperature correction equation...
        s_A = exp(T_A/ T_ref - T_A / T(i));  % Arrhenius factor
        s_L_ratio = (1 + exp(T_AL/ T_ref - T_AL/ T_L)) / ...
                   (1 + exp(T_AL / T(i)   - T_AL/ T_L));
        s_H_ratio = (1 + exp(T_AH/ T_H - T_AH/ T_ref)) / ...
                   (1 + exp(T_AH/ T_H - T_AH / T(i)));
        c_T = s_A * ((T(i) <= T_ref) * s_L_ratio + (T(i) > T_ref) * s_H_ratio); 
 
% 		%% Temperature correction
% 		% In case you want to use the simple temperature correction equation...
%         c_T   = exp(T_A/ T_ref - T_A / T(i));  % simple Arrhenius correction factor
        
		%% Correction of physiology parameters for temperature :
		p_AmT  = c_T * p_Am;
		V_dotT = c_T * V_dot;
		p_MT   = c_T * p_M;
		k_JT   = c_T * k_J;
 
        %% Metabolic acceleration � abj model
        if E_H(i) < E_Hb
            s_M = 1;
        elseif (E_Hb <= E_H(i) && E_H(i) < E_Hj)
            s_M = (V(i)^(1/3))/L_b;
        else
            s_M = L_j / L_b;
        end
        
        %% Shape factor � abj model
        if E_H(i) < E_Hb
            del_M_conti = del_M1; % shape coefficient for standard length of larvae
        elseif (E_Hb <= E_H(i) && E_H(i) < E_Hj)
            del_M_conti = del_M1; % shape coefficient for standard length of larvae
        else
            del_M_conti = del_M1; % shape coefficient for standard length
        end
        
%         %% Shape factor � abj model (Pecquerie 2024)
%         if E_H(i) < E_Hb
%             del_M_pecqu = del_M1; % shape coefficient for standard length of larvae
%         elseif (E_Hb <= E_H(i) && E_H(i) < E_Hj)
%             del_M_pecqu = (del_M1 * (E_Hj - E_H(i)) + del_M2 * (E_H(i) - E_Hb)) / (E_Hj - E_Hb); % shape coefficient for standard length of larvae
%         else
%             del_M_pecqu = del_M2; % shape coefficient for standard length
%         end
        
%         %% Shape factor � abj model (Jusup et al 2011)
%         if E_H(i) < E_Hb
%            del_M_jusup = del_M1; % shape coefficient for standard length of larvae
%         elseif (E_Hb <= E_H(i) && E_H(i) < E_Hy)
%            del_M_jusup = ( del_M1*(E_H2 - E_Hb) + del_M2*(E_H(i)- E_Hb) ) / (E_H(i) + E_H2 - 2*E_Hb);
%         else
%            del_M_jusup = del_M2; % shape coefficient for standard length
%         end
    
        del(i) = del_M_conti;
%         del(i) = del_M_pecqu;
%         del(i) = del_M_jusup;
        acc(i)  = s_M;
        
		% Only two parameters are accelerated by s_M: p_Am and V_dot
		p_AmT  = s_M * p_AmT;
		V_dotT = s_M * V_dotT;
		
        %% Fluxes , j/d
        if E_H < E_Hb
            p_A_flux = 0;
        else
            p_A_flux = f * p_AmT * V(i).^(2/3); % J/d assimilation rate 
        end

        p_M_flux = p_MT * V(i); % J/d , Volume-related somatic maintenance
        p_C_flux = (E(i) / V(i)) * (E_G * V_dotT * V(i)^(2/3) + p_M_flux) ./ (kap * E(i) ./ V(i) + E_G); % Energy for utilisation
        p_J_flux = k_JT * E_H(i)  ; % J/d, Maturity maintenance !!

        if or( kap * p_C_flux < p_M_flux, ((1 - kap) * p_C_flux - p_J_flux <0))
            fprintf('starvation = %d\n', t(i))
%             disp('starvation')
            return
        else		
            dE = p_A_flux - p_C_flux ; % J/d, 
            dV = ((kap * p_C_flux) - p_M_flux) / E_G;      % cm^3/d,       % cm^3/d,
			
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
		F(i+1)   = E_R(i+1) * kap_R / E_0;
        
%     %% Spawning rule
%     if (IGS > 20)
%         fprintf('Spawns and dies = %d\n', t(i))
%     end              
        
        t(i+1)   = t(i) + dt ;
        
        t_vec = repmat(temp(j), n_iter,1);  % C, Temperature
        f_vec = repmat(f_res(k), n_iter,1); % f, Functional response
        
        end
        
        % del1 = continuo; % del2 = pecquerie; % del3 = jusup
        out_mat = table(t,E,V,E_H,E_R,F,acc,del,t_vec,f_vec,...
                        'VariableNames',...
                        {'t','E','V','E_H','E_R','Fec','acc','delta','temp','f'});
        writetable(out_mat, strcat(subdir, '/DEB_out','T',num2str(temp(j)),'f',num2str(f_res(k)),'.txt'))
    end
end
