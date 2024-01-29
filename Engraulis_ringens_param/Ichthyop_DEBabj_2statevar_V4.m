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
T_K    = 273.15;                        % Kelvin degrees

%% FORCING VARIABLES
temp  = 10 : 2 : 30;   % Temperaturas en Cº para testear
f_res = 0 : 0.1 : 1; % Functional response to test
% We assume abundant food / ad libitum for now
talla = 3; % Talla de la larva en cm antes de iniciar una condicion de 'f' diferente a 1

%% PARAMETER VALUES

% Primary parameters
T_ref = 20 + T_K;      % K, Reference temperature (not to be changed) [Pethybridge et al 2013]
T_A   = 10000;          % K, Arrhenius temperature [Pethybridge et al 2013]

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
del_M   = 0.1889;% -

%del_length = 0; % 1 = Total Length | 0 = Standard Length
%if del_length == 1
%    del_M = del_M_t;
dirname = strcat('C:/Users/jflores/Desktop/DEB_outV4', 'L_w',num2str(talla));
mkdir(dirname);
subdir = dirname;
%else
%    del_M = del_M_s;
%    mkdir('DEB_out_s');
%    subdir = 'DEB_out_s';
%end

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
L_w    = zeros(n_iter,1);
f_vari = zeros(n_iter,1);

t(1)   = t_0;    % d,    Time vector initialization 
E(1)   = E_0;    % J,    Initial reserve
V(1)   = V_0;    % cm^3, Initial structure
E_R(1) = E_R0;   % J,    Reproduction buffer   ------ I put an indice to try to run the script
E_H(1) = E_H0;	 % J, Maturity
L_w(1) = (V(1)^(1/3))/del_M; % Talla inicial

for j = 1:size(temp,2)
    
    for k = 1:size(f_res,2)
        
        T = repmat(temp(j) + T_K, n_iter,1); % K, Temperature
               
        for i = 1:n_iter-1
            % Condicional para que luego de una talla determinada, f cambia
            % a una condicion diferente de 1
            if L_w(i) <= talla % Talla de la larva en cm antes de iniciar una condicion de 'f' diferente a 1
                f = 1;
            else
                f = f_res(k);
            end
            f_vari(i) = f;
            
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
        c_T   = exp(T_A/ T_ref - T_A ./ T(i));  % simple Arrhenius correction factor
        
		% Correction of physiology parameters for temperature :
		p_AmT  = c_T * p_Am;
		V_dotT = c_T * V_dot;
		p_MT   = c_T * p_M;
		k_JT   = c_T * k_J;
 
        % Scaled functional response
        % f = X(i) / (X(i) + K); % -, scaled functional response
        % f = 0.9999; % We assume food to satiation
		%% Metabolic acceleration – abj model
		if E_H(i) < E_Hb
			s_M = 1;
		elseif (E_Hb <= E_H(i) && (E_H(i) < E_Hj))
% 			s_M = V(i)^(1/3)/L_b;
            s_M = ((V(i)^(1/3))/del_M)/L_b; % Falta division del_M?
		else
			s_M = L_j / L_b;
		end
		
		% Only two parameters are accelerated by s_M, p_Am and v
		p_AmT  = s_M * p_AmT;
		V_dotT = s_M * V_dotT;
		
        % Fluxes , j/d
        if E_H < E_Hb % antes DEBstd : V(i) < V_b 
            p_A_flux = 0;
        else
            p_A_flux = f * p_AmT * V(i).^(2/3); % J/d assimilation rate 
        end

        p_M_flux = p_MT * V(i); % J/d , Volume-related somatic maintenance
		p_C_flux = (E(i) / V(i)) * (E_G * V_dotT * V(i)^(2/3) + p_M_flux) ./ (kap * E(i) ./ V(i) + E_G); % Energy for utilisation
		% p_C_flux = E(i) .* ( E_G * p_AmT / E_m .* V(i)^(-1/3) + p_MT) ./ ( kap .* E(i) ./ V(i) + E_G); % antes DEBstd
		p_G_flux = max(0, kap * p_C_flux - p_M_flux);% J/d growth
		% p_G_flux = kap * p_C_flux - p_M_flux;% antes DEBstd
		p_J_flux = k_JT * E_H(i)  ; % J/d, Maturity maintenance !!
		% p_J_flux = (1- kap) / kap * p_MT * min(V(i), V_p)  ; % antes DEBstd
		
		% if or( (p_G_flux < 0), ((1 - kap) * p_C_flux - p_J_flux <0)) % Antes DEBstd
        if or( kap * p_C_flux < p_M_flux, ((1 - kap) * p_C_flux - p_J_flux <0))
%             disp('starvation')
%             return
            E(i) = 0; % When the value of E = 0, I can identify starvation in the output.
        else
		
            dE = p_A_flux - p_C_flux ; % J/d, 
            dV = ((kap * p_C_flux) - p_M_flux) / E_G;       % cm^3/d,
			
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
		
        % Physical length
        L_w(i+1) = (V(i+1)^(1/3))/del_M ; % cm, Physical length
        
% 		% Spawning rule 1: Every 30 days
% 		if (i/360 - round(i/360) == 0)
% 			E_R(i+1) = 0; % E_R regresa a cero
% 		end
%                
% 		% Spawning rule 2: One spawning peak in September
% 		if (i == 3240 || i == 7560 || i == 11880 || i == 16200) % iterator indices for the September months
% 			E_R(i+1) = E_R(i) * 0.10; % El E_R es el 10% de la cantidad anterior
% 		end

		% Spawning rule 3: Two spawning peaks in September & March
		if (i+1 == 3240 || i+1 == 7560 || i+1 == 11880 || i+1 == 16200 || ... % iterator indices for the September months (el desove se hace entre i y i+1)
			i+1 == 1080 || i+1 == 5400 || i+1 == 9720  || i+1 == 14040)       % iterator indices for the March months
			F(i+1)   = E_R(i+1) * kap_R / E_0;
			E_R(i+1) = 0; % E_R regresa a cero
		end
        
        t(i+1)   = t(i) + dt ;
        
        t_vec = repmat(temp(j), n_iter,1);  % C, Temperature
        f_vec = repmat(f_res(k), n_iter,1); % f, Functional response
        
        end
        
        %% OBSERVABLE VARIABLES
    
%         % Physical length
%         L_w = (V.^(1/3))./del_M ; % cm, Physical length

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
		
        out_mat = table(t,E,V,E_H,E_R,Ww,L_w,F,t_vec,f_vec,f_vari,...
                        'VariableNames',...
                        {'t','E','V','E_H','E_R','Ww','L_w','F','temp','f','f_vari'});
        writetable(out_mat, strcat(subdir, '/DEB_out','T',num2str(temp(j)),'f',num2str(f_res(k)),'.txt'))
    end
end
