function [Parasite_Data,FF_Table] = ...
    ParasiteDrag(Design_Input,Airfoil,WingGeo_Data,ATMOS,Count)
%%  Parasite Drag Summary
% This function performs the Raymer Component Drag Buildup Method to
% determine total and component parasite drag coefficients using the
% aircraft configuration geometry information from the Design Input,
% Airfoil, and WingGeo_Data tables.  Additionally, it leverages the 
% standard atmosphere properties from the ATMOS table. The output table 
% for this function includes the total parasite drag coefficient (CDo), and
% a breakdown of the contribution to CDo for the fuselage (f), wing (w), 
% horizontal stabilizers (h1, h2), vertical stablizers (v1,v2), misc base 
% drag (misc), and leakage and proturbance drag (lp).  This function also
% outputs the total wetted area for the entire aircraft (Swet_total). 
% Finally, the FF_Table consolidates the form factors (FF) which account
% for zero lift pressure drag for each component for evaluation.
% Note that no FF should be < 1.0.  All these calculation are done for each
% configruation in the Design Input spreadsheet.

%% Outputs:
%
% Parasite Data:
%   Table containing total parasite drag and a component breakdown of
%   contributions to that total (columns), for each input case (rows)
%
% FF_Table:
%   Table containing total form factor and a component breakdown of
%   contributions to that total (columns), for each input case (rows)

%% Preallocate variables of interest
CDo = zeros(Count, 1); % Total parasite drag coefficient
CDo_f = zeros(Count, 1); % Fusalage parasite drag coefficient contribution
CDo_w = zeros(Count, 1); % Wing parasite drag coefficient contribution
CDo_h1 = zeros(Count, 1); % Horizontal stabilizer 1 parasite drag coefficient contribution
CDo_h2 = zeros(Count, 1); % Horizontal stabilizer 2 parasite drag coefficient contribution
CDo_v1 = zeros(Count, 1); % Vertical stabilizer 1 parasite drag coefficient contribution
CDo_v2 = zeros(Count, 1); % Vertical stabilizer 2 parasite drag coefficient contribution
CDo_misc = zeros(Count, 1); % Misc. parasite drag coefficient contribution
CDo_lp = zeros(Count, 1); % Leakage and purturbance parasite drag coefficient contribution
Swet_tot = zeros(Count, 1); % Total wetted area [m^2]

FF_f = zeros(Count, 1); % Fusalage form factor
FF_w = zeros(Count, 1); % Wing form factor
FF_h1 = zeros(Count, 1); % Horizontal stabilizer 1 form factor
FF_h2 = zeros(Count, 1); % Horizontal stabilizer 2 form factor
FF_v1 = zeros(Count, 1); % Vertical stabilizer 1 form factor
FF_v2 = zeros(Count, 1); % Vertical stabilizer 2 form factor



%% Variables Pulled From Tables


%% Loop through different configurations
for n = 1:Count
% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////

    Length_f = Design_Input.Length_f(n);
Rho_air = ATMOS.rho(n);
Nu_air = ATMOS.nu(n); 
Velocity = Design_Input.V_o(n);
Mach = Velocity/ATMOS.a(n);
A_max = Design_Input.Amax_f(n); 
Swet_f = Design_Input.Swet_f(n);
Q_f = Design_Input.Q_f(n);
Abase_f = Design_Input.Abase_f(n);
v1_x_c= Airfoil.X_thick_v1(n);
v1_t_c= Airfoil.Thick_v1(n);
v1_sweep= Design_Input.Sweep_v1;(n)
v1_Swet= Design_Input.Swet_v1(n);
v2_x_c= Airfoil.X_thick_v2(n);
v2_t_c= Airfoil.Thick_v2(n);
v2_sweep= Design_Input.Sweep_v2(n);
v2_Swet= Design_Input.Swet_v2(n);
Surface_finish=0.635*(10^(-5)); %Changes k value

Re_h = 5*(10^(5));
MAC_h1=Design_Input.MAC_h1(n);
MAC_h2=Design_Input.MAC_h2(n);
MAC_v1=Design_Input.MAC_v1(n);
MAC_v2=Design_Input.MAC_v2(n);

h1_x_c= Airfoil.X_thick_h1(n);
h1_t_c= Airfoil.Thick_h1(n);
h1_taper= Design_Input.Taper_h1(n);
h1_sweep = Design_Input.Sweep_h1(n);
    % h1_Q = Design_Input.Q_h1;
    % h1_wet = Design_Input.Swet.h1;
    % h_sref = Design_Input.Sref_w;
    % h2_x_c= Airfoil.X_thick_h2;
    % h2_t_c= Airfoil.Thick_h2;
    % h2_taper= Design_Input.Taper_h2;
    % h2_sweep = Design_Input.Sweep_h2;
    % h2_Q = Design_Input.Q_h2;
    % h2_wet = Design_Input.Swet.h2;
Sref_w = Design_Input.Sref_w(n);
Taper_w = Design_Input.Taper_w(n);
AR_w = Design_Input.AR_w(n);
X_thick_w = Airfoil.X_thick_w(n);
Thick_w = Airfoil.Thick_w(n);
Sweep_w = Design_Input.Sweep_w(n);
Q_w = Design_Input.Q_w(n);
Swet_w = Design_Input.Swet_w(n); 
    %% Fuselage Contribution To CDo
    k_surface = .00000634; % Value for smooth paint(.00000634) Ref paint on aluminum .00001015

    Re_f_L = (Velocity*Length_f)/Nu_air; %Re for fuselage %note possibly multiply air rho

    Re_f_cutoff = (Length_f/k_surface); %Re model using surface roughness
    Re_f_cutoff = (Re_f_cutoff)^1.053;
    Re_f_cutoff = Re_f_cutoff*38.21;

    Re_f_eff = min([Re_f_cutoff, Re_f_L]); %Re to use in rest of calcs
    
    Cf_f = log10(Re_f_eff); %Leaving off mach correction
    Cf_f = Cf_f^2.58; 
    Cf_f = .455/Cf_f;
    
    F_int = (Length_f )/(sqrt((4*A_max )/pi));
    FF_f(n)  = ((.9)+(5/(F_int)^1.5)+(F_int)/400); %Fuselage Form Factor

    CDo_f(n)  = Cf_f*FF_f(n) *Swet_f *Q_f /Sref_w ; %Contribution of Fuselage to CDo

    %% Wing Contribution to CDo
    b_w = sqrt(AR_w *Sref_w );
    Cr_w = (2*Sref_w )/(b_w *(1+Taper_w ));
    MAC_w = (2*Cr_w *(1+Taper_w +(Taper_w )^2)/(3*(1+Taper_w )));
    Re_w = (Velocity *MAC_w )/(Nu_air ); %Wing Re
    Cf_w = .074/(Re_w^0.2); %Wing Flat Plate Coef of Friction for Turbulent Flow
    FF_w(n)  = (1 + ((.6/X_thick_w )*(Thick_w )) + (100*Thick_w ^4))*(1.35*(Mach^.18)*(cosd(Sweep_w )^.28)); %Wing Form Factor
    if FF_w(n) <1
        FF_w(n) =1;
    end

    CDo_w(n)  = (Cf_w*FF_w(n) *Q_w *Swet_w )/(Sref_w ); %Contribution of Wing to CDo
    %% Horizontal tail values
    h1_x_c= Airfoil.X_thick_h1(n);
    h1_t_c= Airfoil.Thick_h1(n);
    h1_taper= Design_Input.Taper_h1(n);
    h1_sweep = Design_Input.Sweep_h1(n);
    h1_Q = Design_Input.Q_h1(n);
    v1_Q = Design_Input.Q_v1(n);
    v2_Q = Design_Input.Q_v2(n);

    h1_wet = Design_Input.Swet_h1(n);
    h_sref = Design_Input.Sref_w(n);
    h2_x_c= Airfoil.X_thick_h2(n);
    h2_t_c= Airfoil.Thick_h2(n);
    h2_taper= Design_Input.Taper_h2(n);
    h2_sweep = Design_Input.Sweep_h2(n);
    h2_Q = Design_Input.Q_h2(n);
    h2_wet = Design_Input.Swet_h2(n);
    v2_wet = Design_Input.Swet_v2(n);
    v1_wet = Design_Input.Swet_v1(n);


    %% Horizontal Tail #1 Contribution to CDo
    if Design_Input.Swet_h1(n)~=0 % If this component exists:
        Re_h1 = (Velocity *MAC_h1 )/Nu_air ;
        Cf_h1 = 0.074/(Re_h1 ^.2); %Flat Plate Coef of Friction for Turbulent Flow
        FF_h1(n)  = (1+(.6/h1_x_c )*h1_t_c +100*(h1_t_c )^4)*(1.35*(Mach ^.18))*(cos(h1_sweep ))^.28; %Horz Tail Form Factor
        if FF_h1(n) <1
            FF_h1(n) =1;
        end

        CDo_h1(n)  = (Cf_h1 *FF_h1(n) *h1_Q *h1_wet )/Sref_w ; %Contribution of Horz Tail 1 to CDo
    end

    %% Horizontal Tail #2 Contribution to CDo
    if Design_Input.Swet_h2 ~=0 % If this component exists:
        Re_h2 = (Velocity *MAC_h2 )/Nu_air ;
        Cf_h2 = 0.074/(Re_h2^.2); %Flat Plate Coef of Friction for Turbulent Flow
        FF_h2(n)  = (1+(.6/h2_x_c)*h2_t_c+100*(h2_t_c)^4)*(1.35*(Mach ^.18))*(cos(h2_sweep))^.28; %Horz Tail Form Factor
        if FF_h2(n) <1
            FF_h2(n) =1;
        end
        CDo_h2(n)  = Cf_h2*FF_h2(n)*h2_Q*h2_wet/Sref_w ; %Contribution of Horz Tail 2 to CDo 
    end

    %% Vertical Tail #1 Contribution to CDo
    if Design_Input.Swet_v1 ~=0 % If this component exists:
        Re_v1 = (Velocity *MAC_v1 )/Nu_air ;
        Cf_v1 = 0.074/(Re_v1)^0.2; %Flat Plate Coef of Friction for Turbulent Flow
        FF_v1(n)  = (1+(0.6/(v1_x_c ))*(v1_t_c )+100*(v1_t_c )^4)*(1.35*(Mach ^0.18)*cos(v1_sweep )^0.28); %Horz Tail Form Factor
        if FF_v1(n) <1
            FF_v1(n) =1;
        end

        CDo_v1(n)  = Cf_v1*FF_v1(n) *v1_Q*v1_Swet /Sref_w ; %Contribution of Vert Tail 1 to CDo 
    end

    %% Vertical Tail #2 Contribution to CDo
    if Design_Input.Swet_v2 ~=0 % If this component exists:
        Re_v2 = (Velocity *MAC_v2 )/Nu_air ;
        Cf_v2 = 0.074/(Re_v2)^0.2; %Flat Plate Coef of Friction for Turbulent Flow
        FF_v2(n)  = (1+(0.6/(v2_x_c ))*(v2_t_c )+100*(v2_t_c )^4)*(1.35*(Mach ^0.18)*cos(v2_sweep )^0.28); %Horz Tail Form Factor
        if FF_v2(n) <1
            FF_v2(n) =1;
        end
        CDo_v2(n)  = Cf_v2*FF_v2(n) *v2_Q*v2_Swet /Sref_w ; %Contribution of Vert Tail 1 to CDo
    end

    %% Misc. and L&P Contributions to CDo
    if Design_Input.Abase_f(n)~=0 % If this component exists:
        D_q_base = (Mach -0.161)^2;
        D_q_base = (.419*D_q_base) + .139;
        D_q_base = D_q_base*Abase_f ;

        %NEED CLARIFICATION ON WHAT Sref IS
        CDo_misc(n)  = D_q_base*(1/Sref_w ); %Contribution of base drag to CDo
    end
    
    %% Leakage and Proturbance Contribution to CDo
    PercentLeakage = .1; % Your choice depending on how bad you think your fabrication will be 
    CDo_lp(n) = PercentLeakage*(CDo_f(n) +CDo_v2(n) +CDo_v1(n) +CDo_h2(n) +CDo_h1(n) +CDo_w(n) ); %Increase in parasite drag due to leakage and protuberance usually 3-15% of total CDo, but here just taking fuselage contribution
 
    %%Total Parasite Drag and Wetted Area
    Swet_tot(n)  = h1_wet + h2_wet  + Swet_w  + Swet_f + v1_wet + v2_wet ; %Total Wetted Area
    CDo(n) = CDo_f(n) +CDo_misc(n)+CDo_v2(n)+CDo_v1(n)+CDo_h2(n)+CDo_h1(n)+CDo_w(n) + CDo_lp(n); %Total Parasite Drag Coefficient
% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
end
    
%% Oraganize into tables for output
Parasite_Data = table(CDo, CDo_f, CDo_w, CDo_h1, CDo_h2, CDo_v1, CDo_v2, CDo_misc, CDo_lp, Swet_tot);
FF_Table = table(FF_f, FF_w, FF_h1, FF_h2, FF_v1, FF_v2);

end

