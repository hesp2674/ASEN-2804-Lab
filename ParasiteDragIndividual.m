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
Length_f = Design_Input.Length_f;
Rho_air = ATMOS.rho;
Nu_air = ATMOS.nu; 
Velocity = Design_Input.V_o;
Mach = Velocity/ATMOS.a;
A_max = Design_Input.Amax_f; 
Swet_f = Design_Input.Swet_f;
Q_f = Design_Input.Q_f;
Abase_f = Design_Input.Abase_f;
Sref_w = Design_Input.Sref_w;
Taper_w = Design_Input.Taper_w;
AR_w = Design_Input.AR_w;
X_thick_w = Design_Input.X_thick_w;
Thick_w = Design.Input.Thick_w;
Sweep_w = Design.Input.Sweep_w;
Q_w = Design.Input.Q_w;
Swet_w = Design.Input.Swet_f; 
%% Loop through different configurations
for n = 1:Count
% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////
    %% Fuselage Contribution To CDo
    k_surface = 00000634; % Value for smooth paint
    
    Re_f_L = (Velocity(n)*Length_f(n))/Nu_air; %Re for fuselage

    Re_f_cutoff = (Length_f(n)/k_surface); %Re model using surface roughness
    Re_f_cutoff = (Re_f_cutoff)^1.053;
    Re_f_cutoff = Re_f_cutoff*38.21;

    Re_f_eff = min([Re_f_cutoff, Re_f_L]); %Re to use in rest of calcs
    
    Cf_f = log10(Re_f_eff); %Leaving off mach correction
    Cf_f = Cf_f^2.58; 
    Cf_f = .455/Cf_f;
    
    F_int = (Length_f(n))/(sqrt((4*A_max)/pi));
    FF_f(n) = ((.9)+(5/(F_int)^1.5)+(F_int)/400); %Fuselage Form Factor

    CDo_f(n) = Cf_f*FF_f(n)*Swet_f(f)*Q_f; %Contribution of Fuselage to CDo

    %% Wing Contribution to CDo
    b_w = sqrt(AR_w(n)*Sref_w(n));
    Cr_w = (2*Sref_w(n))/(b_w*(1+Taper_w(n)));
    MAC_w = (2*Cr_w*(1+Taper_w(n)+(Taper_w(n))^2)/(3*(1+Taper_w(n))));
    Re_w = (Velocity(n)*MAC_w(n))/(Nu_air); %Wing Re
    Cf_w = .074/(Re_w^.2); %Wing Flat Plate Coef of Friction for Turbulent Flow
    FF_w(n) = (1 + ((.6/X_thick_w)*(Thick_w)) + (100*Thick_w^4))*(1.35*(Mach^.18)*(cosd(Sweep_w)^.28)); %Wing Form Factor

    CDo_w(n) = (Cf_w*FF_w(n)*Q_w*Swet_w)/(Sref_w); %Contribution of Wing to CDo

    %% Horizontal Tail #1 Contribution to CDo
    if Design_Input.Swet_h1(n)~=0 % If this component exists:
        Re_h1 = ; %Horz Tail Re
        Cf_h1 = ; %Flat Plate Coef of Friction for Turbulent Flow
        FF_h1(n) = ; %Horz Tail Form Factor
        CDo_h1(n) = ; %Contribution of Horz Tail 1 to CDo 
    end

    %% Horizontal Tail #2 Contribution to CDo
    if Design_Input.Swet_h2(n)~=0 % If this component exists:
        Re_h2 = ; %Horz Tail Re
        Cf_h2 = ; %Flat Plate Coef of Friction for Turbulent Flow
        FF_h2(n) = ; %Horz Tail Form Factor
        CDo_h2(n) = ; %Contribution of Horz Tail 2 to CDo 
    end

    %% Vertical Tail #1 Contribution to CDo
    if Design_Input.Swet_v1(n)~=0 % If this component exists:
        Re_v1 = ; %Horz Tail Re
        Cf_v1 = ; %Flat Plate Coef of Friction for Turbulent Flow
        FF_v1(n) = (1+(0.6/())); %Horz Tail Form Factor
        CDo_v1(n) = ; %Contribution of Vert Tail 1 to CDo 
    end

    %% Vertical Tail #2 Contribution to CDo
    if Design_Input.Swet_v2(n)~=0 % If this component exists:
        Re_v2 = ; %Horz Tail Re
        Cf_v2 = ; %Flat Plate Coef of Friction for Turbulent Flow
        FF_v2(n) = ; %Horz Tail Form Factor
        CDo_v2(n) = ; %Contribution of Vert Tail 2 to CDo
    end

    %% Misc. and L&P Contributions to CDo
    if Design_Input.Abase_f(n)~=0 % If this component exists:
        D_q_base = (Mach(n)-1.61)^2;
        D_q_base = (.419*D_q_base) + .139;
        D_q_base = D_q_base*Abase_f(n);

        %NEED CLARIFICATION ON WHAT Sref IS
        CDo_misc(n) = D_q_base*(1/Sref); %Contribution of base drag to CDo
    end
    
    %% Leakage and Proturbance Contribution to CDo
    PercentLeakage = ; % Your choice depending on how bad you think your fabrication will be 
    CDo_lp(n) = ; %Increase in parasite drag due to leakage and protuberance usually 3-15% of total CDo, but here just taking fuselage contribution
 
    %%Total Parasite Drag and Wetted Area
    Swet_tot(n) = ; %Total Wetted Area
    CDo(n) = ; %Total Parasite Drag Coefficient
% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
end
    
%% Oraganize into tables for output
Parasite_Data = table(CDo, CDo_f, CDo_w, CDo_h1, CDo_h2, CDo_v1, CDo_v2, CDo_misc, CDo_lp, Swet_tot);
FF_Table = table(FF_f, FF_w, FF_h1, FF_h2, FF_v1, FF_v2);

end

