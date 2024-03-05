% Henry

function [WingGeo_Data] = WingGeo(Design_Input,Count)
%% Wing Geometry Function Summary:
% This function creates a variable table that calculates the main wing's
% specific geometry from the user defined primary design inputs of Sref,
% Aspect Ratio, Taper Ratio, and Leading Edge Sweep which are provided from
% the Design Input spreadsheet file.  From these primary design parameters,
% this function will output the wingspan (b_w), root chord length (c_r),
% tip chord length (c_t), mean aerodynamic chord length (MAC_w), and the
% location of the aerodynamic center of the MAC in x and y coordinates
% (x_MAC_w and y_MAC_w) with the origin being the leading edge of the root
% chord.  These calculations are done for every configuration in the Design
% Input file (rows of the WingGeo_Data table are for each configuration).

%% Preallocate variables of interest
b_w = zeros(Count, 1); % Wingspan [m]
cr_w = zeros(Count, 1); % Root chord [m]
ct_w = zeros(Count, 1); % Tip chord [m]
MAC_w = zeros(Count, 1); % Mean aerodynamic chord [m]
y_MAC_w = zeros(Count, 1); % MAC y-location [m]
x_MAC_w = zeros(Count, 1); % MAC y-location [m]

%% Loop through different configurations
for n = 1:Count
% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////
    % Find Wing Span from AR and Sref
    b_w(n)= sqrt(Design_Input.AR_w(n) * Design_Input.Sref_w(n));
    
    %Find Wing Root and Tip Chord from S and Taper Ratio
    cr_w(n)= 2*Design_Input.Sref_w(n)/(b_w(n)*(1+Design_Input.Taper_w(n)));
    ct_w(n)= cr_w(n)*Design_Input.Taper_w(n);
    
    %Calculate Wing Mean Aerodynamic Chord (MAC)
    MAC_w(n)= (2/3)*cr_w(n)*((1+Design_Input.Taper_w(n)+Design_Input.Taper_w(n)^2)/(1+Design_Input.Taper_w(n)));

    %Find x and y location for wing MAC from leading edge of root chord
    y_MAC_w(n)=  (b_w(n)/6)*((1+2*Design_Input.Taper_w(n))/(1+Design_Input.Taper_w(n)));
    x_MAC_w(n)= y_MAC_w(n)*tan(Design_Input.Sweep_w(n)) + 0.25*MAC_w(n);
% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
end

%% Oraganize into tables for output
WingGeo_Data = table(b_w, cr_w, ct_w, MAC_w, y_MAC_w, x_MAC_w);   

end
