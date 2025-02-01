%% A function to calculate the effective diameter (circular diameter with equivalent resistance to the rectangular cross-section) 

function D_eff = calc_effective_diameter(height, aspect_ratio)

f = 0;

for kk=1:50
    f= f + 1/(2*kk-1)^5*(cosh((2*kk-1)*pi*aspect_ratio) - 1)/sinh((2*kk-1)*pi*aspect_ratio);
end

D_eff = height*(128/pi*(aspect_ratio/12 - 16/pi^5*f))^0.25;

% Just calculate the hydraulic diameter instead
% D_eff = 2 * height * aspect_ratio / (1 + aspect_ratio);
disp(['D_eff = ', num2str(D_eff)]);

end
