function kin = define_CRX()
% https://crx.fanucamerica.com/wp-content/uploads/2022/10/CRX-10iA-L-data-sheet.pdf
    ex = [1;0;0];
    ey = [0;1;0];
    ez = [0;0;1];
    zv = [0;0;0];
    kin.H = [ez ex ex ey ex ey];
    % kin.P = [245*ez zv 710*ez zv 540*ey+150*ex zv 160*ey]/1000;
    kin.P = [zv zv 710*ez zv 540*ey+150*ex zv zv]/1000; % p_01 = p_6T = 0
    % syms L1 L2 L3 real
    % kin.P = [zv zv L1*ez zv L2*ey zv L3*ex zv zv]; % p_01 = p_6T = 0
    kin.joint_type = zeros(1,6);
end