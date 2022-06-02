%% pyruvate + H+ -> acetaldehyde + "carbon dioxide" 
%PDC1 consensus
v_PDC=p.PDC_ExprsCor.*((p.PDC1_kcat.*(f.PDC1).*(PYR./p.PDC1_Kpyr).^p.PDC1_hill)./...
    (1+(PYR./p.PDC1_Kpyr).^p.PDC1_hill+PI./p.PDC1_Kpi));
%     (1+(PYR./p.PDC1_Kpyr).^p.PDC1_hill));%+PI./p.PDC1_Kpi));

%% acetaldehyde + H+ + NADH = ethanol + NAD+;  
    %ADH teusink
    v_ADH=-p.ADH_ExprsCor.*(p.ADH.VmADH./(p.ADH.KiADHNAD.*p.ADH.KmADHETOH).*(NAD.*ETOH-NADH.*ACE./p.ADH.KeqADH)./...
        (1+NAD./p.ADH.KiADHNAD+p.ADH.KmADHNAD.*ETOH./(p.ADH.KiADHNAD*p.ADH.KmADHETOH)+p.ADH.KmADHNADH.*ACE./(p.ADH.KiADHNADH.*p.ADH.KmADHACE)...
    +NADH./p.ADH.KiADHNADH+NAD.*ETOH./(p.ADH.KiADHNAD.*p.ADH.KmADHETOH)+p.ADH.KmADHNADH.*NAD.*ACE./(p.ADH.KiADHNAD.*p.ADH.KiADHNADH.*p.ADH.KmADHACE)...
    +p.ADH.KmADHNAD.*ETOH.*NADH./(p.ADH.KiADHNAD.*p.ADH.KmADHETOH.*p.ADH.KiADHNADH)+NADH.*ACE./(p.ADH.KiADHNADH.*p.ADH.KmADHACE)+...
    NAD.*ETOH.*ACE./(p.ADH.KiADHNAD*p.ADH.KmADHETOH.*p.ADH.KiADHACE)+ETOH.*NADH.*ACE./(p.ADH.KiADHETOH*p.ADH.KiADHNADH*p.ADH.KmADHACE)));

%% "sn-glycerol 3-phosphate" + water -> glycerol + phosphate; 
%HOR2 consensus
%oldFormula % v_HOR2=(((p.HOR2_kcat.*f.HOR2)./p.HOR2_Kglyc3p.*(1+p.HOR2_Kpi./PI).*GLYC3P)./...
%oldFormula %     (1+GLYC3P./p.HOR2_Kglyc3p));
v_HOR2=(((p.HOR2_kcat.*f.HOR2)./p.HOR2_Kglyc3p.*GLYC3P)./...
    ((1+PI./p.HOR2_Kpi).*(1+GLYC3P./p.HOR2_Kglyc3p)));

%% Branches

v_ACE=p.VmaxACE*ACE./(ACE+p.KmACE);

ETOHe=f.ETOH_e;
v_ETOHtransport=p.kETOHtransport*(ETOH-ETOHe);

GLYCEROLe=f.GLYCEROL_e;
v_GLYCEROLtransport=p.GlycerolTransport.*(GLYCEROL-GLYCEROLe);

v_PDH=p.PDH_Vmax.*(PYR.^p.PDH_n)./(p.PDH_K50.^p.PDH_n+PYR.^p.PDH_n);

v_mitoNADH=p.mitoNADHVmax*(NADH/(NADH+p.mitoNADHKm));