function [solution,fluxes,min,max,optimal]=FVA_iDunsB(model,rxns_input)
%Model iDuns loaded ready for TestSolutionChlorella
%Jenna Armstrong & Cristal Zuniga
%05.13.2024
%
%
%clear all
%initCobraToolbox    

model=model;
DW = 13*10^(-12); % avg. dry weight of log phase chlamy cell = 13 pg (This study)
CPerStarch300 = 1800; % derived from starch300 chemical formula
starchPAdD=0.042; % Measured in dark g/gDW d
starchPAd=0.113; %Measured in light g/gDW d
starchDegAerLight = starchPAd/(CPerStarch300*24);% estimated fron the rates calculated with experimental data mmol/g/h
starchDegAnLight = (3/2)*starchDegAerLight; % approx. SS rate of aerobic starch degradation in light 3/2 of anaerobic rate in green algae (Nakamura and Miyachi, 1982)
starchDegAnDark = starchPAdD/(CPerStarch300*24); % estimated fron the rates calculated with experimental data mmol/g/h
starchDegAerDark = (3/2)*starchDegAnDark; % approx. SS rate of aerobic starch degradation in light 3/2 of anaerobic rate in green algae (Nakamura and Miyachi, 1982)
dimensionalConversion = 3.836473679; % from emitted microE/m^2/s to incident mmol/gDW/hr
effectiveConversion = 0.06; % from incident mmol/gDw/hr to effective mmol/gDw/hr [Hiramaya 1995]

rxn_List1 = {'GRDPth';'GRTTh';'GPPS';'POAT';'POATh';'PHEOXh';'AKG_na_th';...
            'AKGMALth';'AKG_na_tm';'AKGCITtm';'AKGICITtm';'AKGMALtm';...
            'OAAAKGtm';'ALATA_Lm';'ASPTAh';'ACOTA';'ALATA_L';'ASPTAm';...
            'CYSTAm';'CYSATm';'GLUS';'HSTPT';'LDAPAT';'UNK3';'VALTAim';...
            'POAT';'POATh';'TYRTAh';'ASPTA';'ICDHyr';'ICDHhr';'ILETA2';...
            'BCTA_glu';'BCTA_val_h';'ACAROtu';'CAROtu';'GCAROtu';'BCAROH';...
            'CHYA2';'LCYA';'LCYB';'LCYD';'LCYG';'PDS2';'ZDS';'CBCI'};

rxn_List2 = {'LLOX13';'FA70ACPHi';'FA70h';'FACOAL70';'TAGAH13013070';'TAGAH1307070';'TAGAH13021070';'TAGAH15015070';'TAGAH1507070';...
    'TAGAH15013070';'TAGAH15021070';'TAGAH21021070';'TAGAH2107070';'TAGAH17013070';'TAGAH1707070';'TAGAH17021070';'TAGAH17015070';...
    'TAGAH17017070';'TAGAH21017070';'PHEth';'PHETRSm';'POATm';'EX_phe-L_e';'GGDPtu';'GGDR';'GGPS_h';'PSY';'CHLASG';'CHLBSG';'DEDOLDPS';...
    'FT';'NONDPS2h';'HGGTh';'FA140ACPHi';'AACPS1';'FAH140';'FAC140t_c';'FACOAL140';'FA150ACPHi';'FA150h';'FACOAL150';'TAGAH130130150';...
    'TAGAH13070150';'TAGAH130210150';'TAGAH150150150';'TAGAH15070150';'TAGAH150130150';'TAGAH150210150';'TAGAH210210150';'TAGAH21070150';...
    'TAGAH170130150';'TAGAH17070150';'TAGAH170210150';'TAGAH170150150';'TAGAH170170150';'TAGAH210170150';'FA160h';'FA180h';'FACOAE160';...
    'FACOAE180';'LCYE';'3OAS120';'3OAS140';'ACP1819ZD9DS';'EAR100x';'EAR100y';'EAR120x';'EAR120y';'EAR180x';'FA100ACPHi';'FA120ACPHi';...
    'FA160ACPHi';'FA180ACPHi';'ACACT5m';'ACACT6m';'ACOAR4m';'ACOAR5m';'AACPS3';'AACPS6';'ACACT5';'ACACT6';'ACOA100OR';'ACOA120OR';...
    'FACOAL160';'FAH120';'FAH160';'ASQDCADS1829Z12Z160';'ASQDCADS1839Z12Z15Z160';'FACOAL180';'G3PAT180h';'MAGAH160';'MAGAH180';...
    'SQDGATCA18111Z160';'SQDGATCA1819Z160';'TAGAH16018111Z160';'TAGAH16018111Z180';'TAGAH1601819Z160';'TAGAH1601819Z180';'TAGAH1801819Z160';...
    'TAGAH1801819Z180';'TAGAH18111Z18111Z160';'TAGAH18111Z18111Z180';'TAGAH18111Z1819Z160';'TAGAH18111Z1819Z180';'TAGAH1819Z18111Z160';...
    'TAGAH1819Z18111Z180';'TAGAH1819Z1819Z160';'TAGAH1819Z1819Z180';'LPLPS1AGPE180';'ADSL1r';'ADSS';'ADSSh';'SPHPL';'RETPALMIH';...
    'EX_dca_e';'DCATDc';'TAGAH1835Z9Z12Z1835Z9Z12Z160';'TAGAH1835Z9Z12Z1835Z9Z12Z180';'TAGAH160160180';'CHLASE';'LNLCth';'FACOAE1829Z12Z';...
    'DESAT18_9';'FAH1829Z12Z';'FACOAL1821';'1AGPEAT1801829Z12Z';'1AGPEAT18111Z1829Z12Z';'1AGPEAT1819Z1829Z12Z';'LPLPS1AGPE1829Z12Z';...
    'PLPSA21801829Z12Z';'PLPSA218111Z1829Z12Z';'PLPSA21819Z1829Z12Z';'LNLCCOADS';'1AGPEAT1829Z12Z1829Z12Z';'MGDGH2h';'3OAS180';...
    'ACOADAGAT1835Z9Z12Z1835Z9Z12Z1829Z12Z';'TAGAH1835Z9Z12Z1835Z9Z12Z1829Z12Z';'G3PAT182';'AGPATCOA1829Z12Z1835Z9Z12Z';'ACP1619ZD9DS';...
    'PLDAGAT1829Z12Z1835Z9Z12Z1835Z9Z12Z2';'PLDAGAT18111Z1619Z1835Z9Z12Z2';'PLDAGAT1819Z1619Z1835Z9Z12Z2';'3OAS70';'3OAR70';'3HAD70';...
    'EAR70x';'EAR70y';'3OAS90';'AGPAT13070h';'AGPAT15070h';'AGPAT21070h';'AGPAT17070h';'PPNDH';'3OAS160';'EAR140x';'EAR140y';'AGPAT180140h';...
    'AGPAT1819Z140h';'3OAS150';'3OAR150';'3HAD150';'EAR150x';'EAR150y';'3OAS170';'G3PAT150h';'AGPAT150h';'AGPAT170150h';'3HAD160';'3OAR160';...
    'EAR160x';'EAR160y';'AGPAT160h';'AGPAT18111Z160h';'AGPAT1819Z160h';'G3PAT160h';'CHLAtu';'CHLADMT';'CHLAE';'CHLASP';'GGCHLDAR'};

rxn_List3 = {'ADNtm';'AHCi';'AHCm';'ADNDA';'ADNK';'ADNKm';'ADNUC';'AMP5N';'SEAHCYSHYD';'ARG-Ltm';'ORNtm';'ARGDI';'ARGTRS';'NOS';'NOS1';...
    'ARGSL';'ARGDCm';'ARGN';'ACAROtu';'CAROtu';'GCAROtu';'BCAROH';'BCAROKE';'CHYA1';'CHYA2';'LCYA';'LCYB';'LCYD';'LCYE';'LCYG';'PDS2';...
    'ZDS';'CAROMO';'CBCI';'H2O2tx';'H2O2tm';'ASPO1';'ASPOm';'PAO';'SPDH';'CITALALDOR';'CITALOR';'DDMCITALOR';'DMCITALOR';'GTHP';'G3PO';...
    'AMO';'GLYCTO1p';'AO';'MAOX';'CCP2m';'PTOR';'CPPPGO';'PPPGO';'URO';'DHORDi';'DHORDim';'SULO';'SULOm';'PDX5POi';'PYAM5PO';'PYDXDH';...
    'PYDXNO';'PYDXO';'GLYO1x';'PHEOXh';'TYROXh';'ASCBPOX_u';'GALOc';'ABSALDO';'FA160h';'FACOAE160';'FA160ACPHi';'AACPS3';'FACOAL160';...
    'FAH160';'MAGAH160';'TAGAH16018111Z160';'TAGAH1601819Z160';'TAGAH1801819Z160';'TAGAH18111Z18111Z160';'TAGAH18111Z1819Z160';...
    'TAGAH1819Z18111Z160';'TAGAH1819Z1819Z160';'RETPALMIH';'TAGAH1835Z9Z12Z1835Z9Z12Z160';'HISth';'HISt';'AABHH';'HIDP';'HDC';'HDH';...
    'HDHh';'HISTRS';'HISDr';'HXANt';'HXPRT';'INSH';'XANDH';'ILEth';'ILEtm';'ILETRS';'ILETA2';'ILETA2h';'6MPURPRT';'6TINS5MPOR';'I5NT';...
    'INSt2';'LEUth';'LEUtm';'LEUTL';'BCTA_glu';'BCTA_glu_h';'LIDOt';'LIDOAMH';'LIDODM';'LIDOMO';'LNLCth';'FACOAE1829Z12Z';'DESAT18_9';...
    'FAH1829Z12Z';'FACOAL1821';'1AGPEAT1801829Z12Z';'1AGPEAT18111Z1829Z12Z';'1AGPEAT1819Z1829Z12Z';'LPLPS1AGPE1829Z12Z';...
    'PLPSA21801829Z12Z';'PLPSA218111Z1829Z12Z';'PLPSA21819Z1829Z12Z';'LLOX13';'LNLCCOADS';'1AGPEAT1829Z12Z1829Z12Z';'MGDGH2h';...
    'ACOADAGAT1835Z9Z12Z1835Z9Z12Z1829Z12Z';'TAGAH1835Z9Z12Z1835Z9Z12Z1829Z12Z';'G3PAT182';'AGPATCOA1829Z12Z1835Z9Z12Z';...
    'PLDAGAT1829Z12Z1835Z9Z12Z1835Z9Z12Z2';'PLDAGAT18111Z1619Z1835Z9Z12Z2';'PLDAGAT1819Z1619Z1835Z9Z12Z2';'FACOAE1839Z12Z15Z';'LNLNCACOAL';...
    'LNLNCALO';'MGDGH3h';'LNLNCAth';'DAPDC';'LYSTRSm';'SACCD2';'DAGAH1801819Z';'MAGAH180';'METtm';'METAT';'METATm';'METTRS';'UNK3';'MS';...
    'METSm';'AMPMS_h';'BTS4h';'LIPOASh';'LIPOAS2';'MG2th';'MG2t';'MG2tu';'CHLADMT';'CHLDADMT';'MPML';'CYSTRS';'MM3';'BDMT';'DOLASNT';...
    'DOLPGT1';'DOLPGT2';'DOLPGT3';'DOLPMT1';'DOLPMT2';'G12MT1';'G12MT2';'G12MT3';'G12MT4';'G13MT';'G16MT';'GLCNACPT';'GLCNACT';'MG1A';...
    'MG1B';'MG2A';'MG2B';'MG3A';'MG3B';'MM1';'MM2';'NNAM';'NNDMBRT';'DMT_c';'FA180h';'FACOAE180';'ACP1819ZD9DS';'EAR180x';'FA180ACPHi';...
    'AACPS6';'FACOAL180';'G3PAT180h';'TAGAH16018111Z180';'TAGAH1601819Z180';'TAGAH1801819Z180';'TAGAH18111Z18111Z180';'TAGAH18111Z1819Z180';...
    'TAGAH1819Z18111Z180';'TAGAH1819Z1819Z180';'LPLPS1AGPE180';'TAGAH1835Z9Z12Z1835Z9Z12Z180';'TAGAH160160180';'FA1819Zh';'FACOAE1819Z';...
    'FA1819ZACPH';'AACP1819ZS';'FAH1819Z';'DAGAH1601819Z';'DAGAH18111Z1819Z';'DAGAH1819Z1819Z';'FACOAL181';'MAGAH1819Z';'TAGAH16018111Z1819Z';...
    'TAGAH1601819Z1819Z';'TAGAH1801819Z1819Z';'TAGAH18111Z18111Z1819Z';'TAGAH18111Z1819Z1819Z';'TAGAH1819Z18111Z1819Z';'TAGAH1819Z1819Z1819Z';...
    'LPLPS1AGPE1819Z';'PLPSA21801819Z';'PLPSA218111Z1819Z';'PLPSA21819Z1819Z';'MGDGH1h';'TAGAH1835Z9Z12Z1835Z9Z12Z1819Z';'PHEth';'PHETRSm';...
    'POAT';'POATh';'POATm';'PNTORtm';'PBAL';'APPT';'PBALm';'PDXPP';'PYDXNK';'PYDXOR';'RIBFLVth';'ARPT';'DBRC';'ACP1_FMN';'RIBFS';'TRPth';...
    'SERH';'TRPS2h';'TRPTL';'TRPO2';'TRPt';'FA140ACPHi';'AACPS1';'FAH140';'FAC140t_c';'FACOAL140';'TYRt2h';'TYRt2m';'TYRTA';'TYRTAh';...
    'TYRTAm';'TYRTRS';'PAPSTYRST-2';'TAT';'TYRt';'UREAtm_1';'URCBm';'ATAH';'URCB';'UREA';'VALth';'VALtm';'VALTAim';'VALTL';'VALTLh';'BCTA_val';...
    'BCTA_val_h';'VALt4';'ANXANASCOR';'ANXANOR';'BCRPTXANH';'CXHY';'LUTH';'NEOXANS';'VIOXANOR';'ZAXANOR';'ZHY';'6TXAN5MPAML';'GUAD';'XAND';...
    'XPPT';'XTSNRH';'CVIOXANS';'CVIOXAND';'XANXND';'VIOXANtu';'VIOXANth';'APLh';'APLm';'PDHam1hi';'PDHam1mi';'PDHam2hi';'PDHam2mi';'PDCh';...
    'PDCm';'ACAS_2ahbut';'MOD';'MOD_2mbdhl';'MOD_2mhop';'MOD_3mhtpp';'MOD_3mop';'MOD_4mop';'TMPPP2h';'TMPPH_h';'TMDPK';'TMDPPK';'TMNh';...
    'TMPPPh';'THMt';'ACCOAC_1';'BTNC';'BACCLh';'BTNPLh';'BIOTINt_h';'CHRPLh';'HBZOPTh';'CRETINOLOR1';'CRETINOLOR2';'CRETINOLPMTACT';'RETACIH';...
    'RETINOLACACT';'RETINOLSAT';'TRETINOLOR1';'TRETINOLOR2';'TRETINOLPMTACT'};

%rxn_List = rxn_List2;
%rxn_List = rxn_List3;
rxn_List = rxns_input;
%Green->photo
%Yellow->mix
%Orange->hetero
react_flux = model.rxns;
fva_min = rxn_List;
fva_max = rxn_List;
fba_op = rxn_List;

model = changeRxnBounds(model,{...
    'PRISM_solar_litho',...
    'PRISM_design_growth',...
    'PRISM_solar_exo',...
    'PRISM_fluorescent_warm_18W'...
    'PRISM_incandescent_60W',...
    'PRISM_fluorescent_cool_215W',...
    'PRISM_metal_halide',...
    'PRISM_high_pressure_sodium',...
    'PRISM_growth_room',...
    'PRISM_white_LED',...
    'PRISM_red_LED_array_653nm',...
    'PRISM_red_LED_674nm',...
    },[646.066560000000,417.591200000000,15.9417600000000,8.09577000000000,44.6333200000000,17.5256000000000,36.1503700000000,58.4672200000000,4.59151000000000,96.6281100000000,3.65300000000000,51.8841087821545],'b');

%% Photoautotrophy Green
modelLna = model;
% The single PRISM reaction being used has to be commented-out below.
modelLna = changeRxnBounds(modelLna,{...
   'PRISM_solar_litho',...
    'PRISM_design_growth',...
   'PRISM_solar_exo',...
    'PRISM_fluorescent_warm_18W'...
...    'PRISM_incandescent_60W',...
    'PRISM_fluorescent_cool_215W',...
    'PRISM_metal_halide',...
    'PRISM_high_pressure_sodium',...
    'PRISM_growth_room',...
    'PRISM_white_LED',...
    'PRISM_red_LED_array_653nm',...
    'PRISM_red_LED_674nm',...
},0,'b');
modelLna = changeRxnBounds(modelLna,{'EX_o2_e'},-1000,'l');
modelLna = changeRxnBounds(modelLna,{'EX_hco3_e'},-13.54,'l');
modelLna = changeRxnBounds(modelLna,{'EX_ac_e'},0,'l');
modelLna = changeRxnBounds(modelLna,{'EX_no3_e'},-10,'l');
modelLna = changeRxnBounds(modelLna,{'EX_nh4_e'},0,'l');
modelLna = changeRxnBounds(modelLna,{'EX_starch_h'},0,'b');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGRA',starchDegAerLight,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGR2A',0,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGRB',starchDegAerLight,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGR2B',0,'u');
modelLna = changeRxnBounds(modelLna,{'PCHLDR'},0,'b'); % Crystal Structure of the Nitrogenase-like Dark Operative Protochlorophyllide Oxidoreductase Catalytic Complex [Brocker 2010]
modelLna = changeRxnBounds(modelLna,{'G6PADHh','G6PBDHh'},0,'b'); % Purication of chloroplast G6PDH had not been successful, partly because this isoform is inactivated by light because of redox modulation by thioredoxin, and partly because it easily aggregates unspecically during purication in C. vulgaris [Honjoh, 2003].
modelLna = changeRxnBounds(modelLna,{'FBAh'},0,'b'); % light inactivates FBAh [Willard and Gibbs, 1968]
modelLna = changeRxnBounds(modelLna,{'H2Oth'},0,'u'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark (Packer 1970)
modelLna = changeRxnBounds(modelLna,{'Biomass_Duns_yellow-','Biomass_Duns_orange-'},0,'b');
modelLna = changeObjective(modelLna,'Biomass_Duns_green-');
%Solution
solutionLna = optimizeCbModel(modelLna);

%"Photoautotrophy"
A=solutionLna.f;
for i=1:length(react_flux)
    react_flux{i,2} = solutionLna.x(i);
end

%FVA
[minFlux1, maxFlux1] = fluxVariability(modelLna, 90, 'max', rxn_List);
for i = 1:length(rxn_List)
    fva_min{i,2}=minFlux1(i);
    fva_max{i,2}=maxFlux1(i);
    ind = strmatch(rxn_List{i},model.rxns,'exact');
    fba_op{i,2} = solutionLna.x(ind);
end

%% Photo Yellow
modelLna = model;
% The single PRISM reaction being used has to be commented-out below.
modelLna = changeRxnBounds(modelLna,{...
   'PRISM_solar_litho',...
    'PRISM_design_growth',...
   'PRISM_solar_exo',...
    'PRISM_fluorescent_warm_18W'...
...    'PRISM_incandescent_60W',...
    'PRISM_fluorescent_cool_215W',...
    'PRISM_metal_halide',...
    'PRISM_high_pressure_sodium',...
    'PRISM_growth_room',...
    'PRISM_white_LED',...
    'PRISM_red_LED_array_653nm',...
    'PRISM_red_LED_674nm',...
},0,'b');
modelLna = changeRxnBounds(modelLna,{'EX_o2_e'},-1000,'l');
modelLna = changeRxnBounds(modelLna,{'EX_hco3_e'},-13.54,'l');
modelLna = changeRxnBounds(modelLna,{'EX_ac_e'},0,'l');
modelLna = changeRxnBounds(modelLna,{'EX_no3_e'},-10,'l');
modelLna = changeRxnBounds(modelLna,{'EX_nh4_e'},0,'l');
modelLna = changeRxnBounds(modelLna,{'EX_starch_h'},0,'b');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGRA',starchDegAerLight,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGR2A',0,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGRB',starchDegAerLight,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGR2B',0,'u');
modelLna = changeRxnBounds(modelLna,{'PCHLDR'},0,'b'); % Crystal Structure of the Nitrogenase-like Dark Operative Protochlorophyllide Oxidoreductase Catalytic Complex [Brocker 2010]
modelLna = changeRxnBounds(modelLna,{'G6PADHh','G6PBDHh'},0,'b'); % Purication of chloroplast G6PDH had not been successful, partly because this isoform is inactivated by light because of redox modulation by thioredoxin, and partly because it easily aggregates unspecically during purication in C. vulgaris [Honjoh, 2003].
modelLna = changeRxnBounds(modelLna,{'FBAh'},0,'b'); % light inactivates FBAh [Willard and Gibbs, 1968]
modelLna = changeRxnBounds(modelLna,{'H2Oth'},0,'u'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark (Packer 1970)
modelLna = changeRxnBounds(modelLna,{'Biomass_Duns_green-','Biomass_Duns_orange-'},0,'b');
modelLna = changeObjective(modelLna,'Biomass_Duns_yellow-');
%Solution
solutionLna = optimizeCbModel(modelLna);
B=solutionLna.f;
for i=1:length(react_flux)
    react_flux{i,3} = solutionLna.x(i);
end


%FVA
[minFlux2, maxFlux2] = fluxVariability(modelLna, 90, 'max', rxn_List);

for i = 1:length(rxn_List)
    fva_min{i,3}=minFlux2(i);
    fva_max{i,3}=maxFlux2(i);
    ind = strmatch(rxn_List{i},model.rxns,'exact');
    fba_op{i,3} = solutionLna.x(ind);
end

%% Photoautotrophy Orange
modelLna = model;
% The single PRISM reaction being used has to be commented-out below.
modelLna = changeRxnBounds(modelLna,{...
   'PRISM_solar_litho',...
    'PRISM_design_growth',...
   'PRISM_solar_exo',...
    'PRISM_fluorescent_warm_18W'...
...    'PRISM_incandescent_60W',...
    'PRISM_fluorescent_cool_215W',...
    'PRISM_metal_halide',...
    'PRISM_high_pressure_sodium',...
    'PRISM_growth_room',...
    'PRISM_white_LED',...
    'PRISM_red_LED_array_653nm',...
    'PRISM_red_LED_674nm',...
},0,'b');
modelLna = changeRxnBounds(modelLna,{'EX_o2_e'},-1000,'l');
modelLna = changeRxnBounds(modelLna,{'EX_hco3_e'},-13.54,'l');
modelLna = changeRxnBounds(modelLna,{'EX_ac_e'},0,'l');
modelLna = changeRxnBounds(modelLna,{'EX_no3_e'},-10,'l');
modelLna = changeRxnBounds(modelLna,{'EX_nh4_e'},0,'l');
modelLna = changeRxnBounds(modelLna,{'EX_starch_h'},0,'b');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGRA',starchDegAerLight,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGR2A',0,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGRB',starchDegAerLight,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGR2B',0,'u');modelLna = changeRxnBounds(modelLna,{'PCHLDR'},0,'b'); % Crystal Structure of the Nitrogenase-like Dark Operative Protochlorophyllide Oxidoreductase Catalytic Complex [Brocker 2010]
modelLna = changeRxnBounds(modelLna,{'G6PADHh','G6PBDHh'},0,'b'); % Purication of chloroplast G6PDH had not been successful, partly because this isoform is inactivated by light because of redox modulation by thioredoxin, and partly because it easily aggregates unspecically during purication in C. vulgaris [Honjoh, 2003].
modelLna = changeRxnBounds(modelLna,{'FBAh'},0,'b'); % light inactivates FBAh [Willard and Gibbs, 1968]
modelLna = changeRxnBounds(modelLna,{'H2Oth'},0,'u'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark (Packer 1970)
modelLna = changeRxnBounds(modelLna,{'Biomass_Duns_green-','Biomass_Duns_yellow-'},0,'b');
modelLna = changeObjective(modelLna,'Biomass_Duns_orange-');
%Solution
solutionLna = optimizeCbModel(modelLna);
C=solutionLna.f;
for i=1:length(react_flux)
    react_flux{i,4} = solutionLna.x(i);
end

%FVA
[minFlux3, maxFlux3] = fluxVariability(modelLna, 90, 'max', rxn_List);

for i = 1:length(rxn_List)
    fva_min{i,4}=minFlux3(i);
    fva_max{i,4}=maxFlux3(i);
    ind = strmatch(rxn_List{i},model.rxns,'exact');
    fba_op{i,4} = solutionLna.x(ind);
end

%% Mixotrophy green
modelLMix = model;
% The single PRISM reaction being used has to be commented-out below.
modelLMix = changeRxnBounds(modelLMix,{...
    'PRISM_solar_litho',...
    'PRISM_solar_exo',...
...    'PRISM_incandescent_60W',...
    'PRISM_fluorescent_cool_215W',...
    'PRISM_metal_halide',...
    'PRISM_high_pressure_sodium',...
    'PRISM_growth_room',...
    'PRISM_white_LED',...
    'PRISM_red_LED_array_653nm',...
    'PRISM_red_LED_674nm',...
    'PRISM_fluorescent_warm_18W',...
    'PRISM_design_growth',...    
},0,'b');
modelLMix = changeRxnBounds(modelLMix,{'EX_o2_e',},-10,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_mal-L_e'},-0.35,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_ac_e'},0,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_starch_h'},0,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_no3_e'},-10,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_nh4_e'},0,'b');
modelLMix = changeRxnBounds(modelLMix,{'EX_hco3_e'},-11.4,'l');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGRA',starchDegAerLight,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGR2A',0,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGRB',starchDegAerLight,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGR2B',0,'u');
modelLMix = changeRxnBounds(modelLMix,{'PCHLDR'},0,'b'); % Crystal Structure of the Nitrogenase-like Dark Operative Protochlorophyllide Oxidoreductase Catalytic Complex [Brocker 2010]
modelLMix = changeRxnBounds(modelLMix,{'G6PADHh','G6PBDHh'},0,'b'); % Purication of chloroplast G6PDH had not been successful, partly because this isoform is inactivated by light because of redox modulation by thioredoxin, and partly because it easily aggregates unspecically during purication in C. vulgaris [Honjoh, 2003].
modelLMix = changeRxnBounds(modelLMix,{'FBAh'},0,'b'); % light inactivates FBAh [Willard and Gibbs, 1968]
modelLMix = changeRxnBounds(modelLMix,{'H2Oth'},0,'u'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark [Packer 1970]
modelLMix = changeRxnBounds(modelLMix,{'Biomass_Duns_yellow-','Biomass_Duns_orange-'},0,'b');
modelLMix = changeObjective(modelLMix,'Biomass_Duns_green-');
solutionLMix = optimizeCbModel(modelLMix);
D=solutionLMix.f;
for i=1:length(react_flux)
    react_flux{i,5} = solutionLMix.x(i);
end

%FVA
[minFlux4, maxFlux4] = fluxVariability(modelLMix, 90, 'max', rxn_List);

for i = 1:length(rxn_List)
    fva_min{i,5}=minFlux4(i);
    fva_max{i,5}=maxFlux4(i);
    ind = strmatch(rxn_List{i},model.rxns,'exact');
    fba_op{i,5} = solutionLMix.x(ind);
end

%%
%Mixotrophy yellow
modelLMix = model;
% The single PRISM reaction being used has to be commented-out below.
modelLMix = changeRxnBounds(modelLMix,{...
    'PRISM_solar_litho',...
    'PRISM_solar_exo',...
...    'PRISM_incandescent_60W',...
    'PRISM_fluorescent_cool_215W',...
    'PRISM_metal_halide',...
    'PRISM_high_pressure_sodium',...
    'PRISM_growth_room',...
    'PRISM_white_LED',...
    'PRISM_red_LED_array_653nm',...
    'PRISM_red_LED_674nm',...
    'PRISM_fluorescent_warm_18W',...
    'PRISM_design_growth',...    
},0,'b');
modelLMix = changeRxnBounds(modelLMix,{'EX_o2_e',},-10,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_mal-L_e'},-0.35,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_ac_e'},0,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_starch_h'},0,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_no3_e'},-10,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_nh4_e'},0,'b');
modelLMix = changeRxnBounds(modelLMix,{'EX_hco3_e'},-11.4,'l');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGRA',starchDegAerLight,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGR2A',0,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGRB',starchDegAerLight,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGR2B',0,'u');
modelLMix = changeRxnBounds(modelLMix,{'PCHLDR'},0,'b'); % Crystal Structure of the Nitrogenase-like Dark Operative Protochlorophyllide Oxidoreductase Catalytic Complex [Brocker 2010]
modelLMix = changeRxnBounds(modelLMix,{'G6PADHh','G6PBDHh'},0,'b'); % Purication of chloroplast G6PDH had not been successful, partly because this isoform is inactivated by light because of redox modulation by thioredoxin, and partly because it easily aggregates unspecically during purication in C. vulgaris [Honjoh, 2003].
modelLMix = changeRxnBounds(modelLMix,{'FBAh'},0,'b'); % light inactivates FBAh [Willard and Gibbs, 1968]
modelLMix = changeRxnBounds(modelLMix,{'H2Oth'},0,'u'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark [Packer 1970]
modelLMix = changeRxnBounds(modelLMix,{'Biomass_Duns_green-','Biomass_Duns_orange-'},0,'b');
modelLMix = changeObjective(modelLMix,'Biomass_Duns_yellow-');
solutionLMix = optimizeCbModel(modelLMix);
E=solutionLMix.f;
for i=1:length(react_flux)
    react_flux{i,6} = solutionLMix.x(i);
end

%FVA
[minFlux5, maxFlux5] = fluxVariability(modelLMix, 90, 'max', rxn_List);

for i = 1:length(rxn_List)
    fva_min{i,6}=minFlux5(i);
    fva_max{i,6}=maxFlux5(i);
    ind = strmatch(rxn_List{i},model.rxns,'exact');
    fba_op{i,6} = solutionLMix.x(ind);
end

%% Mixotrophy orange
modelLMix = model;
% The single PRISM reaction being used has to be commented-out below.
modelLMix = changeRxnBounds(modelLMix,{...
    'PRISM_solar_litho',...
    'PRISM_solar_exo',...
...    'PRISM_incandescent_60W',...
    'PRISM_fluorescent_cool_215W',...
    'PRISM_metal_halide',...
    'PRISM_high_pressure_sodium',...
    'PRISM_growth_room',...
    'PRISM_white_LED',...
    'PRISM_red_LED_array_653nm',...
    'PRISM_red_LED_674nm',...
    'PRISM_fluorescent_warm_18W',...
    'PRISM_design_growth',...    
},0,'b');
modelLMix = changeRxnBounds(modelLMix,{'EX_o2_e',},-10,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_mal-L_e'},-0.35,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_ac_e'},0,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_starch_h'},0,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_no3_e'},-10,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_nh4_e'},0,'b');
modelLMix = changeRxnBounds(modelLMix,{'EX_hco3_e'},-11.4,'l');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGRA',starchDegAerLight,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGR2A',0,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGRB',starchDegAerLight,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGR2B',0,'u');
modelLMix = changeRxnBounds(modelLMix,{'PCHLDR'},0,'b'); % Crystal Structure of the Nitrogenase-like Dark Operative Protochlorophyllide Oxidoreductase Catalytic Complex [Brocker 2010]
modelLMix = changeRxnBounds(modelLMix,{'G6PADHh','G6PBDHh'},0,'b'); % Purication of chloroplast G6PDH had not been successful, partly because this isoform is inactivated by light because of redox modulation by thioredoxin, and partly because it easily aggregates unspecically during purication in C. vulgaris [Honjoh, 2003].
modelLMix = changeRxnBounds(modelLMix,{'FBAh'},0,'b'); % light inactivates FBAh [Willard and Gibbs, 1968]
modelLMix = changeRxnBounds(modelLMix,{'H2Oth'},0,'u'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark [Packer 1970]
modelLMix = changeRxnBounds(modelLMix,{'Biomass_Duns_green-','Biomass_Duns_yellow-'},0,'b');
modelLMix = changeObjective(modelLMix,'Biomass_Duns_orange-');
solutionLMix = optimizeCbModel(modelLMix);
F=solutionLMix.f;
for i=1:length(react_flux)
    react_flux{i,7} = solutionLMix.x(i);
end

%FVA
[minFlux6, maxFlux6] = fluxVariability(modelLMix, 90, 'max', rxn_List);

for i = 1:length(rxn_List)
    fva_min{i,7}=minFlux6(i);
    fva_max{i,7}=maxFlux6(i);
    ind = strmatch(rxn_List{i},model.rxns,'exact');
    fba_op{i,7} = solutionLMix.x(ind);
end

%%%
%% Dark green
%"Heterotrophy"
%%% dark, aerobic, w/ malic acid, biomass objective
modelDa = model;
modelDa = changeRxnBounds(modelDa,'EX_photonVis_e',0,'l');
modelDa = changeRxnBounds(modelDa,{'PRISM_solar_litho','PRISM_solar_exo','PRISM_incandescent_60W','PRISM_fluorescent_cool_215W','PRISM_metal_halide','PRISM_high_pressure_sodium','PRISM_growth_room','PRISM_white_LED','PRISM_red_LED_array_653nm','PRISM_red_LED_674nm','PRISM_fluorescent_warm_18W','PRISM_design_growth'},0,'b');
modelDa = changeRxnBounds(modelDa,{'EX_o2_e','EX_co2_e'},-5,'l');
modelDa = changeRxnBounds(modelDa,'EX_hco3_e',0,'l');
modelDa = changeRxnBounds(modelDa,'EX_ac_e',0,'l');
modelDa = changeRxnBounds(modelDa,'EX_no3_e',-5,'l');
modelDa = changeRxnBounds(modelDa,'EX_nh4_e',0,'l');
modelDa = changeRxnBounds(modelDa,'EX_h_e',-5,'l');
modelDa = changeRxnBounds(modelDa,'EX_for_e',0,'b');
modelDa = changeRxnBounds(modelDa,'EX_succ_e',0,'b');
modelDa = changeRxnBounds(modelDa,{'EX_starch_h'},0,'b');
modelDa = changeRxnBounds(modelDa,'EX_mal-L_e',-0.3025,'l');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGRA',0,'u');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGR2A',starchDegAerDark,'u');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGRB',0,'u');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGR2B',starchDegAerDark,'u');
modelDa = changeRxnBounds(modelDa,{'ATPSh'},0,'b'); % Inactive in dark It is suggested that there exists in dark-adapted algae a permanent proton gradient which stimulates the charge recombination process. This proton gradient results from the hydrolysis of a pool of ATP by membranebound ATPases and collapses after the addition of TNBT. The long lifetime of this proton gradient (several hours) indicates that the ATP probably comes from the mitochondria [Joliot and Joliot, 1980].
 modelDa = changeRxnBounds(modelDa,{'GAPDH_nadp_hi'},0,'b'); % Only active in light [Buchanan 1980]
%modelDa = changeRxnBounds(modelDa,{'MDH(nadp)hi','MDHC(nadp)hr'},0,'b'); % inactive in dark [Buchanan 1980]
modelDa = changeRxnBounds(modelDa,{'PPDKh'},0,'b'); % inactive in dark [Buchanan 1980]
modelDa = changeRxnBounds(modelDa,{'PRUK'},0,'b'); % inactive in dark [Buchanan 1980; Villarejo 1995]
modelDa = changeRxnBounds(modelDa,{'RBPCh','RBCh'},0,'b'); % That Rubisco was first identified as a light-activated enzyme in early labeling studies with air-grown Chlorella cells provides evidence in support of this conclusion [Pedersen et al., 1966].
modelDa = changeRxnBounds(modelDa,{'SBP'},0,'b'); % inactive in dark [Buchanan 1980, Koziol 2013] 
modelDa = changeRxnBounds(modelDa,{'H2Oth'},0,'l'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark [Packer 1970]
modelDa = changeRxnBounds(modelDa,{'Biomass_Duns_orange-','Biomass_Duns_yellow-'},0,'b');
modelDa = changeObjective(modelDa,'Biomass_Duns_green-');
solutionDa = optimizeCbModel(modelDa);
G=solutionDa.f;
for i=1:length(react_flux)
    react_flux{i,8} = solutionDa.x(i);
end

%FVA
[minFlux7, maxFlux7] = fluxVariability(modelDa, 90, 'max', rxn_List);

for i = 1:length(rxn_List)
    fva_min{i,8}=minFlux7(i);
    fva_max{i,8}=maxFlux7(i);
    ind = strmatch(rxn_List{i},model.rxns,'exact');
    fba_op{i,8} = solutionDa.x(ind);
end

%% Dark yellow
%%% dark, aerobic, w/ malic acid, biomass objective
modelDa = model;
modelDa = changeRxnBounds(modelDa,'EX_photonVis_e',0,'l');
modelDa = changeRxnBounds(modelDa,{'PRISM_solar_litho','PRISM_solar_exo','PRISM_incandescent_60W','PRISM_fluorescent_cool_215W','PRISM_metal_halide','PRISM_high_pressure_sodium','PRISM_growth_room','PRISM_white_LED','PRISM_red_LED_array_653nm','PRISM_red_LED_674nm','PRISM_fluorescent_warm_18W','PRISM_design_growth'},0,'b');
modelDa = changeRxnBounds(modelDa,{'EX_o2_e','EX_co2_e'},-5,'l');
modelDa = changeRxnBounds(modelDa,'EX_hco3_e',0,'l');
modelDa = changeRxnBounds(modelDa,'EX_ac_e',0,'l');
modelDa = changeRxnBounds(modelDa,'EX_no3_e',-5,'l');
modelDa = changeRxnBounds(modelDa,'EX_nh4_e',0,'l');
modelDa = changeRxnBounds(modelDa,'EX_h_e',-5,'l');
modelDa = changeRxnBounds(modelDa,'EX_for_e',0,'b');
modelDa = changeRxnBounds(modelDa,'EX_succ_e',0,'b');
modelDa = changeRxnBounds(modelDa,{'EX_starch_h'},0,'b');
modelDa = changeRxnBounds(modelDa,'EX_mal-L_e',-0.3025,'l');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGRA',0,'u');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGR2A',starchDegAerDark,'u');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGRB',0,'u');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGR2B',starchDegAerDark,'u');
modelDa = changeRxnBounds(modelDa,{'ATPSh'},0,'b'); % Inactive in dark It is suggested that there exists in dark-adapted algae a permanent proton gradient which stimulates the charge recombination process. This proton gradient results from the hydrolysis of a pool of ATP by membranebound ATPases and collapses after the addition of TNBT. The long lifetime of this proton gradient (several hours) indicates that the ATP probably comes from the mitochondria [Joliot and Joliot, 1980].
 modelDa = changeRxnBounds(modelDa,{'GAPDH_nadp_hi'},0,'b'); % Only active in light [Buchanan 1980]
%modelDa = changeRxnBounds(modelDa,{'MDH(nadp)hi','MDHC(nadp)hr'},0,'b'); % inactive in dark [Buchanan 1980]
modelDa = changeRxnBounds(modelDa,{'PPDKh'},0,'b'); % inactive in dark [Buchanan 1980]
modelDa = changeRxnBounds(modelDa,{'PRUK'},0,'b'); % inactive in dark [Buchanan 1980; Villarejo 1995]
modelDa = changeRxnBounds(modelDa,{'RBPCh','RBCh'},0,'b'); % That Rubisco was first identified as a light-activated enzyme in early labeling studies with air-grown Chlorella cells provides evidence in support of this conclusion [Pedersen et al., 1966].
modelDa = changeRxnBounds(modelDa,{'SBP'},0,'b'); % inactive in dark [Buchanan 1980, Koziol 2013] 
modelDa = changeRxnBounds(modelDa,{'H2Oth'},0,'l'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark [Packer 1970]
modelDa = changeRxnBounds(modelDa,{'Biomass_Duns_green-','Biomass_Duns_orange-'},0,'b');
modelDa = changeObjective(modelDa,'Biomass_Duns_yellow-');
solutionDa = optimizeCbModel(modelDa);
H=solutionDa.f;
for i=1:length(react_flux)
    react_flux{i,9} = solutionDa.x(i);
end
%%

%FVA
[minFlux8, maxFlux8] = fluxVariability(modelDa, 90, 'max', rxn_List);

for i = 1:length(rxn_List)
    fva_min{i,9}=minFlux8(i);
    fva_max{i,9}=maxFlux8(i);
    ind = strmatch(rxn_List{i},model.rxns,'exact');
    fba_op{i,9} = solutionDa.x(ind);
end


%%Dark Orange
%%% dark, aerobic, w/ malic acid, biomass objective
modelDa = model;
modelDa = changeRxnBounds(modelDa,'EX_photonVis_e',0,'l');
modelDa = changeRxnBounds(modelDa,{'PRISM_solar_litho','PRISM_solar_exo','PRISM_incandescent_60W','PRISM_fluorescent_cool_215W','PRISM_metal_halide','PRISM_high_pressure_sodium','PRISM_growth_room','PRISM_white_LED','PRISM_red_LED_array_653nm','PRISM_red_LED_674nm','PRISM_fluorescent_warm_18W','PRISM_design_growth'},0,'b');
modelDa = changeRxnBounds(modelDa,{'EX_o2_e','EX_co2_e'},-5,'l');
modelDa = changeRxnBounds(modelDa,'EX_hco3_e',0,'l');
modelDa = changeRxnBounds(modelDa,'EX_ac_e',0,'l');
modelDa = changeRxnBounds(modelDa,'EX_no3_e',-5,'l');
modelDa = changeRxnBounds(modelDa,'EX_nh4_e',0,'l');
modelDa = changeRxnBounds(modelDa,'EX_h_e',-5,'l');
modelDa = changeRxnBounds(modelDa,'EX_for_e',0,'b');
modelDa = changeRxnBounds(modelDa,'EX_succ_e',0,'b');
modelDa = changeRxnBounds(modelDa,{'EX_starch_h'},0,'b');
modelDa = changeRxnBounds(modelDa,'EX_mal-L_e',-0.3025,'l');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGRA',0,'u');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGR2A',starchDegAerDark,'u');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGRB',0,'u');
% modelDa = changeRxnBounds(modelDa,'STARCH300DEGR2B',starchDegAerDark,'u');
modelDa = changeRxnBounds(modelDa,{'ATPSh'},0,'b'); % Inactive in dark It is suggested that there exists in dark-adapted algae a permanent proton gradient which stimulates the charge recombination process. This proton gradient results from the hydrolysis of a pool of ATP by membranebound ATPases and collapses after the addition of TNBT. The long lifetime of this proton gradient (several hours) indicates that the ATP probably comes from the mitochondria [Joliot and Joliot, 1980].
 modelDa = changeRxnBounds(modelDa,{'GAPDH_nadp_hi'},0,'b'); % Only active in light [Buchanan 1980]
%modelDa = changeRxnBounds(modelDa,{'MDH(nadp)hi','MDHC(nadp)hr'},0,'b'); % inactive in dark [Buchanan 1980]
modelDa = changeRxnBounds(modelDa,{'PPDKh'},0,'b'); % inactive in dark [Buchanan 1980]
modelDa = changeRxnBounds(modelDa,{'PRUK'},0,'b'); % inactive in dark [Buchanan 1980; Villarejo 1995]
modelDa = changeRxnBounds(modelDa,{'RBPCh','RBCh'},0,'b'); % That Rubisco was first identified as a light-activated enzyme in early labeling studies with air-grown Chlorella cells provides evidence in support of this conclusion [Pedersen et al., 1966].
modelDa = changeRxnBounds(modelDa,{'SBP'},0,'b'); % inactive in dark [Buchanan 1980, Koziol 2013] 
modelDa = changeRxnBounds(modelDa,{'H2Oth'},0,'l'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark [Packer 1970]
modelDa = changeRxnBounds(modelDa,{'Biomass_Duns_green-','Biomass_Duns_yellow-'},0,'b');
modelDa = changeObjective(modelDa,'Biomass_Duns_orange-');
solutionDa = optimizeCbModel(modelDa);
I=solutionDa.f;
for i=1:length(react_flux)
    react_flux{i,10} = solutionDa.x(i);
end

%FVA
[minFlux9, maxFlux9] = fluxVariability(modelDa, 90, 'max', rxn_List);

for i = 1:length(rxn_List)
    fva_min{i,10}=minFlux9(i);
    fva_max{i,10}=maxFlux9(i);
    ind = strmatch(rxn_List{i},model.rxns,'exact');
    fba_op{i,10} = solutionDa.x(ind);
end

%%
solution= [A,B,C;D,E,F;G,H,I];
fluxes = react_flux;
min=fva_min;
max=fva_max;
optimal=fba_op;

end






