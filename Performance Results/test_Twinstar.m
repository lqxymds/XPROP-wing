clc
clearvars
close all

%% Test Twinstar

alpha = 2.34;
Vinf = 20;
CLdes = 0.2;

% Wing
wing = Wing('Twinstar',50);

% Propeller
prop = Propeller('APC9x6SF');
prop.y0 = 0.2127; %mm

wingStateName = ['V = ',num2str(Vinf),', alpha = ',num2str(alpha)];
wingState = WingState(wing,wingStateName,Vinf,alpha);
wingState.CLdes = CLdes; % Set design lift coefficient

% Slipstream
slipstream = Slipstream('APC9x6SF','J','0.7');
slipstream.vtVinf = -slipstream.vtVinf;
slipstream.vaVinf = 1.*slipstream.vaVinf;
slipstream.vtVinf = 1.*slipstream.vtVinf;
% slipstream = Slipstream('empty');
slipstream.y0 = prop.y0;
slipstream = slipstream.setRadius(prop.R);

% Numerical results
result_num = PropWingResults('LL in slipstream',wing,wingState,slipstream);
result_num.prop = prop;

% Kroo optimum
result_kroo = PropWingResults('Kroo',wing,wingState,slipstream);
result_kroo.type = 'Kroo optimum';
result_kroo.wngres.CDvisc = 0;
result_kroo.wngres.CDtot = 0;

%% 2D polar

CDvPolar = getCDviscPolar(wing,wingState);

% figure
% plot(CDvPolar.cdv,CDvPolar.CL,...
%     'Linewidth',2)
% grid on


%% Display CL CD
result_num.displaycompare('wing','CDvisc','CDind','CDtot',...
    result_kroo)

rho = wingState.rho;
S = wing.getS;
n = 5000/60;
R = prop.R;

% Drag
CDtot = result_num.wngres.CDtot;
D = 0.5 * rho * Vinf^2 * S * CDtot;

% Thrust
r = slipstream.r;
Ve = Vinf * (1 + slipstream.vaVinf);
T = 0.5 * rho * 2*pi *trapz(r,(Ve.^2 - Vinf^2).*r);
Ct = 0.02;
T2 = rho * (R*2)^4 * n^2 * Ct;

%% Plot
figure
hold all
result_num.pltcompare('wing','yn','cl',...
    result_kroo)
result_num.pltfinish