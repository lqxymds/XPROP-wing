clc
clearvars
close all

%% Compare Tornado with Experiment

% Make sure to select the correct figure

alpha = 4;
Vinf = 50;
Tc = 0.168;
J = 0.8;

% Propeller
beta = 25;
prop = Propeller('PROWIM',beta);
% prop.y0 = 0.325; % Prowim pylon

% Wing
wing = Wing('PROWIM',50);
% wing = Wing('PROWIMPylon',50);
c = wing.c(1);

% Wing operating state
wingStateName = ['V = ',num2str(Vinf),', alpha = ',num2str(alpha)];
wingState = WingState(wing,wingStateName,Vinf,alpha);

slipstream = Slipstream('PROWIM','J','0.8');
slipstream.y0 = prop.y0;
slipstream = slipstream.setRadius(prop.R);
slipstream.vtVinf = -slipstream.vtVinf;

%% numerical result
num = PropWingResults('LL in slipstream',wing,wingState,slipstream);
% num.wngres.cl = num.wngres.cl2; % use F = rho * Vinf * Gamma (instead of F = rho*V(y)*Gamma)
num.wngres.l = num.wngres.cl .* 0.5 * wingState.rho * Vinf.^2 .*c;
tor = PropWingResults('Tornado',wing,wingState,slipstream);
tor.wngres.l = tor.wngres.cl .* 0.5 * wingState.rho * Vinf.^2 .*c;

%% experimental result
exp = PropWingResults('Import',wing,wingState,num2str(Tc),num2str(alpha));
% exp = PropWingResults('Import',wing,wingState,num2str(J),num2str(alpha)); % Prowim Pylon
exp.prop = prop;
yn = exp.wngres.yn;
Cl = exp.wngres.cl;
V = interp1(num.wngres.yn,num.wngres.V,yn,'linear','extrap');
exp.wngres.gamma = exp.wngres.cl .* Vinf.^2 .* c ./ ( 2 * V);
exp.wngres.l = exp.wngres.cl .* 0.5 * wingState.rho * Vinf.^2 .*c;
exp.type = 'Experimental result, \Gamma = C_l V_{\infty}^2 c / (2 V)';
exp2 = exp;
exp2.wngres.gamma = exp.wngres.cl .* Vinf .* c ./ (2);
exp2.wngres.l = exp.wngres.cl .* 0.5 * wingState.rho * Vinf.^2 .*c;
exp2.type = 'Experimental result, \Gamma = C_l V_{\infty} c / (2)';
y = yn .* wing.B/2;

%% Kroo optimum
CLdes = num.wngres.CL;
wingState.CLdes = CLdes;
result_kroo = PropWingResults('Kroo',wing,wingState,slipstream);
result_kroo.type = 'Kroo optimum';
result_kroo.wngres.CDvisc = 0;
result_kroo.wngres.CDtot = 0;
result_kroo.wngres.l = result_kroo.wngres.cl .* 0.5 * wingState.rho * Vinf.^2 .*c;
y_kroo = result_kroo.wngres.yn .* wing.B/2;
Cl_kroo = result_kroo.wngres.cl;

%% Plot
% hold all
% plot(y,Cl,'o')

% plot(y_kroo,Cl_kroo,'s')
% legend('Experiment','Kroo')

%%
figure
hold on
exp.pltcompare('wing','yn','cl',...
    exp2,...
    result_kroo,...
    num,...
    tor)
exp.pltfinish
exp.pltprop

