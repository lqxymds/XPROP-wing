clc
clearvars
% close all

%% Test PROWIM

alpha = 4;
Vinf = 80;
Tc = 0.168;
Tc2 = 0.168;%0.42;

% Wing
wing = Wing('PROWIM',50);

% Propeller
beta = 25;
prop = Propeller('PROWIM',beta);
% prop.reversed = 0;

wingStateName = ['V = ',num2str(Vinf),', alpha = ',num2str(alpha)];
wingState = WingState(wing,wingStateName,Vinf,alpha);

% Slipstream 1
orientation = 'inboard up';
% slipstream = Slipstream('propeller',prop,propState);
% slipstream = SetRakshithSlipstream();
slipstream = Slipstream('PROWIM','J','0.8');
% slipstream = Slipstream('PROWIM','J','0.85 isolated');
% slipstream = Slipstream('PROWIM','J','0.92');
% slipstream.vtVinf = -slipstream.vtVinf;
slipstream.y0 = prop.y0;
slipstream = slipstream.setRadius(prop.R);
% slipstream.vaVinf = 0 .* slipstream.vaVinf;
% slipstream.vtVinf = 0 .* slipstream.vtVinf;
% slipstream.vtVinf = slipstream.vtVinf./2;
slipstream = slipstream.setOrientation(orientation);

% Slipstream 2
slipstream2 = Slipstream('PROWIM','J','0.92');
slipstream2.y0 = prop.y0;
slipstream2 = slipstream2.setRadius(prop.R);
% slipstream2.vtVinf = slipstream2.vtVinf./2;
slipstream2 = slipstream2.setOrientation(orientation);

% Numerical results
results = PropWingResults('LL in slipstream',wing,wingState,slipstream);
results.prop = prop;
results2 = results;
results2.wngres.cl = results.wngres.cl;
results2.type = 'LL in slipstream w/o Vp in C_l';
% results.wngres.cl = results.wngres.cl .* Vinf ./ results.wngres.V;

% results2 = PropWingResults('LL in slipstream 2',wing,wingState,slipstream);
% results3 = PropWingResults('LL Rakshith',wing,wingState,slipstream);
% results4 = PropWingResults('LL Rakshith 2',wing,wingState,slipstream);


% Experimental results

exp = PropWingResults('Import',wing,wingState,num2str(Tc),num2str(alpha));
exp.prop = prop;
exp.type = 'Experiment';
yn = exp.wngres.yn;
V = interp1(results.wngres.yn,results.wngres.V,yn,'linear','extrap');
exp.wngres.gamma = 0.5 * exp.wngres.cl * wing.c(1) * Vinf.^2 ./ V;
exp2 = PropWingResults('Import',wing,wingState,num2str(Tc2),num2str(alpha));
exp2.prop = prop;
exp2.type = 'Experiment';
exp2.wngres.gamma = 0.5 * exp2.wngres.cl * wing.c(1) * Vinf.^2 ./V;


woProp = PropWingResults('LL in slipstream',wing,wingState,Slipstream('empty'));
woProp.type = 'LL in slipstream, w/o prop';

% prop results
% prpres = PropResults(prop,'BEM',propState);

%% Plot

figure
hold all
exp2.pltcompare('wing','yn','gamma',...
    results,...
    results2,...
    ...exp,...
    ...results2,...
    ...results3,...
    ...results4,...
    woProp)
exp2.pltfinish()

figure
hold all
plot(slipstream.r,slipstream.vaVinf)
%plot(slipstream2.r,slipstream2.vaVinf)
legend('J = 0.85','J = 0.92')
grid on

% figure
% hold all
% plot(slipstream.r,slipstream.vtVinf)
% plot(slipstream2.r,slipstream2.vtVinf)
% legend('J = 0.85','J = 0.92')
% grid on

% figure()
% hold all
% results.pltcompare('wing','y','ai')%,...
%     ...results2,...
%     results3)%,...
%     results4)
% results.pltfinish()

%%

r = slipstream2.r;
va = slipstream2.vaVinf;
vt = slipstream2.vtVinf;

dr = diff(r);
dva = diff(va);

figure
hold all
plot(dr)
plot(dva)
grid on




