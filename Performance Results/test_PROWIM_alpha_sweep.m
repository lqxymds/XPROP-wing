clc
clearvars
close all

%% Test PROWIM alpha sweep

% Settings
alpha = 0:2:10;
Vinf = 50;
Tc = 0.168;

% Wing
wing = Wing('PROWIM',50);

% Propeller
beta = 25;
prop = Propeller('PROWIM',beta);

wingStateName = ['V = ',num2str(Vinf),', alpha sweep'];
wingState = WingState(wing,wingStateName,Vinf,alpha);

% Slipstream
orientation = 'inboard up';
% slipstream = Slipstream('PROWIM','J','0.8');
slipstream = Slipstream('PROWIM','J','0.85 VLM4D');
% slipstream = Slipstream('PROWIM','J','0.85 isolated');
% slipstream = Slipstream('PROWIM','J','0.92');
slipstream = slipstream.setRadius(prop.R);
slipstream.y0 = prop.y0;
slipstream = slipstream.setOrientation(orientation);

%% Experimental results
exp_propoff = PropWingResults('Import',wing,wingState,num2str(Tc),'Sweep prop off');
exp_propoff.type = 'Experiment - prop off';
CD0_off = exp_propoff.wngres.CD(1);
exp_propoff.wngres.CDind = exp_propoff.wngres.CD - CD0_off;
exp_propon = PropWingResults('Import',wing,wingState,num2str(Tc),'Sweep prop on');
exp_propon.type = 'Experiment - prop on';
CD0_on = exp_propon.wngres.CD(1);
exp_propon.wngres.CDind = exp_propon.wngres.CD - CD0_on;

%% Numerical results
tor_propon = PropWingResults('Tornado',wing,wingState,slipstream);
tor_propon.type = 'VLM - prop on';
tor_propoff = PropWingResults('Tornado',wing,wingState,Slipstream('empty'));
tor_propoff.type = 'VLM - prop off';


%% Lifting line results
ll_propon = PropWingResults('LL in Slipstream',wing,wingState,slipstream);
ll_propon.type = 'LL - prop on';

ll_propoff = PropWingResults('LL in Slipstream',wing,wingState,Slipstream('empty'));
ll_propoff.type = 'LL - prop off';

%% Kroo optimum prop on
wingState.CLdes = tor_propon.wngres.CL;
kroo_propon = PropWingResults('Kroo',wing,wingState,slipstream);
kroo_propon.type = 'Kroo - prop on';

kroo_propoff = PropWingResults('Kroo',wing,wingState,Slipstream('empty'));
kroo_propoff.type = 'Kroo - prop off';

% %% CL
% for n = 1:length(tor_propon.wngres.alpha)
%     Cl = tor_propon.wngres.cl(:,n);
%     yn = tor_propon.wngres.yn(:,n);
%     V = tor_propon.wngres.V(:,n);
%     CL2(n) = trapz(yn,Cl .* (V./Vinf));
%     AoA(n) = tor_propon.wngres.alpha(n);
% end
% 
% %% CD trefftz
% CDtr = tor_propon.wngres.CDtrefftz;
% CL = tor_propon.wngres.CL;

%% Plot prop on
figure
hold all
exp_propon.pltcompare('wing','CDind','CL',...
    tor_propon,...
    ll_propon,...
    kroo_propon);
exp_propoff.pltfinish

%% Plot prop off
figure
hold all
exp_propoff.pltcompare('wing','CDind','CL',...
    tor_propoff,...
    ll_propoff,...
    kroo_propoff);
exp_propoff.pltfinish

%%
figure
hold all
exp_propon.pltcompare('wing','alpha','CDind',...
    tor_propon,...
    ll_propon,...
    kroo_propon);
exp_propoff.pltfinish

%%
% figure
% hold all
% plot(num_propon.wngres.yn(:,3),num_propon.wngres.cl(:,3))