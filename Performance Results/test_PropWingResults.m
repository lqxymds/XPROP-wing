clc
clearvars
close all

%% Test PropWingResults class

% Wing
wing = Wing('PROWIM',50);

% Propeller
beta = 25;
prop = Propeller('PROWIM',beta);

Vinf = 20;
alpha = 3;
stateName = strcat('V = ',num2str(Vinf),', alpha = ',num2str(alpha));
wingState = WingState(wing,stateName,Vinf,alpha);

% Slipstream
orientation = 'inboard up';
slipstream = Slipstream('PROWIM','J','0.85 VLM4D');
slipstream = slipstream.setRadius(prop.R);
slipstream.y0 = prop.y0;
slipstream = slipstream.setOrientation(orientation);

results = PropWingResults('LL in slipstream',wing,wingState,slipstream);
woProp = WingResults('LL',wing,wingState);

%% Plot
figure
hold all
results.plt('wing','y','gamma');
woProp.plt('y','gamma');
results.pltfinish

