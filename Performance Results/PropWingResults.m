classdef PropWingResults
%   This file is part of PROP-WING TOOLS, tools for modelling the 
%   aerodynamics of propeller wing interaction.
%   This is distributed under the GNU General Public License v3 or higher.
%   Copyright (C) 2016  Kitso Epema
%
%   PropWingResults
%
%   Create a PropWingResult with:
%
%   result = PropWingResults() for an empty propeller wing result
%
%   result = PropWingResults('Import',wing,wingState) to import propeller 
%               wing results from the folder:
%               PROP-WING TOOLS/Prop-Wing/Examples
%               The filename should look like this:
%               Get[Name without spaces]Performance.m
%
%   result = PropWingResults('Import',wing,wingState,Tagname,Tagvalue)
%               Same as above, but now with a specification of the result 
%               that needs to be imported. In the performance file the 
%               Tagname and Tagvalue are used to get the correct operating 
%               conditions to import the results for. It should be 
%               structured as a function. See the examples for
%               the exact details.
%
%   result = PropWingResults(type,wing,wingState,slipstream) with type the
%               calculation method. The calculation method options are type=
%               - 'Lifting Line in Slipstream'
%               - 'Lifting Line in Slipstream 2'
%               - 'Tornado'
%               - 'Kroo'
%               and slipstream an slipstream object
%

    properties
        wingName        % Wing name
        wing            % Wing
        propName        % Propeller name
        prop            % Propeller
        wingStateName   % Wing state name
        wingState       % Wing operating state
        propStateName   % Propeller state name
        propState       % Propeller operating state
        wngres          % Wing results
        prpres          % Propeller results
        type
        slipstream 
    end
    
    methods
        function this = PropWingResults(type,wing,wingState,varargin)
            if nargin < 3
                wingState = WingState();
            end
            if nargin < 2
                wing = Wing('PROWIM');
            end
            if nargin < 1
                type = 'empty';
            end
            
            if strcmp(type,'LL in slipstream')
                type = 'Lifting Line in Slipstream';
            end
            if strcmp(type,'LL in Slipstream')
                type = 'Lifting Line in Slipstream';
            end
            if strcmp(type,'LL in Slipstream 2')
                type = 'Lifting Line in Slipstream 2';
            end
            if strcmp(type,'import')
                type = 'Import';
            end
            if strcmp(type,'tornado')
                type = 'Tornado';
            end
            if strcmp(type,'kroo')
                type = 'Kroo';
            end
            
            this.type = type;
            if ~strcmp(this.type,'empty')
                this.wingName = wing.Name;
                this.wing = wing;
                
                this.wingStateName = wingState.Name;
                this.wingState = wingState;
            end
            
            calc_methods = {'Lifting Line in Slipstream'
                            'Lifting Line in Slipstream 2'
                            'Tornado'
                            'VLM Alba'
                            'Kroo'};
            if max(strcmp(type,calc_methods))
               calculation_method = type;
               type = 'Calculation';
            end
                
            switch type
                case 'empty'
                case 'Import'
                    Name = this.wingName;
                    Name = Name(~isspace(Name));
                    str = strcat('Get',Name,'Performance');
                    fh = str2func(str);
                    
                    if nargin > 3
                        Tagname = varargin{1};
                        Tagvalue = varargin{2};
                        wingResult = fh(Tagname,Tagvalue);
                    else
                        wingResult = fh(varargin);
                    end
                    this = this.saveWingResults(wingResult);
                case 'Calculation'
                    
                    second_input = varargin{1};
                    if isa(second_input,'Slipstream')
                        slipstream = second_input;
                    elseif isa(second_input,'Propeller')
                        display(second_input)
                        error('This method is not yet implemented for a Propeller, try a slipstream instead.')
                    else
                        display(second_input)
                        error('Type of input not as expected')
                    end
                    
                    this.slipstream = slipstream;
                    n_alpha_entries = length(wingState.alpha); 
                    
                    if length(varargin) > 1
                        for n = 2:length(varargin)
                            settings(n-1) = varargin{n};
                        end
                    else
                        settings = [];
                    end
                    
                    if n_alpha_entries > 1
                        temp_wingState = wingState;
                        for n = 1:n_alpha_entries
                            temp_wingState.alpha = wingState.alpha(n);
                            if isempty(wingState.CLdes)
                                
                            else
                                temp_wingState.CLdes = wingState.CLdes(n);
                            end
                            temp_wingResult = this.performCalculation(calculation_method,temp_wingState,settings);
                        
                            if n == 1
                                wingResult = temp_wingResult;
                            else
                                wingResult = this.addWingResult(wingResult,temp_wingResult);
                            end
                        end
                    elseif n_alpha_entries == 1
                        wingResult = this.performCalculation(calculation_method,wingState,settings);
                    else
                        display(wingState.alpha)
                        error('Did not find expected number of angle of attack entries.')
                    end
                    
                    this = this.saveWingResults(wingResult);                   
                otherwise
                    display(type)
                    error('Method for analysis not found')
            end
        end
        function wingResult = performCalculation(this,calculation_method,temp_wingState,settings)
            
            switch calculation_method
                case 'Lifting Line in Slipstream'
                    wingResult = LLinSlipstream(this.wing,this.slipstream,temp_wingState);
                    
                case 'Lifting Line in Slipstream 2'
                    CASE = settings(1);
                    wingResult = LLinSlipstream2(this.wing,this.slipstream,temp_wingState,CASE);
                    
                case 'Tornado'
                    if isempty(settings)
                        wingResult = runTornado(this.wing,this.slipstream,temp_wingState);
                    elseif length(settings) == 1
                        latticetype = settings{1};
                        wingResult = runTornado(this.wing,this.slipstream,temp_wingState,latticetype);
                    elseif length(settings) == 2
                        latticetype = settings(1);
                        CASE = settings(2);
                        wingResult = runTornado(this.wing,this.slipstream,temp_wingState,latticetype,CASE);
                    else
                        display(settings)
                        error('Number of entries in settings not as expected.')
                    end
                case 'VLM Alba'
                    if isempty(settings)
                        wingResult = runVLMAlba(this.wing,this.slipstream,temp_wingState);
                    elseif length(settings) == 1
%                         latticetype = settings{1};
                        wingResult = runVLMAlba(this.wing,this.slipstream,temp_wingState,settings);
                    end
                case 'Kroo'
                    wingResult = KrooPropWing(this.wing,this.slipstream,temp_wingState);
                    
                otherwise
                    display(calculation_method)
                    error('Method for analysis not found')
            end
        end
        function this = saveWingResults(this,result)
            resultNames = fieldnames(result);
            
            for n = 1:length(resultNames)
                try
                    field = char(resultNames(n));
                    this.wngres.(field) = result.(field);
                catch
                    display(field)
                    warning('Field in results not used in PropWingResults')
                end
            end
        end
        function result = addWingResult(~,result,add_result)
            resultNames = fieldnames(add_result);
            result_n = length(result.alpha)+1;
            for n = 1:length(resultNames)
                    field = char(resultNames(n));
                    [~,nc] = size(add_result.(field));
                    
                    if strcmp(field,'gamma')
                        
                    end
                    if nc == 1
                        result.(field)(:,result_n) = add_result.(field);
                    elseif nc > 1
                        result.(field)(:,:,result_n) = add_result.(field);
                    end
            end
            
        end
        function this = savePropResults(this,result)
            resultNames = fieldnames(result);
            
            for n = 1:length(resultNames)
                try
                    field = char(resultNames(n));
                    this.prpres.(field) = result.(field);
                catch
                    display(field)
                    warning('Field in results not used in PropWingResults')
                end
            end
        end
        function [varargout] = plt(this,type,Xname,Yname,LineSpec,MKsize,color)
            % [varargout] = plt(this,type,Xname,Yname)
            % [varargout] = plt(this,type,Xname,Yname,LineSpec)
            % [varargout] = plt(this,type,Xname,Yname,LineSpec,MKsize)
            % [varargout] = plt(this,type,Xname,Yname,LineSpec,MKsize,color)
            %
            
            if nargin < 6
                MKsize = 4;
            end
            if nargin < 5
                LineSpec = 'o-';
            end
            switch type
                case 'wing'
                    type = 'wngres';
                case 'prop'
                    type = 'prpres';
                    display(type)
                    warning('Plotting propeller results in PropWingRestuls not yet tested')
                otherwise
                    
                    if strcmp('wingspan',type(1:end-2))
                        entry = str2double(type(end));
                        type = 'wingspan';
                        type2 = 'wngres';
                    else
                        display(type)
                        error('Results type not found.')
                    end
            end
            
            try
                if strcmp(type,'wingspan')
                    Y = this.(type2).(Yname)(:,entry);
                else
                    Y = this.(type).(Yname);
                end
            catch
                display(Y)
                warning('Could not find this.(Yname)')
            end
            if strcmp(Xname,'n')
                if strcmp(type,'wingspan')
                    X = 1:length(this.(type2).(Yname)(:,entry));
                else
                    X = 1:length(this.(type).(Yname));
                end
            else
                try
                    if strcmp(type,'wingspan')
                        X = this.(type2).(Xname)(:,entry);
                    else
                        X = this.(type).(Xname);
                    end
                catch
                    warning('Could not find this.(Xname)')
                end
            end
            
            displayName = [this.type];
            
            h = plot(X,Y,LineSpec,'DisplayName',displayName,...
                'Linewidth',2,...
                'Markersize',MKsize);
            if exist('color','var')
                h.Color = color;
            end
            xlabel(this.getLabelName(Xname))
            ylabel(this.getLabelName(Yname))
            varargout{1} = h;
        end
        function [] = pltcompare(this,type,Xname,Yname,varargin)
            Linespec = {'--o',...
                '-.^',...
                ':s',...
                '--d',...
                '-.*',...
                ':x',...
                '--+',...
                '-.d',...
                ':v'};
            if nargin > 3
                if length(varargin) > length(Linespec)
                    warning('pltcompare cannot handle more than 9 additional datasets')
                else
                    for n = 1:length(varargin)
                        if length(varargin{n}) == 1
                            result(n) = varargin{n};
                        elseif length(varargin{n}) > 1
                            result_list = varargin{n};
                            for m = 1:length(result_list)
                                
                                result(n+m-1) = result_list(m);
                            end
                        end
                    end
                end
            end
            h = this.plt(type,Xname,Yname);
            h.LineWidth = 3;
            ylabel(this.getLabelName(Yname))
            color = get(h,'Color');
            
            if ~isempty(varargin)
                for n = 1:length(result)
                    if n < length(Linespec)
                        h = result(n).plt(type,Xname,Yname,char(Linespec(n)),10);
                        h.LineWidth = 1.5;
                    else
                        display(['Dataset number ',n,' is not plotted.'])
                    end
                end
            end
        end
        function [] = pltprop(this)
            if isempty(this.prop)
            else
                pyl = (this.prop.y0 - this.prop.R)*2/this.wing.B;
                pyr = (this.prop.y0 + this.prop.R)*2/this.wing.B;
                pyc = (this.prop.y0)*2/this.wing.B;
                plot([pyl pyl],ylim,'k',...
                    [pyc pyc],ylim,'k',...
                    [pyr pyr],ylim,'k')
            end
        end
        function [] = pltfinish(this)
            legend(gca,'show')
            set(gca,'fontsize',22,...
                'Linewidth',1.5)
            grid on
            this.pltprop()
        end
        function labelName = getLabelName(~,Xname)
            switch Xname
                case 'y'
                    labelName = 'Span coordinate [m]';
                case 'yn'
                    labelName = 'Normalized span coordinate 2y/b [-]';
                case 'alpha'
                    labelName = 'Angle of attack \alpha [deg]';
                case 'l'
                    labelName = 'Local lift [N/m]';
                case 'cl'
                    labelName = 'Local lift coefficient C_l [-]';
                case 'clcmac'
                    labelName = 'Normalized local lift C_l c / c_{MAC} [-]';
                case 'cdind'
                    labelName = 'Local induced drag coefficient C_{d_{ind}} [-]';
                case 'CL'
                    labelName = 'Lift coefficient C_L [-]';
                case 'CL2'
                    labelName = 'C_L^2 [-]';
                case 'CDind'
                    labelName = 'Induced drag coefficient C_{D_{ind}} [-]';
                case 'gamma'
                    labelName = 'Circulation [m^2/s]';
                case 'c'
                    labelName = 'Chord [m]';
                case 'circVs'
                    labelName = 'Normalized circulation \Gamma / V_{\infty} 2 / b [-]';
                case 'ai'
                    labelName = 'Induced angle of attack \alpha_i [deg]';
                case 'w'
                    labelName = 'Wing induced downwash w_{wing} [m/s]';
                otherwise
                    labelName = Xname;
            end
        end
        function [] = displaycompare(this,type,Var1,Var2,Var3,varargin)
            switch type
                case 'wing'
                    fprintf('\r\n')
                    fprintf('%s\t\t %s\t\t %s\t %s\t \r\n','Var',Var1,Var2,Var3)
                    fprintf('%s\t %2.4f\t %2.4f\t %2.4f\t \n',this.type,this.wngres.(Var1),this.wngres.(Var2),this.wngres.(Var3))
                    
                    if ~isempty(varargin)
                        for n = 1:length(varargin)
                            result = varargin{n};
                            fprintf('%s\t %2.4f\t %2.4f\t %2.4f\t \n',result.type,result.wngres.(Var1),result.wngres.(Var2),result.wngres.(Var3))
                        end
                    end                    
                case 'propeller'
                    error('Display compare not yet implemented for propeller.')
                otherwise
                    display(type)
                    error('Type of display compare not found.')
            end
        end
    end
    methods(Static)
        function this = circFromWTtest(wing,prop_state,type)
            if nargin < 3
                type = 'piv';
            end
            if nargin < 2
                prop_state = 'off';
            end
            method = 'circ';
            this = PropWingResults.fromWTtest(wing,prop_state,method,type);
        end
        function this = momFromWTtest(wing,prop_state,type,p_method)
            if nargin < 4
                p_method = 'isen';
            end
            if nargin < 3
                type = 'piv';
            end
            if nargin < 2
                prop_state = 'off';
            end
            method = 'momentum';
            this = PropWingResults.fromWTtest(wing,prop_state,method,type,p_method);
        end
        function this = fromWTtest(wing,prop_state,method,type,p_method)
            if nargin < 5
                p_method = 'isen';
            end
            if nargin < 4
                type = 'piv';
            end
            if nargin < 3
                method = 'circ';
            end
            if nargin < 2
                prop_state = 'off';
            end
            
            wingName = wing.Name;
            folder = getResultFolderList(wingName,prop_state,'','-');
            folder = removeSlipstreamEdge(folder);
            
            for n = 1:length(folder)
                % span station
                y(n) = str2double(getStationFromFolder(folder{n}))./1000;
                
                % get ojf data
                ojfdata(n) = OJFdataAvg.fromFolder(folder{n});
            end
            
            % parameters from OJFdata
            Vinf = ojfdata.Uinf;
            rho = ojfdata.rho;
            q = 0.5 * rho * Vinf^2;
            
            % parameters from wing
            s = wing.halfspan;
            c = interp1(wing.b,wing.c,y); 
            
            % calculate other values
            this = PropWingResults('empty');
            
            this.wngres.alpha = VectorField.getWTaoa;
            
            this.wngres.y = y;
            this.wngres.yn = y ./ s;
            this.wngres.c = c;
            
            this.wingName = wingName;
            this.wing = wing;
            
            switch method
                case 'circ'
                    for n = 1:length(folder)
                        % Get circulation
                        [circ_raw,~] = openCirculation(folder{n},type);
                        circ_mean(n) = mean(circ_raw);
                        circ_stdev(n) = std(circ_raw);
                        
                        circ_stdev(n) = circ_stdev(n) + 0.03 .* circ_mean(n);
                    end
            
                    % parameters from circulation
                    gamma = circ_mean;

                    this.wngres.gamma = gamma;
                    this.wngres.gamma_stdev = circ_stdev;
                    this.wngres.circVs = gamma ./ (Vinf .* s);
                    this.wngres.circVs_stdev = circ_stdev ./ (Vinf .* s);
                    this.wngres.l = rho .* Vinf .* gamma;
                    this.wngres.l_stdev = rho .* Vinf .* circ_stdev;
                    this.wngres.cl = 2 .* gamma ./ (Vinf .* c);
                    this.wngres.cl_stdev = 2 .* circ_stdev ./ (Vinf .* c);
                    this.type = 'Wind tunnel test (circulation)';
                    
                    this.wngres.clcmac = this.wngres.cl .* c ./wing.MAC;
                    this.wngres.clcmac_stdev = this.wngres.cl_stdev .* c ./wing.MAC;
                    
                case 'momentum'
                    for n = 1:length(folder)
                        % Get momentum
                        [mom_l_raw,mom_d_raw,~] = openMomentum(folder{n},type,p_method);
                        
                        mom_l_mean(n) = mean(mom_l_raw,'omitnan');
                        mom_l_stdev(n) = std(mom_l_raw,'omitnan');
                        
                        mom_d_mean(n) = mean(mom_d_raw,'omitnan');
                        mom_d_stdev(n) = std(mom_d_raw,'omitnan');
                    end
                    
                    this.wngres.l = mom_l_mean;
                    this.wngres.l_stdev = mom_l_stdev;
                    this.wngres.d = mom_d_mean;
                    this.wngres.d_stdev = mom_d_stdev;
                    
                    this.wngres.cl = this.wngres.l ./ (q.*c);
                    this.wngres.cl_stdev = this.wngres.l_stdev ./ (q.*c);
                    this.wngres.gamma = this.wngres.l ./ (rho * Vinf);
                    this.wngres.gamma_stdev = this.wngres.l_stdev ./ (rho * Vinf);
                    this.wngres.circVs = this.wngres.gamma ./ (Vinf .* s);
                    this.wngres.circVs_stdev = this.wngres.gamma_stdev ./ (Vinf .* s);
                    
                    this.wngres.cd = this.wngres.d ./ (q.*c);
                    this.wngres.cd_stdev = this.wngres.d_stdev ./ (q.*c);
                    
                    this.type = 'Wind tunnel test (momentum)';
                    
                    this.wngres.clcmac = this.wngres.cl .* c ./wing.MAC;
                    this.wngres.clcmac_stdev = this.wngres.cl_stdev .* c ./wing.MAC;
                    
                otherwise
                    display(method)
                    error('Could not find this method.')
            end
        end
        
    end
end

