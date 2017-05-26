classdef WingResults
%   This file is part of PROP-WING TOOLS, tools for modelling the 
%   aerodynamics of propeller wing interaction.
%   This is distributed under the GNU General Public License v3 or higher.
%   Copyright (C) 2016  Kitso Epema
%
%   WingResults
%
%   Create a WingResults with
%
%   result = WingResults() for an empty WingResults
%
%   result = WingResults('Import',wing) to import wing results from the
%               folder:
%               PROP-WING TOOLS/Wing/Examples
%               The filename should look like this:
%               Get[Wing name without spaces]Performance.m
%
%   result = WingResults('Import',wing,state,Tagname,Tagvalue)
%               Same as above, but now with a specification of the result 
%               that needs to be imported. In the performance file the 
%               Tagname and Tagvalue are used to get the correct operating 
%               conditions to import the results for. It should be 
%               structured as a function. See the examples for
%               the exact details.
%
%   result = WingResults(type,wing,state) for a computed wing result. Types
%               of methods that are available, type =
%               'Lifting Line'
%               'Lifting Line Half Wing'
%               'Tornado'
%               'AVL'
%
%   Additional functionality
%       - plt
%       - pltcompare
    
    properties
        wingName    % Wing name
        wing        % Wing
        stateName   % State name
        state       % State of the wing
        type        % Type of analysis
        alpha       % Angle of attack
        CL          % Total lift coefficient
        CLll
        CLanderson
        CLrakshith
        CLrakshith2
        CD          % Total drag coefficient
        CDll
        CDanderson
        CDrakshith
        CDrakshith2
        CDind       % Induced drag coefficient
        CDvisc
        cl          % local lift coefficient
        cd          % local total drag coefficient
        cdind       % local induced drag coefficient
        cdvisc      % local viscous drag coefficient
        gamma       % local circulation
        circVs      % Normalized circulation ( Gamma / Vinf s)
        y           % spanwise coordinate [m]
        theta       % alternative span coordinate [rad]
        yn          % local normalized spanwise coordinate
        c           % local chord
        ai          % induced angle of attack
        w           % Downwash
        aeff        % Effective angle of attack
        A           % Coefficients An
        B           % matrix B
        rhs         % Right hand side
        source      % source (in case of Import)
        note        % note (used in case of Import) 
        
        % TODO: remove variables that are not used.
        % TODO: rename to WingResult (no s, as in not plural)
    end
    
    methods
        function this = WingResults(type,wing,state,varargin)
            %TODO, put state in the varargin such that with import no state
            %has to be given.
            if nargin < 3
                state = WingState();
            end
            if nargin < 2
                wing = Wing('ellipse');
            end
            if nargin < 1                
                type = 'empty';
            end
            
            if strcmp(type,'LL')
                type = 'Lifting Line';
            end
            if strcmp(type,'avl')||strcmp(type,'Avl')
                type = 'AVL';
            end
            if strcmp(type,'LL halfwing')
                type = 'Lifting Line Half Wing';
            end
            if strcmp(type,'tornado')
                type = 'Tornado';
            end
            if strcmp(type,'import')
                type = 'Import';
            end
            
            this.type = type;
            if ~strcmp(this.type,'empty')
                this.wingName = wing.Name;
                this.wing = wing;
                this.stateName = state.Name;
                this.state = state;
                this.alpha = state.alpha;
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
                        result = fh(Tagname,Tagvalue);
                    else
                        result = fh();
                    end
                    
                    this = this.saveResults(result);
                case 'Lifting Line Sanker'
                    result = LiftingLineSanker(wing,state);
                    this = this.saveResults(result);
                case 'Lifting Line'
                    result = LiftingLine(wing,state);
                    this = this.saveResults(result);
                case 'Lifting Line Half Wing'
                    result = LiftingLineHalfWing(wing,state);
                    this = this.saveResults(result);
                case 'AVL'
                    result = runAVL(wing,state);
                    this = this.saveResults(result);
                case 'Tornado'
                    result = runTornado(wing,slipstream,wingState); 
                    this = this.saveResults(result);
                otherwise
                    display(type)
                    error('Type of import not found. Try "Import" or "Lifting Line"')
            end
        end
        function this = saveResults(this,result)
            resultNames = fieldnames(result);
            
            for n = 1:length(resultNames)
                try
                    field = char(resultNames(n));
                    this.(field) = result.(field);
                catch
                    display(field)
                    warning('Field in the results not used in WingResults')
                end
            end
        end
        function [varargout] = plt(this,Xname,Yname,LineSpec,MKsize,color)
            if nargin < 5
                MKsize = 4;
            end
            if nargin < 4
                LineSpec = 'o-';
            end 
            try
                X = this.(Xname);
            catch
                display(Xname)
                error('Could not find this.(Xname)')
            end
            try 
                Y = this.(Yname);                        
            catch
                display(Yname)
                error('Could not find this.(Yname)')
            end
%             displayName = [this.getLabelName(Yname),', ','Wing: ',this.wingName];
%             displayName = [this.getLabelName(Yname),', ','Analysis type: ',this.type];
            displayName = [this.type];
            
            h = plot(X,Y,LineSpec,'DisplayName',displayName,...
                'Linewidth',2,...
                'Markersize',MKsize);
%             if exist('color','var')
%                 h.Color = color;
%             end
            xlabel(this.getLabelName(Xname))
            ylabel(this.getLabelName(Yname))
            grid on
            varargout{1} = h;
        end
        function [] = pltcompare(this,Xname,Yname,varargin)
            Linespec = {'--o',...
                '-.^',...
                ':s',...
                '--d',...
                '-.*',...
                ':x'};
            if nargin > 3
                if length(varargin) > length(Linespec)
                    warning('pltcompare cannot handle more than 3 additional datasets')
                else
                    for n = 1:length(varargin)
                        result(n) = varargin{n};
                    end
                end
            end
            h = this.plt(Xname,Yname);
            h.LineWidth = 3;
            ylabel(this.getLabelName(Yname))
            color = get(h,'Color');

            for n = 1:length(varargin)
                if n < length(Linespec)
                    h = result(n).plt(Xname,Yname,char(Linespec(n)),10,color);
                    h.LineWidth = 1.5;
                else
                    display(['Dataset number ',n,' is not plotted.'])
                end
            end
        end
        function [] = pltfinish(~)            
            legend(gca,'show')
            set(gca,'fontsize',18,...
               'Linewidth',1.5)
            grid on 
        end 
        function labelName = getLabelName(~,Xname)
            
            switch Xname
                case 'y'
                    labelName = 'Normalized span coordinate, y/b/2 [-]';
                case 'cl'
                    labelName = 'C_l [-]';
                case 'gamma'
                    labelName = 'Circulation';
                case 'c'
                    labelName = 'Chord [m]';
                case 'circVs'
                    labelName = 'Normalized circulation \Gamma / V_{\infty} 2 / b [-]';
                otherwise
                    labelName = Xname;
            end
                    
        end
        
    end
end

