%/************************************************************/%
%/*    NAME: Blake Cole                                      */%
%/*    ORGN: MIT                                             */%
%/*    FILE: xfoil.m                                         */%
%/*    DATE: 19 DEC 2018                                     */%
%/*    NOTE: Adapted from Louis Edelman's implementation of  */%
%/*          Gus Brown's Xfoil Interface. Modified for MacOS */%
%/************************************************************/%

function [pol,foil] = xfoil(coord,Alpha,Re,Mach,varargin)
% Run XFoil and return the results.
% [polar,foil] = xfoil(coord,alpha,Re,Mach,{extra commands})
%
% Inputs:
%    coord: Normalised foil co-ordinates (n by 2 array, of x & y
%           from the TE-top passed the LE to the TE bottom)
%           or a filename of the XFoil co-ordinate file
%           or a NACA 4 or 5 digit descriptor (e.g. 'NACA0012')
%    alpha: Angle-of-attack, can be a vector for an alpha polar
%       Re: Reynolds number (use Re=0 for inviscid mode)
%     Mach: Mach number
% extra commands: Extra XFoil commands
%           The extra XFoil commands need to be proper xfoil commands
%           in a character array. e.g. 'oper/iter 150'
%
% The transition criterion Ncrit can be specified using the
% 'extra commands' option as follows,
% foil = xfoil('NACA0012',10,1e6,0.2,'oper/vpar n 12')
%
%   Situation           Ncrit
%   -----------------   -----
%   sailplane           12-14
%   motorglider         11-13
%   clean wind tunnel   10-12
%   average wind tunnel    9 <= standard "e^9 method"
%   dirty wind tunnel    4-8
%
% A flap deflection can be added using the following command,
% 'gdes flap {xhinge} {yhinge} {flap_defelction} exec'
%
% Outputs:
%  polar: structure with the polar coefficients (alpha,CL,CD,CDp,CM,
%          Top_Xtr,Bot_Xtr)
%   foil: stucture with the specific aoa values (s,x,y,UeVinf,
%          Dstar,Theta,Cf,H,cpx,cp) each column corresponds to a different
%          angle-of-attack.
%         If only one left hand operator is specified, only the polar will be parsed and output
%
% If there are different sized output arrays for the different incidence
% angles then they will be stored in a structured array, foil(1),foil(2)...
%
% If the output array does not have all alphas in it, that indicates a convergence failure in Xfoil.
% In that event, increase the iteration count with 'oper iter ##;
%
% Examples:
%    % Single AoA with a different number of panels
%    [pol foil] = xfoil('NACA0012',10,1e6,0.0,'panels n 330')
%
%    % Change the maximum number of iterations
%    [pol foil] = xfoil('NACA0012',5,1e6,0.2,'oper iter 50')
%
%    % Deflect the trailing edge by 20deg at 60% chord and run multiple incidence angles
%    [pol foil] = xfoil('NACA0012',[-5:15],1e6,0.2,'oper iter 150','gdes flap 0.6 0 5 exec')
%    % Deflect the trailing edge by 20deg at 60% chord and run multiple incidence angles and only
%    parse or output a polar.
%    pol = xfoil('NACA0012',[-5:15],1e6,0.2,'oper iter 150','gdes flap 0.6 0 5 exec')
%    % Plot the results
%    figure;
%    plot(pol.alpha,pol.CL); xlabel('alpha [\circ]'); ylabel('C_L'); title(pol.name);
%    figure; subplot(3,1,[1 2]);
%    plot(foil(1).xcp(:,end),foil(1).cp(:,end)); xlabel('x');
%    ylabel('C_p'); title(sprintf('%s @ %g\\circ',pol.name,foil(1).alpha(end)));
%    set(gca,'ydir','reverse');
%    subplot(3,1,3);
%    I = (foil(1).x(:,end)<=1);
%    plot(foil(1).x(I,end),foil(1).y(I,end)); xlabel('x');
%    ylabel('y'); axis('equal');
%
% Some default values
if ~exist('coord','var'), coord = 'NACA0012'; end
if ~exist('Alpha','var'), Alpha = 0; end
if ~exist('Re','var'),    Re = 1e6; end
if ~exist('Mach','var'),  Mach = 0; end

Nangle = length(Alpha); % number of alphas swept
[wd,fname,~] = fileparts(which(mfilename)); % working directory

% SPECIFY FOIL PROFILE ----------------------------------------------------
if ischar(coord)  % check for (1) NACA string or (2) foil coordinate file
    if isempty(regexpi(coord,'^NACA *[0-9]{4,5}$'))
        coord_file = coord;
    end
    
elseif ismatrix(coord) % check for (3) foil coordinate vector [n x 2]
    coord_file = [fname '.foil'];
    
    % gate 1: overwrite existing file?
    if exist(coord_file,'file')
        answer = questdlg(['Overwrite existing file: "',coord_file,'"?'],...
            'Warning', ...
            'NO','YES','NO');
        % handle response
        switch answer
            case 'NO'
                error('Exiting program.');
            case 'YES'
                delete(coord_file); %overwrite existing file
                disp('Existing file overwritten.');
        end
    end
    
    % gate 2: correct format?
    if size(coord,2) ~= 2
        if size(coord,1) == 2
            warndlg('Coordinate vector transposed. Verify format.','Input Warning');
            coord = coord'; % [2 x n]' = [n x 2]
        else
            error('Foil coordinate vector must be [nx2] or [2xn]. Exiting program.');
        end
    end
    
    % generate foil coordinate file:
    fid = fopen(coord_file,'w');
    if (fid<=0)
        error('Unable to create xfoil coordinate file (%s) ',coord_file);
    else
        fprintf(fid,'%s\n',fname);
        fprintf(fid,'%9.5f   %9.5f\n',coord');
        fclose(fid);
    end
else
    error([mfilename ':io'],'Unable to read foil coordinate input.');
end

% WRITE XFOIL COMMAND FILE ------------------------------------------------
cmd_file = [fname '.inp'];
fid = fopen([wd filesep cmd_file],'w');
if (fid<=0)
    error('Unable to create xfoil command file (%s)',cmd_file);
else
    % write: foil profile
    if ischar(coord)
        if ~isempty(regexpi(coord,'^NACA *[0-9]{4,5}$'))  % NACA string supplied
            fprintf(fid,'NACA %s\n',coord(5:end));
        else  % user-supplied coordinate file
            fprintf(fid,'LOAD %s\n',coord_file);
        end
    elseif ismatrix(coord) % user-supplied coordinate vector
        fprintf(fid,'LOAD %s\n',coord_file);
    end
    
    % write: blake's additions for modern computers
    fprintf(fid,'\n\nPPAR\n');
    fprintf(fid,'N %u\n',120);    % increase n panels
    fprintf(fid,'T %g\n\n\n',1);  % uniform LE/TE panel density
    
    % write: specify additional xfoil commands
    for ii = 1:length(varargin)
        txt = varargin{ii};
        txt = regexprep(txt,'[ \\\/]+','\n');
        fprintf(fid,'%s\n\n',txt);
    end
    
    % write: enter OPER submenu
    fprintf(fid,'\n\nOPER\n');
    % write: reynolds number and mach number (from function arguments)
    fprintf(fid,'RE %g\n',Re);
    fprintf(fid,'MACH %g\n',Mach);
    fprintf(fid,'ITER %u\n',500); % increase max iterations
    
    % write: switch to viscous mode
    if (Re>0)
        fprintf(fid,'VISC\n');
    end
    
    % write: polar accumulation mode enabled
    fprintf(fid,'pacc\n\n\n');
    
    % write: loop through alphas, create DUMP file & CPx file for each
    [file_dump, file_cpwr] = deal(cell(1,Nangle)); % preallocate cells
    
    for ii = 1:Nangle
        % assign output filenames to preallocated cells
        file_dump{ii} = sprintf('%s_a%06.3f_dump.dat',fname,Alpha(ii));
        file_cpwr{ii} = sprintf('%s_a%06.3f_cpwr.dat',fname,Alpha(ii));
        
        % write: run x-foil for each alpha, store in file
        fprintf(fid,'ALFA %g\n',Alpha(ii));
        fprintf(fid,'DUMP %s\n',file_dump{ii});
        fprintf(fid,'CPWR %s\n',file_cpwr{ii});
    end
    
    % write: polar output filename
    file_pwrt = sprintf('%s_pwrt.dat',fname);
    fprintf(fid,'PWRT\n%s\n',file_pwrt);
    fprintf(fid,'PLIS\n'); % why? perhaps this gets captured by xfoil.out
    fprintf(fid,'\nquit\n'); % quit xfoil
    fclose(fid); % close input file. ready to run.
end

% RUN: XFOIL WITH INPUT FILE (WRITTEN ABOVE) ------------------------------
xfoil_app = '/Applications/Xfoil.app/Contents/Resources/xfoil';
cmd = sprintf(['cd %s && ',xfoil_app,' < xfoil.inp > xfoil.out'],wd);
[status,result] = system(cmd);
if (status~=0)
    disp(result);
    error('Xfoil execution failed! CMD: %s',cmd);
end
% read: DUMP file from xfoil
%    #    s        x        y     Ue/Vinf    Dstar     Theta      Cf       H
jj = 0;
ind = 1;
foil.Alpha = zeros(1,Nangle); % Preallocate alphas

% Find the number of panels with an inital run
Nout = nargout; % Number of outputs checked. If only one left hand operator then only do polar
if Nout >1 % Only do the foil calculations if more than one left hand operator is specificed
    %WHY? what if user supplies: [foil] = xfoil()?
    for ii = 1:Nangle
        jj = jj + 1;
        fid = fopen([wd filesep file_dump{ii}],'r');
        if (fid<=0)
            error('Unable to read xfoil output file (%s)',file_dump{ii});
        else
            % store data
            D = textscan(fid,'%f%f%f%f%f%f%f%f','Delimiter',' ',...
                'MultipleDelimsAsOne',1,'CollectOutput',1,'HeaderLines',1);
            % close and delete file
            fclose(fid);
            delete([wd filesep file_dump{ii}]);
            
            if ii == 1
                Npanel = length(D{1}); % determine number of panels for first alpha
                % preallocate data structures
                [foil.s, foil.x, foil.y, ...
                    foil.UeVinf, foil.Dstar, foil.Theta,...
                    foil.Cf, foil.H] = deal(zeros(Npanel,Nangle));
            end
            
            % store data
            % note: not sure what this logic gate does.  seems weird.
            if ((jj>1) && (size(D{1},1)~=length(foil(ind).x))...
                    && sum(abs(foil(ind).x(:,1)-size(D{1},1)))>1e-6 )
                ind = ind + 1;
                jj = 1;
            end
            
            foil.s(:,jj) = D{1}(:,1);
            foil.x(:,jj) = D{1}(:,2);
            foil.y(:,jj) = D{1}(:,3);
            foil.UeVinf(:,jj) = D{1}(:,4);
            foil.Dstar(:,jj) = D{1}(:,5);
            foil.Theta(:,jj) = D{1}(:,6);
            foil.Cf(:,jj) = D{1}(:,7);
            foil.H (:,jj)= D{1}(:,8);
        end
        foil.Alpha(1,jj) = Alpha(jj);
        
        % read: CPx file from xfoil
        fid = fopen([wd filesep file_cpwr{ii}],'r');
        
        if (fid<=0)
            error('Unable to read xfoil output file (%s)',file_cpwr{ii});
        else
            % store all data in matlab array
            C = textscan(fid, '%10f%10f', 'Delimiter', '', 'WhiteSpace',...
                '', 'ReturnOnError', false, 'HeaderLines' ,1);
            %NOTE: WhiteSpace spec might be unneccesary
            % close and delete file
            fclose(fid);
            delete([wd filesep file_cpwr{ii}]);
            
            % store data
            if ii == 1
                NCp = length(C{1});
                % preallocate data structures
                [foil.xcp, foil.cp] = deal(zeros(NCp,Nangle));
                foil.xcp = C{1}(:,1);
            end
            foil.cp(:,jj) = C{2}(:,1);
        end
    end
end

if Nout <= 1% clear files for default run
    % NOTE: why is this not elseif? Nout is either <=1 or >1
    
    for ii=1:Nangle % Clear out the xfoil dump files not used
        delete([wd filesep file_dump{ii}]);
        delete([wd filesep file_cpwr{ii}]);
    end
end
% Read polar file
%
%       XFOIL         Version 6.97
%
% Calculated polar for: NACA 0012
%
% 1 1 Reynolds number fixed          Mach number fixed
%
% xtrf =   1.000 (top)        1.000 (bottom)
% Mach =   0.000     Re =     1.000 e 6     Ncrit =  9.000
%
%   alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr
%  ------ -------- --------- --------- -------- -------- --------
fid = fopen([wd filesep file_pwrt],'r');
if (fid<=0)
    error('Unable to read xfoil polar file (%s)',file_pwrt);
else
    % Header
    % Calculated polar for: NACA 0012
    P = textscan(fid,' Calculated polar for: %[^\n]','Delimiter',' ',...
        'MultipleDelimsAsOne',true,'HeaderLines',3);
    pol.name = strtrim(P{1}{1});
    % xtrf =   1.000 (top)        1.000 (bottom)
    P = textscan(fid, '%*s%*s%f%*s%f%s%s%s%s%s%s', 1, 'Delimiter', ' '...
        , 'MultipleDelimsAsOne', true, 'HeaderLines', 2, 'ReturnOnError', false);
    %NOTE: includes lots of strings after 1.0000 (bottom)?
    %NOTE: header lines = 3?
    pol.xtrf_top = P{1}(1);
    pol.xtrf_bot = P{2}(1);
    % Mach =   0.000     Re =     1.000 e 6     Ncrit =  12.000
    P = textscan(fid, '%*s%*s%f%*s%*s%f%*s%f%*s%*s%f', 1, 'Delimiter', ' ',...
        'MultipleDelimsAsOne', true, 'HeaderLines', 0, 'ReturnOnError', false);
    pol.Mach = P{1}(1);
    pol.Re = P{2}(1) * 10^P{3}(1);
    pol.Ncrit = P{4}(1);
    % data
    P = textscan(fid, '%f%f%f%f%f%f%f%*s%*s%*s%*s', 'Delimiter',  ' ',...
        'MultipleDelimsAsOne', true, 'HeaderLines' , 4, 'ReturnOnError', false);
    % NOTE: header lines = 3?
    fclose(fid);
    %delete([wd filesep file_pwrt]);
    % store data
    pol.Alpha = P{1}(:,1);
    pol.CL  = P{2}(:,1);
    pol.CD  = P{3}(:,1);
    pol.LD  = (pol.CL)./(pol.CD);
    pol.CDp = P{4}(:,1);
    pol.Cm  = P{5}(:,1);
    pol.Top_xtr = P{6}(:,1);
    pol.Bot_xtr = P{7}(:,1);
end
if length(pol.Alpha) ~= Nangle % Check if xfoil failed to converge
    warning(['One or more alpha values failed to converge. ',...
        'Recorded as NaN.']);
    
    %fill interior gaps with NaN
    dAlpha = diff(Alpha);
    gaps = find(diff(pol.Alpha)~=dAlpha(1));
    
    for i=1:length(gaps)
        pol.Alpha = [pol.Alpha(1:gaps(i));NaN;pol.Alpha(gaps(i)+1:end)];
        pol.CL = [pol.CL(1:gaps(i));NaN;pol.CL(gaps(i)+1:end)];
        pol.CD = [pol.CD(1:gaps(i));NaN;pol.CD(gaps(i)+1:end)];
        pol.LD = [pol.LD(1:gaps(i));NaN;pol.LD(gaps(i)+1:end)];
        pol.CDp = [pol.CDp(1:gaps(i));NaN;pol.CDp(gaps(i)+1:end)];
        pol.Cm = [pol.Cm(1:gaps(i));NaN;pol.Cm(gaps(i)+1:end)];
        pol.Top_xtr = [pol.Top_xtr(1:gaps(i));NaN;pol.Top_xtr(gaps(i)+1:end)];
        pol.Bot_xtr = [pol.Bot_xtr(1:gaps(i));NaN;pol.Bot_xtr(gaps(i)+1:end)];
        
        gaps = gaps+1; % compensate for inserted NaN
    end
    
    %fill bounding gaps with NaN
    if (Nangle-length(pol.Alpha)==2)    % front and back missing
        pol.Alpha = [NaN;pol.Alpha;NaN];
        pol.CL = [NaN;pol.CL;NaN];
        pol.CD = [NaN;pol.CD;NaN];
        pol.LD = [NaN;pol.LD;NaN];
        pol.CDp = [NaN;pol.CDp;NaN];
        pol.Cm = [NaN;pol.Cm;NaN];
        pol.Top_xtr = [NaN;pol.Top_xtr;NaN];
        pol.Bot_xtr = [NaN;pol.Bot_xtr;NaN];
        
    elseif (Nangle-length(pol.Alpha)==1)    %either front or back missing
        if pol.Alpha(1) ~= Alpha(1)
            pol.Alpha = [NaN;pol.Alpha];
            pol.CL = [NaN;pol.CL];
            pol.CD = [NaN;pol.CD];
            pol.LD = [NaN;pol.LD];
            pol.CDp = [NaN;pol.CDp];
            pol.Cm = [NaN;pol.Cm];
            pol.Top_xtr = [NaN;pol.Top_xtr];
            pol.Bot_xtr = [NaN;pol.Bot_xtr];
        else
            pol.Alpha = [pol.Alpha;NaN];
            pol.CL = [pol.CL;NaN];
            pol.CD = [pol.CD;NaN];
            pol.LD = [pol.LD;NaN];
            pol.CDp = [pol.CDp;NaN];
            pol.Cm = [pol.Cm;NaN];
            pol.Top_xtr = [pol.Top_xtr;NaN];
            pol.Bot_xtr = [pol.Bot_xtr;NaN]; 
        end
    end
end

end
