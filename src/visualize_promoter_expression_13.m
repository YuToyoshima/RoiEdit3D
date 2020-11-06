function visualize_promoter_expression_13(varargin)
%%% This script can visualize promoter expression patterns with 3D graphics
%%%
%%% This script can throw a query to wormbase.org in order to obtain
%%% expression pattern using geneIDs.csv.
%%% Please see the instruction on following web site to renew geneIDs.csv.
%%% https://github.com/toricor/Research/tree/master/geneID_to_expression_pattern
%%% 
%%% This script uses additional files.
%%% Please compile as 
%%% mcc('-m','visualize_promoter_expression_12.m',...
%%%     '-a','MyTableModelVisPromExp.class',...
%%%     '-a','merge_emtrace_white.mat',...
%%%     '-a','visualize_promoter_expression_param.txt',...
%%%     '-a','geneIDs.csv'); 
%%%
%%% Written by Yu Toyoshima, 2017/7/10
%%%
%%% Version
%%% v6: first version
%%% v7: add promoter search functionality
%%% v8: add display roiedit3d result (Drug & Drop to roiedit3d.mat to exe)
%%% v9: extend number of markup tabs more than 6. 
%%%     scaling for roiedit3d result.
%%%     disable TeX interpreter for roi name (for underscore)
%%% v10: format of query result from wormbase was changed.
%%% v11: add cooporation to roiedit3d (for auto annotaion)
%%% v12: add median positions of neurons in the Iinolab annotation dataset
%%%      (as one of the reference position data set)
%%%      add editor function
%%% v13: format of query result from wormbase was changed.
%%%      fix a bug for saving the edited result
%%%
%%% TODO
%%%  open & close control tab (to enlarge each panel)
%%%  editor function (save base & mark)
%%%  columns in mark table and display data origin (paper or experiment)
%%%  (add & remove columns)
%%%  search ROIs (highlighting)
%%%  

%%% specify data file name
strdatamat = 'visualize_promoter_expression.mat';
stropttxt  = 'visualize_promoter_expression_param.txt';



%%% load data


%%% load default data file
pathdatamat = fullfile(fileparts(mfilename('fullpath')),strdatamat);
disp(['Search data file (global):',pathdatamat]);
if numel(dir(pathdatamat))==1 % an option file found
    disp('Load data file (global)...');
    l = load(pathdatamat,'pos_base','list_base','names_base','list_markup','names_markup');
else 
    disp('Data file (global) not found');
end

%%% load local data file
pathdatamat = fullfile(pwd,strdatamat);
disp(['Search data file (local):',pathdatamat]);
if numel(dir(pathdatamat))==1 % an option file found
    disp('Load data file (local)...');
    l = load(pathdatamat,'pos_base','list_base','names_base','list_markup','names_markup');
else
    disp('Data file (local) not found');
end

% 
% pos = l.pos;
% names = l.names;
% pos_pha = l.pos_pha;
% names_pha = l.names_pha;
% pos_white = l.pos_white;
% names_white = l.names_white;
% pos_emt = l.pos_emt;
% names_emt = l.names_emt;


% %%% set expression pattern
% exp_dye = {'ASK[LR]'; 'ADL[LR]'; 'ASI[LR]'; 'AWB[LR]'; 'ASH[LR]'; 'ASJ[LR]'; 'PHA[LR]'; 'PHB[LR]';};
% 
% 
% %%% filtering input
% %%% already loaded variables: 
% %%% pos, pos_white, pos_pha, pos_emt, names, names_white, names_pha, and names_emt 
% flag_h20_head = pos(:,1)<pos(strcmp(names,'AVG'),1);
% pos_h20_head = pos(flag_h20_head,:);
% names_h20_head = names(flag_h20_head);
% pos_h20_all = pos;
% names_h20_all = names;
% pos_h20_pha = pos_pha;
% names_h20_pha = names_pha;
% flag_white_head = pos_white(:,1)<pos_white(strcmp(names_white,'AVG'),1);
% pos_white_head = pos_white(flag_white_head,:);
% names_white_head = names_white(flag_white_head);
% pos_white_all = pos_white;
% names_white_all = names_white;
% flag_white_pha = false(size(names_white));
% for p=1:numel(names_white)
%     flag_white_pha(p) = any(strcmp(names_pha,names_white(p)));
% end
% pos_white_pha = pos_white(flag_white_pha,:);
% names_white_pha = names_white(flag_white_pha);


%%% prepare data
strgidcsv = 'geneIDs.csv';
strregexp = '( neuron| nucleus)$';
strrep = '';
numMarkDisp = 30;
strAxesLabel = {'A/P','D/V','L/R'};

% pos_base = [{pos_h20_head},{pos_h20_all},{pos_h20_pha},...
%     {pos_white_head},{pos_white_all},{pos_white_pha},{pos_emt}];
% list_base = {'H20_Head','H20_All','H20_Pharynx','White_Head','White_All','White_Pharynx','EM_Trace'};
% names_base = {names_h20_head,names_h20_all,names_h20_pha,names_white_head,names_white_all,names_white_pha,names_emt};
% list_markup = {'dye','  '};
% names_markup = {exp_dye,{'  '}};
% 
% pos_base = [{pos_h20_head},{pos_h20_all},{pos_h20_pha},{pos_white_head},{pos_white_all},{pos_white_pha},{pos_emt},{zeros(1,3)}];
% list_base = {'H20_Head','H20_All','H20_Pharynx','White_Head','White_All','White_Pharynx','EM_Trace','-'};
% names_base = {names_h20_head,names_h20_all,names_h20_pha,names_white_head,names_white_all,names_white_pha,names_emt,{'-'}};
% list_markup = {'eat-4p','tax-4p','unc-17/cha-1','cho-1p','glr-1p','ser-2p2','glr-1p+ser-2p2','flp12p','dye','dye_ADF','dye_7','asynmetric','ciliated','pharynx','-'};
% names_markup = {exp_eat4,exp_tax4,exp_unc17,exp_cho1,exp_glr1,exp_ser2p2,exp_glr1ser2p2,exp_flp12p,exp_dye,exp_dye_ADF,exp_dye_7,exp_asynmetric,ciliated,exp_pharynx,{'-'}};
%
% d3 = load('D:\download\documents\180918\AnnoDB_paper\elements\test181013_ellipsoid_expression3.mat');
% flag_confocal_head = d3.posPcaMedian(:,1)<d3.posPcaMedian(strcmp(d3.names,'AVG'),1);
% pos_confocal_head = d3.posPcaMedian(flag_confocal_head,:);
% names_confocal_head = d3.names(flag_confocal_head);
% pos_confocal_pha = d3.posPcaMedian(d3.flagPha,:);
% names_confocal_pha = d3.names(d3.flagPha);
% pos_base = [{pos_confocal_head},{d3.posPcaMedian},{pos_confocal_pha},{pos_h20_head},{pos_h20_all},{pos_h20_pha},{pos_white_head},{pos_white_all},{pos_white_pha},{pos_emt}];
% list_base = {'H20_Head','H20_All','H20_pharynx','H20(EM)_Head','H20(EM)_All','H20(EM)_Pharynx','White_Head','White_All','White_Pharynx','EM_Trace'};
% names_base = {names_confocal_head,d3.names,names_confocal_pha,names_h20_head,names_h20_all,names_h20_pha,names_white_head,names_white_all,names_white_pha,names_emt};

pos_base = l.pos_base;
list_base = l.list_base;
names_base = l.names_base;
list_markup = l.list_markup;
names_markup = l.names_markup;

limUIRange = [eps,10,50; eps,6,50];
strColumnNamesBase = {'Name','x','y','z'};
strColumnNamesMark = {'Name'};
strClassifier = 'Promoter';
str_clmap = {'Black','Red','Green','Blue','Cyan','Magenta','Yellow'};
mat_clmap = [
    0,   0,   0;
    1,   0,   0;
    0,   0.8, 0;
    0,   0,   1;
    0,   0.5, 0.7;
    0.7, 0,   0.7;
    0.7, 0.5, 0; 
    ];


%%% load default option file
flag_option_global_found = false;
pathopttxt = fullfile(fileparts(mfilename('fullpath')),stropttxt);
disp(['Search option file (global):',pathopttxt]);
if numel(dir(pathopttxt))==1 % an option file found
    flag_option_global_found = true;
    disp('Load option file (global)...')
    fid = fopen(pathopttxt,'rt');
    eval(fread(fid,inf,'*char'));
    fclose(fid);
else 
    disp('Option file (global) not found');
end

%%% load local option file
flag_option_local_found = false;
pathopttxt = fullfile(pwd,stropttxt);
disp(['Search option file (local):',pathopttxt]);
if numel(dir(pathopttxt))==1 % an option file found
    flag_option_local_found = true;
    disp('load option file (local)...')
    fid = fopen(pathopttxt,'rt');
    eval(fread(fid,inf,'*char'));
    fclose(fid);
else
    disp('Option file (local) not found');
end

if (~flag_option_global_found && ~flag_option_local_found)
    disp('Use default option ...');
end



%%% load default geneID file
pathgidcsv = fullfile(fileparts(mfilename('fullpath')),strgidcsv);
disp(['Search geneID file (global):',pathgidcsv]);
if numel(dir(pathgidcsv))==1 % an option file found
    disp('Load geneID file (global)...');
    geneids = readtable(pathgidcsv,'ReadVariableNames',false);
else 
    disp('GeneIDs file (global) not found');
end

%%% load local geneID file
pathgidcsv = fullfile(pwd,strgidcsv);
disp(['Search options file (local):',pathgidcsv]);
if numel(dir(pathgidcsv))==1 % an option file found
    disp('Load geneID file (local)...');
    geneids = readtable(pathgidcsv,'ReadVariableNames',false);
else
    disp('GeneIDs file (local) not found');
end


%%% load roiedit3d mat file
frame_target = 1;
idxCol.check = 1;
idxCol.Name = 3;
idxCol.frame = 4;
idxCol.mu = 10:12;
for p=nargin:-1:1
    [~,filename] = fileparts(varargin{p});
    disp(['loading roiedit3d mat file: ',varargin{p}]);
    r = load(varargin{p});
    flag_target = cell2mat(r.data_edited(:,idxCol.frame))==frame_target;
    flag_check = cell2mat(r.data_edited(flag_target,idxCol.check));
    name_roiedit3d = cellfun(@char,r.data_edited(flag_target,idxCol.Name),...
        'UniformOutput',false);    
%     pos_roiedit3d = bsxfun(@cell2mat(r.data_edited(flag_target,idxCol.mu));
    pos_roiedit3d = bsxfun(@times,r.scaling,cell2mat(r.data_edited(flag_target,idxCol.mu)));
    pos_base = [{pos_roiedit3d},pos_base];
    list_base = [{filename},list_base];
    names_base = [{name_roiedit3d},names_base];
    list_markup = [{['checked_',filename]},list_markup];
    names_markup = [{name_roiedit3d(flag_check)},names_markup];
end


%%% create UI


% %%% add java classpath for 'MyTabelModel.class'
classpath = javaclasspath('-all');
if ~exist('MIJ','class') && isempty(cell2mat(regexp(classpath,pwd)))
    javaaddpath(pwd,'-end');
end


%%% make figure
hfig = figure(...
    'Tag','figure_vis_prom_exp',...
    'WindowStyle','normal',...    
    'Position',[630,100,900,450],...
    'MenuBar','none',...
    'NumberTitle','off',...
    'Name','Visualize Promoter Expression');
% hfig = figure('WindowScrollWheelFcn',@cbWheel);
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(handle(hfig),'JavaFrame');
hlistener_figure = addlistener(jFrame.fHG2Client.getAxisComponent,'MouseWheelMoved',{@cbWheel,hfig});
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

hpanel(1) = uipanel('Parent',hfig,'Title','Control','Position',[0,0,0.25,1]);
hpanel(2) = uipanel('Parent',hfig,'Title','View','Position',[0.25,0,0.75,1]);
htabgp = uitabgroup(hpanel(1),'TabLocation','left');


%%% plot
haxis = gca(hpanel(2));
[hmarker,hlabel] = plotfun(haxis,[0,0,0],{'dummy'},strAxesLabel);


%%% create size contoroller
htab(1) = uitab(htabgp,'Title','Size');

strUIrange = strAxesLabel;
limUIrange = [
    get(haxis,'xlim');
    get(haxis,'ylim');
    get(haxis,'zlim')
    ];
for p=1:3
    [htext(2*p-1),hedit(2*p+(-1:0)),hslider(2*p+(-1:0))] = ...
        createControls2(htab(1),1-0.2*p,0.05); %#ok<AGROW>
    set(htext(2*p-1),'String',strUIrange{p});
    set(hslider(2*p+(-1:0)),'Min',limUIrange(p,1),'Max',limUIrange(p,2));
    set(hslider(2*p-1),'Value',limUIrange(p,1));
    set(hslider(2*p),  'Value',limUIrange(p,2));
    set(hedit(2*p-1),'String',num2str(limUIrange(p,1)));
    set(hedit(2*p),  'String',num2str(limUIrange(p,2)));
end
set(hslider(1:6),'Callback',@cbRange);
hlistener(1:6) = addlistener(hslider(1:6),'ContinuousValueChange',@cbRange);
set(hedit(1:6),'Callback',@cbRange);
 
strUISize = {'Marker Size','Label Size'};
funUISize = {@cbSizeMarker,@cbSizeLabel};
tmpoffset = 6;         
for p=1:2
    [htext(p+tmpoffset),hedit(p+tmpoffset),hslider(p+tmpoffset)] = ...
        createControls(htab(1),0.35-(0.1+0.05)*p,0.05);     
    set(htext(p+tmpoffset),'String',strUISize{p});
    set(hedit(p+tmpoffset),...
        'String',num2str(limUIRange(p,2)),...
        'Callback',funUISize{p});
    set(hslider(p+tmpoffset),...
        'Min',  limUIRange(p,1),...
        'Value',limUIRange(p,2),...
        'Max',  limUIRange(p,3),...
        'Callback',funUISize{p});
    hlistener(p+tmpoffset) = addlistener(hslider(p+tmpoffset),...
        'ContinuousValueChange',funUISize{p});
end


%%% create base controller
htab(2) = uitab(htabgp,'Title','Base');
tmpdata = cat(2,names_base{1},num2cell(pos_base{1}));
tmprange = 9:10;
[htext(tmprange),hpopup(tmprange),hedit(tmprange),hbutton(tmprange),...
    hmtm(1),hstm(1),hjtable(1)] ...
    = createControls3(htab(2),tmpdata,strColumnNamesBase);
set(htext(9),'String',strClassifier);
set(hpopup(9),...
    'String',list_base,...
    'Value',1,...
    'Callback',{@cbBasePopup,hmtm(1),hedit(9)});
set(hedit(9),...
    'String',list_base(1),...
    'CallBack',@setTitle);
set(hpopup(10),...
    'String',str_clmap,...
    'Value',1,...
    'Callback',@cbMarkTableData);
set(hbutton(9),'Enable','off')
set(hbutton(10),'Callback',@(obj,ev)hmtm(1).addRow({'',0,0,0}));
set(findobj('Tag','button_vpe_Base_save'),'Callback',@(~,~)saveBase(hmtm(1),hedit(9),pathopttxt));
hlistener_jtable(1) ...
    = addlistener(hjtable(1),'KeyReleased',{@cbTableKey,hmtm(1)});
hlistener_mtm(1) ...
    = addlistener(hmtm(1),'TableChanged',{@cbBaseTableData,hmtm(1)});


%%% create Mark controller
tmpv = numel(list_markup)*ones(1,numMarkDisp);
tmpv(1) = 1;
numcolor = size(mat_clmap,1);
for p=1:numMarkDisp
    offset = 8+2*p;
    tmprange = (1:2) + offset;
    htab(p+2) = uitab(htabgp,'Title',['Mark',num2str(p)]);
    [htext(tmprange),hpopup(tmprange),hedit(tmprange),hbutton(tmprange),...
        hmtm(p+1),hstm(p+1),hjtable(p+1)] ...
        = createControls3(htab(p+2),names_markup{tmpv(p)},strColumnNamesMark); %#ok<AGROW>
    set(htext(offset+1),'String',strClassifier);   
    set(hpopup(offset+1),...
        'String',list_markup,...
        'Value',tmpv(p),...
        'Callback',{@cbMarkPopup,hmtm(p+1),hedit(offset+1)});
    set(hedit(offset+1),...
        'Tag',['edit_subset_name_',num2str(p)],...
        'String',list_markup{tmpv(p)},...
        'Callback',@setTitle);
    set(hpopup(offset+2),...
        'String',str_clmap,...
        'Value',mod(p,numcolor)+1,...
        'Callback',@cbMarkTableData);
    set(hbutton(offset+1),'Callback',{@cbSearch,hmtm(p+1),hedit(offset+1)});
    set(hbutton(offset+2),'Callback',@(obj,ev)hmtm(p+1).addRow({''}));
    set(findobj('Tag',['button_vpe_Mark',num2str(p),'_save']),'Callback',@(~,~)saveMark(hmtm(p+1),hedit(offset+1),pathopttxt));
    hlistener_jtable(p+1) ...
        = addlistener(hjtable(p+1),'KeyReleased',{@cbTableKey,hmtm(p+1)}); %#ok<AGROW>
    hlistener_mtm(p+1) ...
        = addlistener(hmtm(p+1),'TableChanged',@cbMarkTableData); %#ok<AGROW>
end


%%% create title
for p=1:numel(hjtable)
    htitle(p) = uicontrol(hpanel(2),...
        'Style','text',...
        'Unit','Normalized',...
        'Position',[0.9,1-0.03*p,0.1,0.03]); %#ok<AGROW>
end


%%% hold handles
userdata.pos_base = pos_base;
userdata.list_base = list_base;
userdata.list_markup = list_markup;
userdata.names_base = names_base;
userdata.names_markup = names_markup;
userdata.geneids = geneids;
userdata.strregexp = strregexp;
userdata.strrep = strrep;
userdata.strColumnNamesBase = strColumnNamesBase;
userdata.strColumnNamesMark = strColumnNamesMark;
userdata.hfig = hfig;
userdata.htabgp = htabgp;
userdata.htab = htab;
userdata.hpanel = hpanel;
userdata.haxis = haxis;
userdata.hmarker = hmarker;
userdata.hlabel = hlabel;
userdata.strAxesLabel = strAxesLabel;
userdata.htext = htext;
userdata.mat_clmap = mat_clmap;
userdata.hedit = hedit;
userdata.hslider = hslider;
userdata.hpopup = hpopup;
userdata.hbutton = hbutton;
userdata.hmtm = hmtm;
userdata.hstm = hstm;
userdata.hjtable = hjtable;
userdata.hlistener = hlistener;
userdata.hlistener_figure = hlistener_figure;
userdata.hlistener_jtable = hlistener_jtable;
userdata.hlistener_mtm = hlistener_mtm;
userdata.htable_semaphore = java.util.concurrent.Semaphore(1,true);
userdata.htitle = htitle;
set(gcf,'UserData',userdata);

cbBaseTableData([],[],hmtm(1));

end


function saveBase(hmtm,hedit,pathopttxt)

names = arrayfun(@char,hmtm.getDataColumn(0).toArray(),'UniformOutput',false);
pos = cat(2,...
    arrayfun(@double,hmtm.getDataColumn(1).toArray()),...
    arrayfun(@double,hmtm.getDataColumn(2).toArray()),...
    arrayfun(@double,hmtm.getDataColumn(3).toArray()));
strLabel = hedit.String;

fid = fopen(pathopttxt,'at+');

fprintf(fid,['\n\r\n\rnames_',strLabel,'={\n\r']);
fprintf(fid,'''%s'';\n\r',names{:});
fprintf(fid,'};\n\r\n\r');

fprintf(fid,['pos_',strLabel,'=[\n\r']);
fprintf(fid,'%f,%f,%f;\n\r',pos');
fprintf(fid,'];\n\r\n\r');

fprintf(fid,['names_base = cat(2,names_base(1:end-1),{names_',strLabel,'},names_base(end));\n\r']);
fprintf(fid,['pos_base = cat(2,pos_base(1:end-1),{pos_',strLabel,'},pos_base(end));\n\r']);
fprintf(fid,['list_base = cat(2,list_base(1:end-1),{''',strLabel,'''},list_base(end));\n\r\n\r']);

fclose(fid);

disp(['Base information was saved to ',pathopttxt]);

end


function saveMark(hmtm,hedit,pathopttxt)

names = arrayfun(@char,hmtm.getDataColumn(0).toArray(),'UniformOutput',false);
strLabel = hedit.String;
strLabelValidVar = matlab.lang.makeValidName(strLabel);

fid = fopen(pathopttxt,'at+');

% fprintf(fid,['\n\r\n\rnames_',strLabel,'={\n\r']);
fprintf(fid,['\n\r\n\rnames_',strLabelValidVar,'={\n\r']);
fprintf(fid,'''%s'';\n\r',names{:});
fprintf(fid,'};\n\r\n\r');

% fprintf(fid,['names_markup = cat(2,names_markup(1:end-1),{names_',strLabel,'},names_markup(end));\n\r']);
fprintf(fid,['names_markup = cat(2,names_markup(1:end-1),{names_',strLabelValidVar,'},names_markup(end));\n\r']);
fprintf(fid,['list_markup = cat(2,list_markup(1:end-1),{''',strLabel,'''},list_markup(end));\n\r\n\r']);

fclose(fid);
disp(['Mark information was saved to ',pathopttxt]);

end


function [hmarker,hlabel] = plotfun(haxis,pos,names,strAxesLabel)
axes(haxis);
hfig = get(get(haxis,'Parent'),'parent');
set(haxis,'nextplot','replace');
for p=1:numel(names)
    hmarker(p) = plot3(...
        pos(p,1),...
        pos(p,2),...
        pos(p,3),...
        'Marker','o',...
        'MarkerFaceColor','w',...
        'MarkerSize',10); %#ok<AGROW>
    hlabel(p) = text(...
        pos(p,1),...
        pos(p,2),...
        pos(p,3),...
        names(p),...
        'FontSize',6,...
        'Clipping','on',...
        'Interpreter','none'); %#ok<AGROW>
    set(haxis,'nextplot','add');
end
axis tight equal vis3d;
box on;
grid on;
xlabel(strAxesLabel{1});
ylabel(strAxesLabel{2});
zlabel(strAxesLabel{3});
camproj(haxis,'perspective');
cameratoolbar(hfig,'Show');
cameratoolbar(hfig,'SetMode','orbit');
cameratoolbar(hfig,'SetCoordSys','none');
set(haxis,...
    'CameraTargetMode','auto',...
    'CameraPosition',get(gca,'CameraPosition').*[1,1,-1],...
    'CameraUpVector',[0,-1,0]);
end


function cbWheel(~,ev,hfig)
userdata = get(hfig,'UserData');
haxis = userdata.haxis;
cp = haxis.CameraPosition;
ct = haxis.CameraTarget;
cuv = haxis.CameraUpVector;
cva = haxis.CameraViewAngle;

modifier = ev.getModifiersEx();
mask_ctrl  = java.awt.event.MouseEvent.CTRL_DOWN_MASK;
mask_shift = java.awt.event.MouseEvent.SHIFT_DOWN_MASK;
flag_ctrl  = bitand(modifier,mask_ctrl )==mask_ctrl;
flag_shift = bitand(modifier,mask_shift)==mask_shift;
scroll = ev.getPreciseWheelRotation()*ev.getScrollAmount();

if flag_ctrl % horizontally scrolling
    transloc = cross(ct-cp,cuv)*cva*0.0001*scroll;
    haxis.CameraPosition = cp + transloc;
    haxis.CameraTarget   = ct + transloc;
elseif flag_shift % vertically scrolling
    transloc = sqrt(sum((ct-cp).^2))*cuv*cva*0.0001*scroll;
    haxis.CameraPosition = cp + transloc;
    haxis.CameraTarget   = ct + transloc;
else % zoom 
    haxis.CameraViewAngle = cva*(1+0.01*scroll);
end

end


function cbRange(obj,~)
userdata = get(gcf,'UserData');
switch obj.Style
    case {'slider'}
        objlist = userdata.hslider(1:6);
        values = cell2mat(get(objlist,'Value'));
    case {'edit'}
        objlist = userdata.hedit(1:6);
        values = str2num(char(get(objlist,'String'))); %#ok<ST2NM>
end
objlist = reshape(objlist,[2,3]);
values = reshape(values,[2,3]);

%%% if min>=max
flag_inv = diff(values)<=0;
if any(flag_inv)
    idx_xyz = find(flag_inv);
    flag_minmax = objlist(:,idx_xyz)==obj;    
    values(~flag_minmax,idx_xyz) = values(flag_minmax,idx_xyz);
    values(1,idx_xyz) = values(1,idx_xyz)-abs(values(1,idx_xyz))*eps;
    values(2,idx_xyz) = values(2,idx_xyz)+abs(values(1,idx_xyz))*eps;
end

%%% force input in min-max
for p=1:6
    minvalue = get(objlist(p),'Min');
    maxvalue = get(objlist(p),'Max');
    values(p) = max(min(values(p),maxvalue),minvalue);
end

%%% update components
for p=1:6
    set(userdata.hslider(p),'Value',values(p));    
    set(userdata.hedit(p),'String',num2str(values(p)));    
end

%%% update view
set(userdata.haxis,'xlim',values(1:2),'ylim',values(3:4),'zlim',values(5:6));

end


function cbSizeMarker(obj,~)
userdata = get(gcf,'UserData');
switch obj.Style
    case {'slider'}
        value = get(obj,'Value');
    case {'edit'}
        value = str2double(char(get(obj,'String')));
end
minvalue = get(userdata.hslider(7),'Min');
maxvalue = get(userdata.hslider(7),'Max');
value = max(min(value,maxvalue),minvalue);
set(userdata.hslider(7),'Value',value);
set(userdata.hedit(7),'String',num2str(value));
set(userdata.hmarker,'MarkerSize',value);
end


function cbSizeLabel(obj,~)
userdata = get(gcf,'UserData');
switch obj.Style
    case {'slider'}
        value = get(obj,'Value');
    case {'edit'}
        value = str2double(char(get(obj,'String')));
end
minvalue = get(userdata.hslider(8),'Min');
maxvalue = get(userdata.hslider(8),'Max');
value = max(min(value,maxvalue),minvalue);
set(userdata.hslider(8),'Value',value);
set(userdata.hedit(8),'String',num2str(value));
set(userdata.hlabel,'FontSize',value);
end


function cbTableKey(obj,ev,hmtm)

% userdata = get(gcf,'UserData');

% if userdata.htable_semaphore.availablePermits() == 0
%     return;
% end
% userdata.htable_semaphore.acquire();

if ev.getKeyCode() == ev.VK_DELETE % delete table rows
    
    hjtable = obj;    
    selections = hjtable.getSelectedRows()';
    
    modelRows = zeros(size(selections));
    for p=1:numel(selections)
        modelRows(p) = com.jidesoft.grid.TableModelWrapperUtils.getActualRowAt(...
            hjtable.getModel,selections(p));
    end
    modelRowsSorted = sort(modelRows,'Descend');
    for p=1:numel(modelRowsSorted)
            hmtm.removeRow(modelRowsSorted(p));            
    end
    disp('Rows were removed.');
        
end

% if userdata.htable_semaphore.availablePermits() == 0
%     userdata.htable_semaphore.release();
% end

end


function cbMarkTableData(~,~)
userdata = get(gcf,'UserData');
hmtm_base = userdata.hmtm(1);
numEntry = hmtm_base.getRowCount();
numTable = numel(userdata.hmtm);
clmap_mark = zeros(numEntry,3);
numcolor = size(userdata.mat_clmap,1);
for p=1:numTable
    mat_clmap = userdata.mat_clmap(mod(userdata.hpopup(p*2+8).Value-1,numcolor)+1,:);
    flag_mark = hmtm_base.regexp(userdata.hmtm(p).getDataColumn(0),0)>=0;
    clmap_mark(flag_mark,:) = bsxfun(@plus,clmap_mark(flag_mark,:),mat_clmap);
end
clmap_mark = min(max(clmap_mark,0),1);

for p=1:numEntry
    set(userdata.hmarker(p),'MarkerEdgeColor',clmap_mark(p,:));
    set(userdata.hlabel(p),'Color',clmap_mark(p,:));
end

setTitle();
end


function cbBaseTableData(~,~,hmtm)
userdata = get(gcf,'UserData');
names = arrayfun(@char,hmtm.getDataColumn(0).toArray(),'UniformOutput',false);
pos = cat(2,...
    arrayfun(@double,hmtm.getDataColumn(1).toArray()),...
    arrayfun(@double,hmtm.getDataColumn(2).toArray()),...
    arrayfun(@double,hmtm.getDataColumn(3).toArray()));
[userdata.hmarker,userdata.hlabel] ...
    = plotfun(userdata.haxis,pos,names,userdata.strAxesLabel);
limUIrange = [
    get(userdata.haxis,'xlim');
    get(userdata.haxis,'ylim');
    get(userdata.haxis,'zlim')
    ];
warning('off','MATLAB:hg:uicontrol:ValueMustBeInRange');
for p=1:3    
    set(userdata.hslider(2*p+(-1:0)),'Min',limUIrange(p,1),'Max',limUIrange(p,2));
    set(userdata.hslider(2*p-1),'Value',limUIrange(p,1));
    set(userdata.hslider(2*p),  'Value',limUIrange(p,2));
end
warning('on','MATLAB:hg:uicontrol:ValueMustBeInRange');
cbRange(userdata.hslider(1),[]);
set(gcf,'UserData',userdata);
cbMarkTableData();
end


function cbMarkPopup(obj,~,hmtm,hedit)
userdata = get(gcf,'UserData');
set(hedit,'String',userdata.list_markup{obj.Value});
hmtm.setDataVector(userdata.names_markup{obj.Value},userdata.strColumnNamesMark);
setTitle();
end


function cbBasePopup(obj,~,hmtm,hedit)
userdata = get(gcf,'UserData');
set(hedit,'String',userdata.list_base{obj.Value});
tmpdata = cat(2,userdata.names_base{obj.Value},num2cell(userdata.pos_base{obj.Value}));
hmtm.setDataVector(tmpdata,userdata.strColumnNamesBase);
setTitle();
end


function cbSearch(~,~,hmtm,hedit)
userdata = get(gcf,'UserData');
geneids = userdata.geneids;
gid = geneids.Var1(strcmpi(hedit.String,geneids.Var2));
if isempty(gid); return; end
strquery = ['http://www.wormbase.org/rest/widget/gene/',gid{1},'/expression'];
disp(['Search query: ' strquery])
hobjson = webread(strquery,'content-type','application/json');
% expdata = hobjson.fields.expression_patterns.data;
expdata = hobjson.fields.expressed_in.data;
names = {};
for p=1:numel(expdata)
    % if ~isempty(expdata(p).expressed_in) ...
    %    && isfield(expdata(p).expressed_in,'label')       
    %     names = [names,{expdata(p).expressed_in.label}]; %#ok<AGROW>
    % end
    if ~iscell(expdata)
        if ~isempty(expdata(p).ontology_term) ...
                && isfield(expdata(p).ontology_term,'label')
            names = [names,{expdata(p).ontology_term.label}]; %#ok<AGROW>
        end
    else
        if ~isempty(expdata{p}.ontology_term) ...
                && isfield(expdata{p}.ontology_term,'label')
            names = [names,{expdata{p}.ontology_term.label}]; %#ok<AGROW>
        end
    end
end
names_unique = unique(regexprep(names(:),userdata.strregexp,userdata.strrep));
hmtm.setDataVector(names_unique,userdata.strColumnNamesMark);
end


function setTitle(~,~)
userdata = get(gcf,'UserData');
numTable = numel(userdata.hjtable);
numcolor = size(userdata.mat_clmap,1);    
for p=1:numTable
   set(userdata.htitle(p),...
       'String',userdata.hedit(p*2+7).String,...
       'ForegroundColor',userdata.mat_clmap(mod(userdata.hpopup(p*2+8).Value-1,numcolor)+1,:));
end
end


function [htext,hedit,hslider] = createControls(parent,start,height)
htext = uicontrol(parent,...
    'Style','edit',...
    'Enable','inactive',...
    'BackgroundColor',ones(1,3)*0.8,...    
    'Units','normalized',...
    'Position',[0.05,start+height,0.45,height]);
hedit = uicontrol(parent,...
    'Style','edit',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String','0',...
    'Units','normalized',...
    'Position',[0.50,start+height,0.45,height]);
hslider = uicontrol(parent,...
    'Style','slider',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'Units','normalized',...
    'Position',[0.05,start,0.90,height]);
end


function [htext,hedit,hslider] = createControls2(parent,start,height)
htext = uicontrol(parent,...
    'Style','edit',...
    'Enable','inactive',...
    'BackgroundColor',ones(1,3)*0.8,...   
    'Units','normalized',...
    'Position',[0.05,start+height*2,0.2,height]);
hedit(1) = uicontrol(parent,...
    'Style','edit',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String','0',...
    'Units','normalized',...
    'Position',[0.25,start+height*2,0.35,height]);
hedit(2) = uicontrol(parent,...
    'Style','edit',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String','1',...
    'Units','normalized',...
    'Position',[0.6,start+height*2,0.35,height]);
hslider(1) = uicontrol(parent,...
    'Style','slider',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'Value',0,...
    'Units','normalized',...
    'Position',[0.05,start+height,0.90,height]);
hslider(2) = uicontrol(parent,...
    'Style','slider',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'Value',1,...
    'Units','normalized',...
    'Position',[0.05,start,0.90,height]);
end


function [htext,hpopup,hedit,hbutton,hmtm,hstm,hjtable]...
    = createControls3(hparent,data,columnNames)

htext(1) = uicontrol(hparent,...
    'String','dummy',...
    'Style','edit',...
    'Enable','inactive',...
    'BackgroundColor',ones(1,3)*0.8,...   
    'Units','normalized',...
    'Position',[0.05,0.9,0.35,0.05]);    
hpopup(1) = uicontrol(hparent,...
    'Style','popup',...
    'String','dummy',...
    'Units','normalized',...
    'Position',[0.40,0.9,0.55,0.05]);
hbutton(1) = uicontrol(hparent,...
    'Style','pushbutton',...
    'String','Search',...
    'Units','normalized',...
    'Position',[0.05,0.85,0.35,0.05]);
hedit(1) = uicontrol(hparent,...
    'Style','edit',...    
    'Units','normalized',...
    'Position',[0.40,0.85,0.55,0.05]);

htext(2) = uicontrol(hparent,...
    'String','Color',...
    'Style','edit',...
    'Enable','inactive',...
    'BackgroundColor',ones(1,3)*0.8,...   
    'Units','normalized',...
    'Position',[0.05,0.75,0.35,0.05]);    
hpopup(2) = uicontrol(hparent,...
    'Style','popup',...
    'String','dummy',...
    'Units','normalized',...
    'Position',[0.40,0.75,0.55,0.05]);

hbutton(2) = uicontrol(hparent,...
    'Style','pushbutton',...
    'String','Add Row',...
    'Units','normalized',...
    'Position',[0.05,0.05,0.45,0.05]);

uicontrol(hparent,...
    'Style','pushbutton',...
    'String','Save',...
    'Tag',['button_vpe_',hparent.Title,'_save'],...
    'Units','normalized',...    
    'Position',[0.50,0.05,0.45,0.05]);

%%% make TableModel and override getColumnClass function
%%% myTableModel.class should be placed in current dir or javaclasspath
hmtm = javaObjectEDT(MyTableModelVisPromExp(data,columnNames));

%%% create SortableTable of JIDE from TableModel and enable sorting & filtering
hstm = javaObjectEDT(com.jidesoft.grid.SortableTableModel(hmtm));
hjtable = javaObjectEDT(com.jidesoft.grid.SortableTable(hstm));
tableHeader = javaObjectEDT(com.jidesoft.grid.AutoFilterTableHeader(hjtable));

tableHeader.setAutoFilterEnabled(true);
tableHeader.setShowFilterName(true);
tableHeader.setShowFilterIcon(true);
hjtable.setTableHeader(tableHeader);
hjtable.setAutoResort(false);

%%% In order to enable sliders, disable AutoResize
hjtable.setAutoResizeMode(hjtable.AUTO_RESIZE_OFF); 

%%% If clicked twice, enable editing
hjtable.setClickCountToStart(2);

%%% display
[~,hjcontainer] = javacomponent(javax.swing.JScrollPane(hjtable),[0,0,1,1],hparent);
set(hjcontainer,'Units','normalized','Position',[0.05,0.1,0.9,0.60]);

%%% enabling convenient features
installer = com.jidesoft.grid.TableHeaderPopupMenuInstaller(hjtable);
pmCustomizer1=com.jidesoft.grid.AutoResizePopupMenuCustomizer;
installer.addTableHeaderPopupMenuCustomizer(pmCustomizer1);
pmCustomizer2=com.jidesoft.grid.TableColumnChooserPopupMenuCustomizer;
installer.addTableHeaderPopupMenuCustomizer(pmCustomizer2);

end


