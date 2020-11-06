%%% roiedit3d - annotation 
%%% by Stephen Wu @ 2016/09/21
%%% 
%%% This is the add-on part to the existing roiedit3d for auto-annotation
%%% and human annotation support based on atlas & majority voting idea
%%% 
%%% updated @ 2016/10/25
%%% 1. use hungarian algorithm for any neuron number matchings
%%% 2. fixed "hist" count if only 1 consistent answer found
%%% 3. added CPD, rank-feedback, parallel-run functions
%%% updated @ 2018/04/18
%%% 1. add hungarian algorithm for final best matches
%%% 2. remove CPD and rank-feedback
%%% 3. upgrade capability of labels/promoters for annotations

% dummy function
function roiedit3d_anno_v3_mod7()
%     create_anno_GUI();

% loadAtlas();
% loadMainData();
% loadOption();
% run_annotation();
% update_table();

init();

end


function init()

hbutton = findobj('Tag','button_anno_label_details');
set(hbutton,'Callback',@(~,~)createFigureAnnoDetail);

hbutton = findobj('Tag','button_anno_run');
set(hbutton,'Callback',@(~,~)run_annotation_main);

hbutton = findobj('Tag','button_anno_approve');
set(hbutton,'Callback',@(~,~)run_approve_anno);

hbutton = findobj('Tag','button_anno_approvelist');
set(hbutton,'Callback',@(~,~)createFigureApproveList);

userdataMain = getUserDataMain();
if isfield(userdataMain,'para_A')
    userdata.para_A = userdataMain.para_A;
else
    userdata.para_A.labels_neuron = {};
    userdata.para_A.labels_name = {};
    userdata.para_A.labels_dataname = {};
    userdata.para_A.approve_neuron = {};
    userdata.para_A.approve_name = {};
    userdata.para_A.N = 100;
    userdata.para_A.N_rank = 3;
    userdata.para_A.CPU = 1;
    userdata.para_A.TF_fHuman = true;
    userdata.para_A.neu_human = {};
    userdata.para_A.TF_labels = false;
    userdata.para_A.w_labels = [];
end
setUserData(userdata);

%%%% TOBE MODIFIED %%%%%%
%%%% Annotation page loading refresh (ex check box of elabel/promoterf etc)
approve_name = cat(1,{'Unnamed'},userdata.para_A.approve_name);
set(findobj('Tag','popup_anno_moveneu'),'String',approve_name);

set(findobj('Tag','check_anno_feat_eat4'),'Value',userdata.para_A.TF_labels);
set(findobj('Tag','check_anno_feat_human'),'Value',userdata.para_A.TF_fHuman);
set(findobj('Tag','edit_anno_nA'),'String',userdata.para_A.N);
set(findobj('Tag','edit_anno_MV_rank'),'String',userdata.para_A.N_rank);
set(findobj('Tag','edit_anno_MV_nCPU'),'String',userdata.para_A.CPU);

end




function createFigureApproveList()

if isempty(findobj('Tag','figure_vis_prom_exp'))
    visualize_promoter_expression_13();
end

if isempty(findobj('Tag','figure_approve_list'))
    figure(...
        'Tag','figure_approve_list',...
        'WindowStyle','normal',...
        'Position',[420,300,200,300],...
        'MenuBar','none',...
        'NumberTitle','off',...
        'Name','Approve Auto-Annotation list');
end

createApproveListPanel();

end


function createFigureAnnoDetail()

if isempty(findobj('Tag','figure_vis_prom_exp'))
    visualize_promoter_expression_13();
end

if isempty(findobj('Tag','figure_labelpromoter_detail'))
    figure(...
        'Tag','figure_labelpromoter_detail',...
        'WindowStyle','normal',...
        'Position',[420,100,300,350],...
        'MenuBar','none',...
        'NumberTitle','off',...
        'Name','Label/Promoter list');
end

createAnnoDetailPanel();

end


function run_annotation_main()

loadAtlas();
loadMainData();
loadOption();
run_annotation();
update_table();

end


function run_approve_anno()

userdataMain = getUserDataMain();
hmtm = userdataMain.hmtm;
idxRowFull = 0:hmtm.getRowCount()-1;
flagAlreadyNamed = ~hmtm.regexpfind('^[0-9]+$',idxRowFull,'Name');

userdata = getUserData();
% nameApprove = cat(1,userdata.para_A.approve_neuron{:});
idxSelectApprove = get(findobj('Tag','popup_anno_moveneu'),'Value');
if idxSelectApprove>1
    nameApprove = userdata.para_A.approve_neuron{idxSelectApprove-1};
    flagApprove = hmtm.regexpfind(makeRegExp(nameApprove),idxRowFull,'Name_estim_unique');
else % Unnamed
    flagApprove = ~flagAlreadyNamed;
end
flagApprove = flagApprove & ~hmtm.regexpfind('^(Null|0|N/A)$',idxRowFull,'Name_estim_unique');

flagOverwrite = true;
if any(flagApprove & flagAlreadyNamed)
    disp('Warning! Approved ROIs are already named');
    idxAlreadyExist = find(flagApprove & flagAlreadyNamed);
    uidAlreadyExist = hmtm.getValuesAt(idxAlreadyExist-1,'uniqueID');
    nameAlreadyExist = cell(hmtm.getValuesAt(idxAlreadyExist-1,'Name'));
    disp([{'row','uID','Name'};[num2cell(idxAlreadyExist),num2cell(uidAlreadyExist),nameAlreadyExist]]);
    ret = questdlg({'Approved ROIs are already named','Do you want to overwrite ROI name?'},'Warning');
    if strcmp(ret,'Cancel')
        return;
    elseif strcmp(ret,'NO')
        flagOverwrite = false;
        flagApprove = flagApprove & ~flagAlreadyNamed;
    end
end

nameNew = cell(hmtm.getValuesAt(flagApprove,'Name_estim_unique'));
if numel(unique(nameNew))~=numel(nameNew)
    disp('Error! Overwrapped names are found in name_estim_unique.');
end

flagAlreadyExist = hmtm.regexpfind(makeRegExp(nameNew),idxRowFull,'Name');
if flagOverwrite
    flagAlreadyExist = flagAlreadyExist & ~flagApprove;
end
if any(flagAlreadyExist)
    disp('Error! Approved ROI names already exist.')
    idxAlreadyExist = find(flagAlreadyExist);
    uidAlreadyExist = hmtm.getValuesAt(flagAlreadyExist,'uniqueID');
    nameAlreadyExist = cell(hmtm.getValuesAt(flagAlreadyExist,'Name'));
    disp([{'row','uID','Name'};[num2cell(idxAlreadyExist),num2cell(uidAlreadyExist),nameAlreadyExist]]);
    return;
end

hmtm.setValuesAt(nameNew,flagApprove,'Name',[false,true]);

% hmtm.setValuesAt(true(sum(flagApprove),1),flagApprove,'Check',[true,false]);
hmtm.setValuesAt(true(sum(flagApprove),1),flagApprove,'Unsure',[true,false]);

end


function strRegExp = makeRegExp(cellstr)
cellstr = strcat(upper(cellstr),'|');
strRegExp = cat(2,cellstr{:});
strRegExp = cat(2,'^(',strRegExp(1:end-1),')$');
end



function userdata = getUserData()
userdata = get(findobj('Tag','tab_annotation'),'UserData');
end


function setUserData(userdata)
set(findobj('Tag','tab_annotation'),'UserData',userdata);
end


function userdata_main = getUserDataMain()
userdata_main = get(findobj('Tag','figure_main'),'UserData');
end


function setUserDataMain(userdata_main)
set(findobj('Tag','figure_main'),'UserData',userdata_main);
end


function loadAtlas()
userdata = getUserData();
path_to_atlas = get(findobj('Tag','edit_anno_path'),'String');
if isfield(userdata,'path_to_atlas') ...
    && ~isempty(path_to_atlas) ...
    && strcmp(userdata.path_to_atlas,path_to_atlas)
disp('Custom Msg: Loaded atlas data found: skip loading,')
    return;
end
d = load(path_to_atlas,'atlas');
userdata.atlas = d.atlas;
userdata.path_to_atlas = path_to_atlas;
setUserData(userdata);
disp('Custom Msg: finished loading atlas data.');
end


function loadMainData()

userdata = getUserData();
userdata_main = getUserDataMain();
hmtm = userdata_main.hmtm;

frame_target = userdata_main.himp.getT();
idx_full = (1:hmtm.getRowCount())-1;
frames = hmtm.getValuesAt(idx_full,'frame');
flag_target = frames==frame_target;

% --- update here for getting promoter data form columns
% huic = findobj('Tag','popup_anno_feat_eat4');
% str_col_eat4 = huic.String(huic.Value);

real_data.frame_target = frame_target;
real_data.pi_k = hmtm.getValuesAt(flag_target,'pi_k');
% real_data.pos = hmtm.getValuesAt(flag_target,'mu');
pos = hmtm.getValuesAt(flag_target,'mu');
zscale = 1;
if get(findobj('Tag','check_main_interpz'),'Value')~=1
    zscale = hmtm.zscale; % if interpolation z is ignored, manually interpolate z position
end
pos(:,3) = pos(:,3)*zscale;
real_data.pos = pos;
real_data.sigma = hmtm.getValuesAt(flag_target,'sigma');
real_data.TF_labels = zeros(size(pos,1),length(userdata.para_A.labels_dataname));
% --- update here for getting promoter data form columns
for i = 1:length(userdata.para_A.labels_dataname)
    real_data.TF_labels(:,i) = hmtm.getValuesAt(flag_target,userdata.para_A.labels_dataname{i});
end

names = cell(hmtm.getValuesAt(flag_target,'Name'));
names_anno_human = cell(hmtm.getValuesAt(flag_target,'Name_anno_human'));
tokens = regexp(names_anno_human,'([^, ]+)','tokens');
% real_data.neu_human = cell(size(tokens));
% for p=1:numel(tokens)
%     real_data.neu_human(p) = {cat(2,names(p),tokens{p}{:})};
% end
tokens_mod = cellfun(@(c)cat(2,c{:}),tokens,'UniformOutput',false);
neu_human = cellfun(@cat,num2cell(2*ones(size(names))),names,tokens_mod,...
    'UniformOutput',false);
neu_human = cellfun(@unique,neu_human,'UniformOutput',false); % unique
% real_data.neu_human = neu_human;
real_data.neu_human = cat(2,neu_human,names); % memos and fixed names

disp('Custom Msg: finished loading data.')

userdata.real_data = real_data;
setUserData(userdata);

end


function loadOption()
userdata = getUserData();

userdata.para_A.N =       round(str2double(get(findobj('Tag','edit_anno_nA'),'String')));
% userdata.para_CPD.TF =                     get(findobj('Tag','check_anno_BM_CPD'),'Value');
% userdata.para_CPD.beta =        str2double(get(findobj('Tag','edit_anno_BM_beta'),'String'));
% userdata.para_CPD.lambda =      str2double(get(findobj('Tag','edit_anno_BM_lambda'),'String'));
% userdata.para_CPD.outlier =     str2double(get(findobj('Tag','edit_anno_BM_outlier'),'String'));
userdata.para_A.N_rank = round(str2double(get(findobj('Tag','edit_anno_MV_rank'),'String')));
% userdata.para_MV.Rank_TF =                 get(findobj('Tag','check_anno_MV_rank'),'Value');
% userdata.para_MV.Rank_w =       str2double(get(findobj('Tag','edit_anno_MV_rank_w'),'String'));
% userdata.para_MV.TF_saveW =                get(findobj('Tag','check_anno_MV_saveW'),'Value');
userdata.para_A.CPU =    round(str2double(get(findobj('Tag','edit_anno_MV_nCPU'),'String')));
% userdata.para_fName =                      get(findobj('Tag','edit_anno_feat_fName'),'String');
% userdata.para_A.TF_fSize =                 get(findobj('Tag','check_anno_feat_shape'),'Value');
% userdata.para_A.w_fSize =       str2double(get(findobj('Tag','edit_anno_feat_shape'),'String'));
% userdata.para_A.TF_fInt =                  get(findobj('Tag','check_anno_feat_int'),'Value');
% userdata.para_A.w_fInt =        str2double(get(findobj('Tag','edit_anno_feat_int'),'String'));
%userdata.para_A.TF_fEat4 =                 get(findobj('Tag','check_anno_feat_eat4'),'Value');
%userdata.para_A.w_fEat4 =       str2double(get(findobj('Tag','edit_anno_feat_eat4'),'String'));
userdata.para_A.TF_fHuman =                get(findobj('Tag','check_anno_feat_human'),'Value');
% userdata.para_A.w_fHuman =      str2double(get(findobj('Tag','edit_anno_feat_human'),'String'));
userdata.para_A.neu_human = userdata.real_data.neu_human;

%%% ----- new parts
userdata.para_A.TF_labels =                 get(findobj('Tag','check_anno_feat_eat4'),'Value');
% userdata.para_A.w_labels =       str2double(get(findobj('Tag','edit_anno_feat_eat4'),'String'));
%userdata.para_A.labels_neuron = 0;
%userdata.para_A.labels_name = 0;

setUserData(userdata);

end


function update_table()
userdata = getUserData();
userdata_main = getUserDataMain();
hmtm = userdata_main.hmtm;

frame_target = userdata.real_data.frame_target;
idx_full = (1:hmtm.getRowCount())-1;
flag_target = hmtm.getValuesAt(idx_full,'frame')==frame_target;

name_result = userdata.name_result;
strNameResult = name_result(:,1);
for p=2:size(name_result,2)
    strNameResult = strcat(strNameResult,{', '},name_result(:,p));
end
hmtm.setValuesAt(strNameResult,flag_target,'Name_estim',[false,true]);
hmtm.setValuesAt(userdata.best_result,flag_target,'Name_estim_unique',[true,false]);

end


function addPromoterList()
%%% userdata.para_A.labels_name{k} : label name
%%% userdata.para_A.labels_neuron{k} : neurons list
%%% userdata.para_A.labels_dataname : column name
userdata = getUserData();
[label,list_mark] = loadPromoterAndExpression();
labels_name = {};
labels_neuron = {};
labels_dataname = {};
w_labels = [];
if isfield(userdata,'para_A') && isfield(userdata.para_A,'labels_name')
    labels_name = userdata.para_A.labels_name;
    labels_neuron = userdata.para_A.labels_neuron;
    labels_dataname = userdata.para_A.labels_dataname;
    w_labels = userdata.para_A.w_labels;
end
idxOverlap = find(~cellfun(@isempty,regexp(labels_name,['^',label,'$'])));
if ~isempty(idxOverlap) % update neurons list only;
    userdata.para_A.labels_neuron(idxOverlap) = {list_mark};
else
    userdata.para_A.labels_name = cat(1,labels_name,label);
    userdata.para_A.labels_neuron = cat(1,labels_neuron,{list_mark});
    
    %%% add choice for column of boolean
    userdataMain = getUserDataMain();
    hci = userdataMain.hmtm.getColumnIdentifiers();
    strCol = cell(hci.toArray());
    listColumnNameBoolean = strCol(cellfun(@islogical,userdataMain.defaultdata));    
    userdata.para_A.labels_dataname = cat(1,labels_dataname,listColumnNameBoolean(1));
    userdata.para_A.w_labels = cat(1,w_labels,1);
end
setUserData(userdata);

createAnnoDetailPanel();

end


function deletePromoterList()
%%% userdata.para_A.labels_name{k} : label name
%%% userdata.para_A.labels_neuron{k} : neurons list
%%% userdata.para_A.labels_dataname : column name

hObjs = findobj('-regexp','Tag','checkbox_anno_detail_*');
flag = false(numel(hObjs),1);
for p=1:numel(hObjs)
    hObj = findobj(hObjs,'Tag',['checkbox_anno_detail_',num2str(p)]);
    flag(p) = get(hObj,'Value')==1; %T/F
end
userdata = getUserData();
userdata.para_A.labels_name = userdata.para_A.labels_name(~flag);
userdata.para_A.labels_neuron = userdata.para_A.labels_neuron(~flag);
userdata.para_A.labels_dataname = userdata.para_A.labels_dataname(~flag);
userdata.para_A.w_labels = userdata.para_A.w_labels(~flag);
setUserData(userdata);

createAnnoDetailPanel();

end


function addApproveList()
%%% userdata.para_A.labels_name{k} : label name
%%% userdata.para_A.labels_neuron{k} : neurons list
%%% userdata.para_A.labels_dataname : column name
userdata = getUserData();
[label,list_mark] = loadPromoterAndExpression();
labels_name = {};
labels_neuron = {};
if isfield(userdata,'para_A') && isfield(userdata.para_A,'approve_name')
    labels_name = userdata.para_A.approve_name;
    labels_neuron = userdata.para_A.approve_neuron;
end
idxOverlap = find(~cellfun(@isempty,regexp(labels_name,['^',label,'$'])));
if ~isempty(idxOverlap) % update neurons list only;
    userdata.para_A.approve_neuron(idxOverlap) = {list_mark};
else
    userdata.para_A.approve_name = cat(1,labels_name,label);
    userdata.para_A.approve_neuron = cat(1,labels_neuron,{list_mark});
end
setUserData(userdata);

createApproveListPanel();

end


function deleteApproveList()
%%% userdata.para_A.labels_name{k} : label name
%%% userdata.para_A.labels_neuron{k} : neurons list
%%% userdata.para_A.labels_dataname : column name

hObjs = findobj('-regexp','Tag','checkbox_approve_*');
flag = false(numel(hObjs),1);
for p=1:numel(hObjs)
    hObj = findobj(hObjs,'Tag',['checkbox_approve_',num2str(p)]);
    flag(p) = get(hObj,'Value')==1; %T/F
end
userdata = getUserData();
userdata.para_A.approve_name = userdata.para_A.approve_name(~flag);
userdata.para_A.approve_neuron = userdata.para_A.approve_neuron(~flag);
setUserData(userdata);

createApproveListPanel();

end


function [label,list_mark] = loadPromoterAndExpression()
hfig = findobj('Tag','figure_vis_prom_exp');
if isempty(hfig)
    visualize_promoter_expression_13;
    hfig = findobj('Tag','figure_vis_prom_exp');
end
panelControl = findobj(get(hfig,'Children'),'Title','Control');
tabgroup = get(panelControl,'Children');
tabSelected = tabgroup.SelectedTab;
obj_edit_subset = findobj(get(tabSelected,'Children'),'-regexp','Tag','edit_subset_name*');
label = obj_edit_subset.String;
ret = regexp(get(obj_edit_subset,'Tag'),'edit_subset_name_([0-9]+)','tokens');
idxTab = str2double(ret{1}{1});
userdata_fig = get(hfig,'userdata');
hmtm_base = userdata_fig.hmtm(1);
hmtm_select = userdata_fig.hmtm(idxTab+1);
flag_mark = hmtm_base.regexp(hmtm_select.getDataColumn(0),0)>=0;
list = cell(hmtm_base.getDataColumn(0).toArray());
list_mark = list(flag_mark);
end


function cbColumnName(hObj,~)
ret = regexp(hObj.Tag,'popup_anno_detail_([0-9]*)','tokens');
idx = str2double(ret{1}{1});
userdata = getUserData();
userdata.para_A.labels_dataname(idx) = hObj.String(hObj.Value);
setUserData(userdata);
end


function cbWeight(hObj,~)
ret = regexp(hObj.Tag,'edit_anno_detail_([0-9]*)','tokens');
idx = str2double(ret{1}{1});
userdata = getUserData();
userdata.para_A.w_labels(idx) = str2double(hObj.String);
setUserData(userdata);
end


function createApproveListPanel()

hf = findobj('Tag','figure_approve_list');
clf(hf);

userdata = getUserData();

% n_l = length(userdata.para_A.approve_name);
n_l = 0;
if isfield(userdata,'para_A') && isfield(userdata.para_A,'approve_name')
    n_l = length(userdata.para_A.approve_name);
end

hp_0  = fun_panel(hf,'','panel_approvelist', [0.05,0.05,0.9,0.9]);

tmp_h_button = 0.07;
%tmp_h = 1 - n_l*(tmp_h_text+0.01) - tmp_h_button - 0.1;
tmp_h = 0.95 - tmp_h_button;
fun_pushbutton(hp_0,'Add','button_approve_add',...
    [0.05,tmp_h,0.4,tmp_h_button],'Callback',@(~,~)addApproveList);
fun_pushbutton(hp_0,'Delete','button_approve_delete',...
    [0.5,tmp_h,0.4,tmp_h_button],'Callback',@(~,~)deleteApproveList);

tmp_h_text = 0.07;
tmp_h = tmp_h - tmp_h_button - 0.05;
fun_text(hp_0,'Unnamed','text_approve_unannotated',...
        [0.05,tmp_h,0.4,tmp_h_text]);
tmp_h = tmp_h - tmp_h_text - 0.01;
for i = 1:n_l
    fun_checkbox(hp_0,userdata.para_A.approve_name{i},...
        ['checkbox_approve_',num2str(i)],[0.05,tmp_h,0.9,tmp_h_text]);
%     fun_text(hp_0,userdata.para_A.approve_name{i},['text_approve_',num2str(i)],...
%         [0.05,tmp_h,0.4,tmp_h_text]);
%     fun_popupmenu(hp_0,'(column)',['popup_approve_',num2str(i)],...
%         [0.45,tmp_h,0.45,tmp_h_text]);
    tmp_h = tmp_h - tmp_h_text - 0.01;
end

approve_name = cat(1,{'Unnamed'},userdata.para_A.approve_name);
set(findobj('Tag','popup_anno_moveneu'),'String',approve_name);

end


function createAnnoDetailPanel()
userdata = getUserData();

hf = findobj('Tag','figure_labelpromoter_detail');
clf(hf);

%%% add choice for column of boolean
userdataMain = getUserDataMain();
hci = userdataMain.hmtm.getColumnIdentifiers();
strCol = cell(hci.toArray());
listColumnNameBoolean = strCol(cellfun(@islogical,userdataMain.defaultdata));

n_l = 0;
if isfield(userdata,'para_A') && isfield(userdata.para_A,'labels_name')
    n_l = length(userdata.para_A.labels_name);
end

hp_0  = fun_panel(hf,'','panel_anno_detail', [0.05,0.05,0.9,0.9]);

tmp_h_button = 0.06;
%tmp_h = 1 - n_l*(tmp_h_text+0.01) - tmp_h_button - 0.1;
tmp_h = 0.95 - tmp_h_button;
fun_pushbutton(hp_0,'Add','button_anno_detail_add',...
    [0.05,tmp_h,0.3,tmp_h_button],'Callback',@(~,~)addPromoterList);
fun_pushbutton(hp_0,'Delete','button_anno_detail_delete',...
    [0.4,tmp_h,0.3,tmp_h_button],'Callback',@(~,~)deletePromoterList);

tmp_h_text = 0.06;
tmp_h = tmp_h - tmp_h_button - 0.05;
for i = 1:n_l
    fun_checkbox(hp_0,userdata.para_A.labels_name{i},['checkbox_anno_detail_',num2str(i)],...
        [0.05,tmp_h,0.35,tmp_h_text]);
    fun_popupmenu(hp_0,listColumnNameBoolean,['popup_anno_detail_',num2str(i)],...
        [0.4,tmp_h,0.4,tmp_h_text],'Callback',@cbColumnName,'Value',...
        find(strcmp(listColumnNameBoolean,userdata.para_A.labels_dataname(i))));
    fun_text(hp_0,'w:',['text_anno_detail_',num2str(i)],...
        [0.8,tmp_h,0.05,tmp_h_text]);
    fun_edit(hp_0,'1',['edit_anno_detail_',num2str(i)],...
        [0.85,tmp_h,0.1,tmp_h_text],'Callback',@cbWeight,'String',...
        userdata.para_A.w_labels(i));
    tmp_h = tmp_h - tmp_h_text - 0.01;
end

end



%% Basic functions

% load atlas
function load_atlas_callback(~,~)

%%% button 1 F
%%% open mat file and obtain atlas data
%%% check validity of the atlas data (to be added later)

enableGUI(false);
userdata = get(gcf,'UserData');

[strFileOpen,strPath] = uigetfile('*.mat');
tmp = load(fullfile(strPath,strFileOpen));

userdata.atlas = tmp.atlas;
set(gcf,'UserData',userdata);

disp('Custom Msg: finished loading atlas data.');
disp(' ');
enableGUI(true);

end

% load a single frame data for annotation
function load_frame_callback(~,~)

%%% button 1 F
%%% open mat file and obtain single frame data

enableGUI(false);
userdata = get(gcf,'UserData');

[strFileOpen,strPath] = uigetfile('*.mat');
tmp = load(fullfile(strPath,strFileOpen));

userdata.real_data = tmp.real_data;
set(gcf,'UserData',userdata);

disp('Custom Msg: finished loading one frame data.');
disp(' ');
enableGUI(true);

end

% ----- NEED TO BE UPDATED!
function enableGUI(flag)

userdata = get(gcf,'UserData');

if flag
    set(userdata.huicontrol(:),'Enable','on');
else
    set(userdata.huicontrol(:),'Enable','off');
end

end
% -----

function h = fun_buttongroup(hp,title,tag,pos,varargin)
h = uibuttongroup(...
    'Parent',hp,...
    'Title',title,...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    varargin{:});
end

function h = fun_panel(hp,title,tag,pos,varargin)
h = uipanel(...
    'Parent',hp,...
    'Title',title,...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    varargin{:});
end

function h = fun_edit(hp,string,tag,pos,varargin)
h = uicontrol(...
    'Parent',hp,...
    'Style','edit',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String',string,...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    varargin{:});
end

function h = fun_label(hp,string,tag,pos,varargin)
h = uicontrol(...
    'Parent',hp,...
    'Style','edit',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String',string,...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    'Enable','inactive',...
    'BackgroundColor',ones(1,3)*0.8,...
    varargin{:});
end
 
function h = fun_text(hp,string,tag,pos,varargin)
h = uicontrol(...
    'Parent',hp,...
    'Style','text',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String',string,...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    varargin{:});
end

function h = fun_checkbox(hp,string,tag,pos,varargin)
h = uicontrol(...
    'Parent',hp,...
    'Style','checkbox',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String',string,...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    varargin{:});
end

function h = fun_pushbutton(hp,string,tag,pos,varargin)
h = uicontrol(...
    'Parent',hp,...
    'Style','pushbutton',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String',string,...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    varargin{:});
end

function h = fun_togglebutton(hp,string,tag,pos,varargin)
h = uicontrol(...
    'Parent',hp,...
    'Style','togglebutton',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String',string,...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    varargin{:});
end

function h = fun_radiobutton(hp,string,tag,pos,varargin)
h = uicontrol(...
    'Parent',hp,...
    'Style','radiobutton',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String',string,...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    varargin{:});
end

function h = fun_popupmenu(hp,string,tag,pos,varargin)
h = uicontrol(...
    'Parent',hp,...
    'Style','popupmenu',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String',string,...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    varargin{:});
end

function h = fun_slider(hp,tag,pos,minv,maxv,value,varargin)
h = uicontrol(...
    'Parent',hp,...
    'Style','slider',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    'Min',minv,...
    'Max',maxv,...
    'Value',value,...
    varargin{:});
end

function h = fun_listbox(hp,string,tag,pos,varargin)
h = uicontrol(...
    'Parent',hp,...
    'Style','listbox',...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'String',string,...
    'Tag',tag,...
    'Units','normalized',...
    'Position',pos,...
    varargin{:});
end

%% Annotation codes

% Start running annotation
function run_annotation(~,~)

% enableGUI(false)

userdata = getUserData();

disp('Custom Msg: running automatic annotation...');

% userdata = get(gcf,'UserData');
% 
% %=== get all input parameters
% userdata.para_A.N =       round(str2num(get(findobj(userdata.huicontrol,'Tag','edit_nA'),'String')));
% userdata.para_CPD.TF =                  get(findobj(userdata.huicontrol,'Tag','cBox_CPD'),'Value');
% userdata.para_CPD.beta =        str2num(get(findobj(userdata.huicontrol,'Tag','edit_beta'),'String'));
% userdata.para_CPD.lambda =      str2num(get(findobj(userdata.huicontrol,'Tag','edit_lambda'),'String'));
% userdata.para_CPD.outlier =     str2num(get(findobj(userdata.huicontrol,'Tag','edit_outlier'),'String'));
% userdata.para_A.N_rank = round(str2num(get(findobj(userdata.huicontrol,'Tag','edit_rank'),'String')));
% userdata.para_MV.Rank_TF =              get(findobj(userdata.huicontrol,'Tag','cBox_rank'),'Value');
% userdata.para_MV.Rank_w =       str2num(get(findobj(userdata.huicontrol,'Tag','edit_w_rank'),'String'));
% userdata.para_MV.TF_saveW =              get(findobj(userdata.huicontrol,'Tag','cBox_saveW'),'Value');
% userdata.para_A.CPU =    round(str2num(get(findobj(userdata.huicontrol,'Tag','edit_nCPU'),'String')));
% userdata.para_fName =                   get(findobj(userdata.huicontrol,'Tag','edit_fName'),'String');
% userdata.para_A.TF_fSize =              get(findobj(userdata.huicontrol,'Tag','cBox_fSize'),'Value');
% userdata.para_A.w_fSize =       str2num(get(findobj(userdata.huicontrol,'Tag','edit_w_fSize'),'String'));
% userdata.para_A.TF_fInt =              get(findobj(userdata.huicontrol,'Tag','cBox_fInt'),'Value');
% userdata.para_A.w_fInt =       str2num(get(findobj(userdata.huicontrol,'Tag','edit_w_fInt'),'String'));
% userdata.para_A.TF_fEat4 =              get(findobj(userdata.huicontrol,'Tag','cBox_fEat4'),'Value');
% userdata.para_A.w_fEat4 =       str2num(get(findobj(userdata.huicontrol,'Tag','edit_w_fEat4'),'String'));
% userdata.para_A.TF_fHuman =              get(findobj(userdata.huicontrol,'Tag','cBox_fHuman'),'Value');
% userdata.para_A.w_fHuman =       str2num(get(findobj(userdata.huicontrol,'Tag','edit_w_fHuman'),'String'));

%=== Prepare atlas and target data
tmp_atlas = reshape(userdata.atlas.xyz,size(userdata.atlas.xyz,1),size(userdata.atlas.xyz,2),[]);
N_a = size(tmp_atlas,3);
Nr = userdata.para_A.N;

ind_rnd = randperm(N_a);
ind_r = ind_rnd(1:Nr);
loc_r = tmp_atlas(:,:,ind_r); % ref. xyz
loc_t = userdata.real_data.pos; % target xyz
    % centralize real data
loc_t = loc_t - repmat(mean(loc_t,1),size(loc_t,1),1,1);

Nn = size(loc_t,1);
Nn_r = size(loc_r,1);

%=== Majority voting loop (set up CPU use)
Nh = userdata.para_A.N_rank;
    % prepare for parallel sectioning for jobs
    % --- refine later to avoid crashing when N_CPU high and low Nt
N_minWork = 100; % min. workload for each CPU
if userdata.para_A.CPU > 1
    % prepare CPU in parallel (local)
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        poolobj = parpool(userdata.para_A.CPU);
    elseif poolobj.NumWorkers < userdata.para_A.CPU
        delete(poolobj)
        poolobj = parpool(userdata.para_A.CPU);
    end
    
    % distribute jobs
    N_sec = userdata.para_A.CPU*4;
    if Nr/N_sec < N_minWork
        N_sec = userdata.para_A.CPU*2;
    end
else
    N_sec = 1;
end
    
% %=== rank-feedback or not? (abandoned in this version)
% if userdata.para_MV.Rank_TF
% %     c_std = 0;
%     % setup rank schedule
%     N_rank = [10,7,5,4,3,2];
%     N_rank(N_rank < Nh) = []; % lower rank not reliable after feedback
%     N_rank = [0,N_rank];
% else
%     N_rank = 0;
% end
% w_R = userdata.para_MV.Rank_w; % rank-feedback weight
    
% %=== CPD options (abandoned in this version)
%     % Set the CPD options
% opt_CPD.TF = userdata.para_CPD.TF; % run CPD or not
% opt_CPD.method='nonrigid'; % use nonrigid registration
% opt_CPD.beta=2;            % the width of Gaussian kernel (smoothness)
% opt_CPD.lambda=3;          % regularization weight
% opt_CPD.fgt=0;              % do not use FGT (default)
% opt_CPD.viz=0;          % show every iteration
% opt_CPD.outliers=0.3;   % use 0.6 noise weight to add robustness 
% opt_CPD.normalize=1;    % normalize to unit variance and zero mean before registering (default)
% opt_CPD.scale=1;        % estimate global scaling too (default) [1 when normalize = 1]
% opt_CPD.rot=1;          % estimate strictly rotational matrix (default)
% opt_CPD.corresp=0;      % do not compute the correspondence vector at the end of registration (default)
% opt_CPD.max_it=1000;     % max number of iterations
% opt_CPD.tol=1e-10;       % tolerance
% opt_CPD.TF_msg = false;   % show msg or not

%=== Setup feature contributions
%     % get eat4 weight matrix
% W_f.eat4 = false(Nn,size(loc_r,1));
% if userdata.para_A.TF_fEat4
%     name_eat4 = {'ADAL';'ADAR';'ADLL';'ADLR';'AFDL';'AFDR';'AIBL';'AIBR';'AIML';'AIMR';'AIZL';'AIZR';'AQR';'ASEL';'ASER';'ASGL';'ASGR';'ASHL';'ASHR';'ASKL';'ASKR';'AUAL';'AUAR';'AWCL';'AWCR';'BAGL';'BAGR';'FLPL';'FLPR';'I2L';'I2R';'I5';'IL1DL';'IL1DR';'IL1VL';'IL1VR';'IL1L';'IL1R';'M3L';'M3R';'MI';'OLLL';'OLLR';'OLQDL';'OLQDR';'OLQVL';'OLQVR';'RIAL';'RIAR';'RICL';'RICR';'RIGL';'RIGR';'RIML';'RIMR';'URYDL';'URYDR';'URYVL';'URYVR'};
%     [~,ind_eat4,ind_atl] = intersect(name_eat4,userdata.atlas.names);
%     if length(ind_eat4) ~= length(name_eat4)
%         warning('Custom Msg: Not all eat4 names appear in atlas')
%     end
%     for i = 1:length(ind_atl)
%         W_f.eat4(userdata.real_data.TF_eat4,ind_atl(i)) = true;
%     end
%     W_f.eat4 = -userdata.para_A.w_fEat4*W_f.eat4;
% end

% assumed new variables:
%   userdata.para_A.labels_neuron - cell to store neuron name list
%   userdata.para_A.labels_name - vector of the label's name
%   userdata.para_A.w_labels - vector of weights
%   userdata.para_A.TF_labels - T/F of using labels or not
%   userdata.real_data.TF_labels - matrix to store if neurons are on for a
%       label (binary vector for each label)
    % get diff promotor/labels weight matrix
W_f.labels = false(Nn,Nn_r);
if userdata.para_A.TF_labels
    tmp_ind = find(logical(userdata.para_A.w_labels));
    for k = tmp_ind(:)'
        name_eat4 = userdata.para_A.labels_neuron{k};
        [~,ind_eat4,ind_atl] = intersect(name_eat4,userdata.atlas.names);
        if length(ind_eat4) ~= length(name_eat4)
            warning(['Custom Msg: Not all ',userdata.para_A.labels_name{k},' names appear in atlas'])
        end
        tmp_W_labels = false(Nn,size(loc_r,1));
        for i = ind_atl(:)'
            tmp_W_labels(logical(userdata.real_data.TF_labels(:,k)),i) = true;
        end
        W_f.labels = W_f.labels - userdata.para_A.w_labels(k)*tmp_W_labels;
    end
end
    % get manual update weight matrix
% --- for testing only --- %
% userdata.para_A.neu_human = cell(Nn,2); % col 1 is memo, col 2 is fixed
% ------------------------ %
W_f.human = zeros(Nn,Nn_r);
if userdata.para_A.TF_fHuman
    for i = 1:Nn %%% for memos
%         if ~isempty(userdata.para_A.neu_human{i})
%             [~,~,ind_atl] = intersect(userdata.para_A.neu_human{i},userdata.atlas.names);
%             W_f.human(i,:) = inf;
%             W_f.human(i,ind_atl) = 0;
%         end
        [~,~,ind_atl] = intersect(userdata.para_A.neu_human{i,1},userdata.atlas.names);        
        if ~isempty(ind_atl)
            W_f.human(i,:) = inf;
            W_f.human(i,ind_atl) = 0;
        end

    end
    
    
    for i = 1:Nn % for fixed name
%         if ~isempty(userdata.para_A.neu_human{i})
%             [~,~,ind_atl] = intersect(userdata.para_A.neu_human{i},userdata.atlas.names);
%             W_f.human(i,:) = inf;
%             W_f.human(i,ind_atl) = 0;
%         end
        [~,~,ind_atl] = intersect(userdata.para_A.neu_human{i,2},userdata.atlas.names);        
        if ~isempty(ind_atl)
            if numel(ind_atl) == 1
                W_f.human(:,ind_atl) = inf;
            end 
            W_f.human(i,:) = inf;
            W_f.human(i,ind_atl) = 0;
        end
    end
    
end

W_f = W_f.labels + W_f.human;

%=== Run Bipartite Graph Matching
% if userdata.para_MV.TF_saveW
%     [m_result,m_cost,Run_time] = BGMatch_fast(loc_t,loc_r,W_f,N_rank,w_R,opt_CPD,N_sec);
% else
%     [m_result,m_cost,Run_time] = BGMatch_slow(loc_t,loc_r,W_f,N_rank,w_R,opt_CPD,N_sec);
% end
% if userdata.para_MV.TF_rerun
%     [m_result,m_cost,Run_time] = BGMatch_rerun(userdata.para_MV.W_dist,W_f,N_sec);
% else
%     [m_result,m_cost,Run_time] = BGMatch_basic(loc_t,loc_r,W_f,N_sec);
% end
[m_result,m_cost,Run_time] = BGMatch_basic(loc_t,loc_r,W_f,N_sec);

%=== Collect top Nh rank voting results
vote_result = zeros(Nh,2,Nn);
name_result = cell(Nn,Nh);
% W_final = inf(Nn,Nn_r); % modified by Yu Toyoshima, 2018/5/14, under the instruction by Wu Stephen
W_final = -inf(Nn,Nn_r); 
% W_final = zeros(Nn,Nn_r); 
for j = 1:Nn
    tmp = squeeze(m_result(j,:));
    tmp0 = unique(tmp(:));
    if length(tmp0) == 1
        % modified by Yu Toyoshima, 2018/5/14, under the instruction by Wu Stephen
        %         tmp = zeros(Nh,2);
        %         tmp(1,1) = length(tmp);
        %         tmp(1,2) = tmp0;
        %         tmp(2:3,:) = nan;
        tmp = nan(Nh,2);
        tmp(1,1) = length(tmp);
        tmp(1,2) = tmp0;
    else
        [a,b] = hist(tmp(:),tmp0);
        tmp = flipud(sortrows([a(:),b(:)]));
        if size(tmp,1) < Nh
            tmp(end+1:Nh,:) = nan;
        end
    end
    vote_result(1:Nh,:,j) = fliplr(tmp(1:Nh,:));
    for k = 1:Nh
        if isnan(vote_result(k,1,j))
            name_result{j,k} = 'N/A';
        elseif logical(vote_result(k,1,j))
            name_result{j,k} = userdata.atlas.names{vote_result(k,1,j)};
            % ===== correction by Stephen @ 2018.04.19
            % W_final(j,vote_result(k,1,j)) = k; %use rank info only
            W_final(j,vote_result(k,1,j)) = vote_result(k,2,j); %use #votes
            % ========================
        else
            name_result{j,k} = 'Null';
        end
    end
end

% record results
userdata.vote_result = vote_result;
userdata.name_result = name_result;

% find best 1-to-1 matches

% [m_best,~] = munkres(W_final); % modified by Yu Toyoshima, 2018/5/14, under the instruction by Wu Stephen
% [m_best,~] = munkres(-W_final);
[m_best,~] = munkres(Nr-W_final);
userdata.best_result = cell(Nn,1);
tmp = logical(m_best);
userdata.best_result(tmp) = userdata.atlas.names(m_best(tmp));
userdata.best_result(~tmp) = {'Null'};

% set(gcf,'UserData',userdata);
setUserData(userdata);

disp('Custom Msg: Done! You are Genius!')
disp(['Custom Msg: run time = ~',num2str(Run_time/60,'%.1f'),'min'])

% if userdata.para_A.CPU > 1
%     % auto-shutdown of the parallel pool
%     delete(poolobj)
% end

% enableGUI(true);

end

%%% perform match in parallel based on number of sections
%%%   -> more than N_cpu as buffer when non-uniform job time expected
%%%   -> best to be multiple of N_cpu
%%% Note: count number of total job based on loc_r
%%%
%%% Note2: currently just add feature weight directly while matching!!!
%%%
%%% fast -> save W for faster rank-feedback (careful: storage issue!)
function [m_result,m_cost,tEnd] = BGMatch_basic(loc_t,loc_r,W_f,N_sec)
    
tStart = tic;

% options = optimset('Display','none');

Nr = size(loc_r,3);
Nn_t = size(loc_t,1);
Nn_r = size(loc_r,1);

%  Calc. W & basic matching
if N_sec > 1 % start parallel run
    % setup job distribution
    ind_job = 1:Nr;
        % filled in empty jobs
    ind_job(end+1:ceil(Nr/N_sec)*N_sec) = 0;
        % arrange into "almost" equal sections
    ind_job = reshape(ind_job,N_sec,[]);
        % length of each job
    job_len = sum(logical(ind_job),2);
        % rearrange loc_r into loc_r_job
    loc_r_job = cell(N_sec,1);
    for i = 1:N_sec
        loc_r_job{i} = loc_r(:,:,ind_job(i,1:job_len(i)));
    end
    % Initialization
%     job_W = cell(N_sec,1);
    job_out = cell(N_sec,1);
    job_cost = cell(N_sec,1);
    % parallel job
    parfor iS = 1:N_sec
        tmp_r = loc_r_job{iS};
        Nr0 = job_len(iS);

        % setup weight matrix
        W0 = zeros(Nn_t,Nn_r,Nr0);
        for j = 1:Nr0
            % L-2 norm for all sample to all ref.
            W0(:,:,j) = sqrt(squeeze(sum((repmat(loc_t,[1,1,Nn_r]) - permute(repmat(tmp_r(:,:,j),[1,1,Nn_t]),[3,2,1])).^2,2)));
        end
%         if opt_CPD.TF
%             for j = 1:Nr0
%                 % perform CPD
%                 [Transform, ~] = cpd_register(loc_t,tmp_r(:,:,j),opt_CPD);
%                 % L-2 norm for all sample to all ref.
%                 W0(:,:,j) = sqrt(squeeze(sum((repmat(loc_t,[1,1,Nn_r]) - permute(repmat(Transform.Y,[1,1,Nn_t]),[3,2,1])).^2,2)));
%             end
%         else
%             for j = 1:Nr0
%                 % L-2 norm for all sample to all ref.
%                 W0(:,:,j) = sqrt(squeeze(sum((repmat(loc_t,[1,1,Nn_r]) - permute(repmat(tmp_r(:,:,j),[1,1,Nn_t]),[3,2,1])).^2,2)));
%             end
%         end
%         job_W{iS} = W0;
        
        % 1-to-1 matching
        m_result0 = zeros(Nn_t,Nr0);
        m_cost0 = zeros(Nr0,1);
        for j = 1:Nr0
            [m_result0(:,j),m_cost0(j)] = munkres(W0(:,:,j) + W_f);
        end
        job_out{iS} = m_result0;
        job_cost{iS} = m_cost0;
    end
    % assemble results
%     W = zeros(Nn_t,Nn_r,Nr);
    m_result = zeros(Nn_t,Nr);
    m_cost = zeros(Nr,1);
    for iS = 1:N_sec
%         W(:,:,ind_job(iS,1:job_len(iS))) = job_W{iS};
        m_result(:,ind_job(iS,1:job_len(iS))) = job_out{iS};
        m_cost(ind_job(iS,1:job_len(iS))) = job_cost{iS};
    end
    
%     % rank-feedback (will skip if N_rank = [0])
%     W_r = zeros(Nn_t,Nn_r);
%     for k = 2:length(N_rank)
%         % get ranks info.
%         W0_r = false(Nn_t,Nn_r);
%         tmp_vote = zeros(N_rank(k),2,Nn_t);
%         for j = 1:Nn_t
%             tmp = squeeze(m_result(j,:)); % use previous results
%             tmp0 = unique(tmp(:));
%             if length(tmp0) == 1 % hist doesn't work for this case
%                 tmp = zeros(N_rank(k),2);
%                 tmp(1,1) = length(tmp);
%                 tmp(1,2) = tmp0;
%             else
%                 [a,b] = hist(tmp(:),tmp0);
%                 tmp = flipud(sortrows([a(:),b(:)]));
%                 if size(tmp,1) < N_rank(k)
%                     tmp(end+1:N_rank(k),:) = 0;
%                 end
%             end
%             tmp_vote(:,:,j) = fliplr(tmp(1:N_rank(k),:));
% 
%             tmp = tmp_vote(:,1,j);
%             W0_r(j,tmp(logical(tmp))) = true;
%         end
%         W_r = W_r + W0_r; % accumulate feedback information
% 
%         % feedback matching with parallel run
%         job_out = cell(N_sec,1);
%         job_cost = cell(N_sec,1);
%         parfor iS = 1:N_sec
%             W0 = job_W{iS};
%             Nr0 = job_len(iS);
%             m_result0 = zeros(Nn_t,Nr0);
%             m_cost0 = zeros(Nr0,1);
%             for j = 1:Nr0
%                 tmp = W0(:,:,j); tmp = tmp(:);
%                 tmp = mean(tmp(~isinf(tmp)));
%                 W_tmp = (1-w)*W0(:,:,j)/tmp - w*W_r;
%                 % 1-to-1 matching
%                 [m_result0(:,j),m_cost0(j)] = munkres(W_tmp + W_f);
%             end
%             job_out{iS} = m_result0;
%             job_cost{iS} = m_cost0;
%         end
%             % assemble results
%         m_result = zeros(Nn_t,Nr);
%         m_cost = zeros(Nr,1);
%         for iS = 1:N_sec
%             m_result(:,ind_job(iS,1:job_len(iS))) = job_out{iS};
%             m_cost(ind_job(iS,1:job_len(iS))) = job_cost{iS};
%         end
%     end
    
else
    % setup weight matrix
    W = zeros(Nn_t,Nn_r,Nr);
    for j = 1:Nr
        % L-2 norm for all sample to all ref.
        W(:,:,j) = sqrt(squeeze(sum((repmat(loc_t,[1,1,Nn_r]) - permute(repmat(loc_r(:,:,j),[1,1,Nn_t]),[3,2,1])).^2,2)));
    end
%     if opt_CPD.TF
%         for j = 1:Nr
%             % perform CPD
%             [Transform, ~] = cpd_register(loc_t,loc_r(:,:,j),opt_CPD);
%             % L-2 norm for all sample to all ref.
%             W(:,:,j) = sqrt(squeeze(sum((repmat(loc_t,[1,1,Nn_r]) - permute(repmat(Transform.Y,[1,1,Nn_t]),[3,2,1])).^2,2)));
%         end
%     else
%         for j = 1:Nr
%             % L-2 norm for all sample to all ref.
%             W(:,:,j) = sqrt(squeeze(sum((repmat(loc_t,[1,1,Nn_r]) - permute(repmat(loc_r(:,:,j),[1,1,Nn_t]),[3,2,1])).^2,2)));
%         end
%     end
    
    % 1-to-1 matching
    m_result = zeros(Nn_t,Nr);
    m_cost = zeros(Nr,1);
    for j = 1:Nr
        [m_result(:,j),m_cost(j)] = munkres(W(:,:,j) + W_f);
    end
    
%     % rank-feedback (will skip if N_rank = [0])
%     W_r = zeros(Nn_t,Nn_r);
%     for k = 2:length(N_rank)
%         % get ranks info.
%         W0_r = false(Nn_t,Nn_r);
%         tmp_vote = zeros(N_rank(k),2,Nn_t);
%         for j = 1:Nn_t
%             tmp = squeeze(m_result(j,:)); % use previous results
%             tmp0 = unique(tmp(:));
%             if length(tmp0) == 1 % hist doesn't work for this case
%                 tmp = zeros(N_rank(k),2);
%                 tmp(1,1) = length(tmp);
%                 tmp(1,2) = tmp0;
%             else
%                 [a,b] = hist(tmp(:),tmp0);
%                 tmp = flipud(sortrows([a(:),b(:)]));
%                 if size(tmp,1) < N_rank(k)
%                     tmp(end+1:N_rank(k),:) = 0;
%                 end
%             end
%             tmp_vote(:,:,j) = fliplr(tmp(1:N_rank(k),:));
% 
%             tmp = tmp_vote(:,1,j);
%             W0_r(j,tmp(logical(tmp))) = true;
%         end
%         W_r = W_r + W0_r; % accumulate feedback information
% 
%         % feedback matching
%         for j = 1:Nr
%             tmp = W(:,:,j); tmp = tmp(:);
%             tmp = mean(tmp(~isinf(tmp)));
%             W_tmp = (1-w)*W(:,:,j)/tmp - w*W_r;
% 
%             % 1-to-1 matching
%             [m_result(:,j),m_cost(j)] = munkres(W_tmp + W_f);
%         end
%     end
end

tEnd = toc(tStart);

end

% --- not yet finished (2018.04.18)
function [m_result,m_cost,tEnd] = BGMatch_rerun(W_dist,W_f,N_sec)
    
tStart = tic;

% options = optimset('Display','none');

Nr = size(loc_r,3);
Nn_t = size(loc_t,1);
Nn_r = size(loc_r,1);

%  Calc. W & basic matching
if N_sec > 1 % start parallel run
    % setup job distribution
    ind_job = 1:Nr;
        % filled in empty jobs
    ind_job(end+1:ceil(Nr/N_sec)*N_sec) = 0;
        % arrange into "almost" equal sections
    ind_job = reshape(ind_job,N_sec,[]);
        % length of each job
    job_len = sum(logical(ind_job),2);
        % rearrange loc_r into loc_r_job
    loc_r_job = cell(N_sec,1);
    for i = 1:N_sec
        loc_r_job{i} = loc_r(:,:,ind_job(i,1:job_len(i)));
    end
    % Initialization
    job_W = cell(N_sec,1);
    job_out = cell(N_sec,1);
    job_cost = cell(N_sec,1);
    % parallel job
    parfor iS = 1:N_sec
        tmp_r = loc_r_job{iS};
        Nr0 = job_len(iS);

        % setup weight matrix
        W0 = zeros(Nn_t,Nn_r,Nr0);
% % % % % % % %         if opt_CPD.TF
% % % % % % % %             for j = 1:Nr0
% % % % % % % %                 % perform CPD
% % % % % % % %                 [Transform, ~] = cpd_register(loc_t,tmp_r(:,:,j),opt_CPD);
% % % % % % % %                 % L-2 norm for all sample to all ref.
% % % % % % % %                 W0(:,:,j) = sqrt(squeeze(sum((repmat(loc_t,[1,1,Nn_r]) - permute(repmat(Transform.Y,[1,1,Nn_t]),[3,2,1])).^2,2)));
% % % % % % % %             end
% % % % % % % %         else
            for j = 1:Nr0
                % L-2 norm for all sample to all ref.
                W0(:,:,j) = sqrt(squeeze(sum((repmat(loc_t,[1,1,Nn_r]) - permute(repmat(tmp_r(:,:,j),[1,1,Nn_t]),[3,2,1])).^2,2)));
            end
% % % % % % % %         end
        job_W{iS} = W0;
        
        % 1-to-1 matching
        m_result0 = zeros(Nn_t,Nr0);
        m_cost0 = zeros(Nr0,1);
        for j = 1:Nr0
            [m_result0(:,j),m_cost0(j)] = munkres(W0(:,:,j) + W_f);
        end
        job_out{iS} = m_result0;
        job_cost{iS} = m_cost0;
    end
    % assemble results
%     W = zeros(Nn_t,Nn_r,Nr);
    m_result = zeros(Nn_t,Nr);
    m_cost = zeros(Nr,1);
    for iS = 1:N_sec
%         W(:,:,ind_job(iS,1:job_len(iS))) = job_W{iS};
        m_result(:,ind_job(iS,1:job_len(iS))) = job_out{iS};
        m_cost(ind_job(iS,1:job_len(iS))) = job_cost{iS};
    end
    
%     % rank-feedback (will skip if N_rank = [0])
%     W_r = zeros(Nn_t,Nn_r);
%     for k = 2:length(N_rank)
%         % get ranks info.
%         W0_r = false(Nn_t,Nn_r);
%         tmp_vote = zeros(N_rank(k),2,Nn_t);
%         for j = 1:Nn_t
%             tmp = squeeze(m_result(j,:)); % use previous results
%             tmp0 = unique(tmp(:));
%             if length(tmp0) == 1 % hist doesn't work for this case
%                 tmp = zeros(N_rank(k),2);
%                 tmp(1,1) = length(tmp);
%                 tmp(1,2) = tmp0;
%             else
%                 [a,b] = hist(tmp(:),tmp0);
%                 tmp = flipud(sortrows([a(:),b(:)]));
%                 if size(tmp,1) < N_rank(k)
%                     tmp(end+1:N_rank(k),:) = 0;
%                 end
%             end
%             tmp_vote(:,:,j) = fliplr(tmp(1:N_rank(k),:));
% 
%             tmp = tmp_vote(:,1,j);
%             W0_r(j,tmp(logical(tmp))) = true;
%         end
%         W_r = W_r + W0_r; % accumulate feedback information
% 
%         % feedback matching with parallel run
%         job_out = cell(N_sec,1);
%         job_cost = cell(N_sec,1);
%         parfor iS = 1:N_sec
%             W0 = job_W{iS};
%             Nr0 = job_len(iS);
%             m_result0 = zeros(Nn_t,Nr0);
%             m_cost0 = zeros(Nr0,1);
%             for j = 1:Nr0
%                 tmp = W0(:,:,j); tmp = tmp(:);
%                 tmp = mean(tmp(~isinf(tmp)));
%                 W_tmp = (1-w)*W0(:,:,j)/tmp - w*W_r;
%                 % 1-to-1 matching
%                 [m_result0(:,j),m_cost0(j)] = munkres(W_tmp + W_f);
%             end
%             job_out{iS} = m_result0;
%             job_cost{iS} = m_cost0;
%         end
%             % assemble results
%         m_result = zeros(Nn_t,Nr);
%         m_cost = zeros(Nr,1);
%         for iS = 1:N_sec
%             m_result(:,ind_job(iS,1:job_len(iS))) = job_out{iS};
%             m_cost(ind_job(iS,1:job_len(iS))) = job_cost{iS};
%         end
%     end
    
else
    % setup weight matrix
    W = zeros(Nn_t,Nn_r,Nr);
    if opt_CPD.TF
        for j = 1:Nr
            % perform CPD
            [Transform, ~] = cpd_register(loc_t,loc_r(:,:,j),opt_CPD);
            % L-2 norm for all sample to all ref.
            W(:,:,j) = sqrt(squeeze(sum((repmat(loc_t,[1,1,Nn_r]) - permute(repmat(Transform.Y,[1,1,Nn_t]),[3,2,1])).^2,2)));
        end
    else
        for j = 1:Nr
            % L-2 norm for all sample to all ref.
            W(:,:,j) = sqrt(squeeze(sum((repmat(loc_t,[1,1,Nn_r]) - permute(repmat(loc_r(:,:,j),[1,1,Nn_t]),[3,2,1])).^2,2)));
        end
    end
    
    % 1-to-1 matching
    m_result = zeros(Nn_t,Nr);
    m_cost = zeros(Nr,1);
    for j = 1:Nr
        [m_result(:,j),m_cost(j)] = munkres(W(:,:,j) + W_f);
    end
    
%     % rank-feedback (will skip if N_rank = [0])
%     W_r = zeros(Nn_t,Nn_r);
%     for k = 2:length(N_rank)
%         % get ranks info.
%         W0_r = false(Nn_t,Nn_r);
%         tmp_vote = zeros(N_rank(k),2,Nn_t);
%         for j = 1:Nn_t
%             tmp = squeeze(m_result(j,:)); % use previous results
%             tmp0 = unique(tmp(:));
%             if length(tmp0) == 1 % hist doesn't work for this case
%                 tmp = zeros(N_rank(k),2);
%                 tmp(1,1) = length(tmp);
%                 tmp(1,2) = tmp0;
%             else
%                 [a,b] = hist(tmp(:),tmp0);
%                 tmp = flipud(sortrows([a(:),b(:)]));
%                 if size(tmp,1) < N_rank(k)
%                     tmp(end+1:N_rank(k),:) = 0;
%                 end
%             end
%             tmp_vote(:,:,j) = fliplr(tmp(1:N_rank(k),:));
% 
%             tmp = tmp_vote(:,1,j);
%             W0_r(j,tmp(logical(tmp))) = true;
%         end
%         W_r = W_r + W0_r; % accumulate feedback information
% 
%         % feedback matching
%         for j = 1:Nr
%             tmp = W(:,:,j); tmp = tmp(:);
%             tmp = mean(tmp(~isinf(tmp)));
%             W_tmp = (1-w)*W(:,:,j)/tmp - w*W_r;
% 
%             % 1-to-1 matching
%             [m_result(:,j),m_cost(j)] = munkres(W_tmp + W_f);
%         end
%     end
end

tEnd = toc(tStart);

end

function [assignment,cost] = munkres(costMat)
maxcost = max(costMat(costMat~=inf));
if maxcost<=0
   disp('warn: neg + inf matrix at munkres') 
end
[assignment,cost] = lapjv(costMat);
siz = size(costMat);
if siz(1)>siz(2) % lapjv returns assignment with smaller dimension
    ret = zeros(1,siz(1));
    ret(assignment) = 1:siz(2);
    assignment = ret;
end
end

%%% original munkres function was replaced with lapjv (faster version)
%%% modified by Yu Toyoshima, 2018/4/20
%{ 
% 
% 
% %%% Hungarian algorithm for 1-to-1 matching
% %%% Note: return zero for unassigned elements
% function [assignment,cost] = munkres(costMat)
% % MUNKRES   Munkres (Hungarian) Algorithm for Linear Assignment Problem. 
% %
% % [ASSIGN,COST] = munkres(COSTMAT) returns the optimal column indices,
% % ASSIGN assigned to each row and the minimum COST based on the assignment
% % problem represented by the COSTMAT, where the (i,j)th element represents the cost to assign the jth
% % job to the ith worker.
% %
% % Partial assignment: This code can identify a partial assignment is a full
% % assignment is not feasible. For a partial assignment, there are some
% % zero elements in the returning assignment vector, which indicate
% % un-assigned tasks. The cost returned only contains the cost of partially
% % assigned tasks.
% 
% % This is vectorized implementation of the algorithm. It is the fastest
% % among all Matlab implementations of the algorithm.
% 
% % Examples
% % Example 1: a 5 x 5 example
% %{
% [assignment,cost] = munkres(magic(5));
% disp(assignment); % 3 2 1 5 4
% disp(cost); %15
% %}
% % Example 2: 400 x 400 random data
% %{
% n=400;
% A=rand(n);
% tic
% [a,b]=munkres(A);
% toc                 % about 2 seconds 
% %}
% % Example 3: rectangular assignment with inf costs
% %{
% A=rand(10,7);
% A(A>0.7)=Inf;
% [a,b]=munkres(A);
% %}
% % Example 4: an example of partial assignment
% %{
% A = [1 3 Inf; Inf Inf 5; Inf Inf 0.5]; 
% [a,b]=munkres(A)
% %}
% % a = [1 0 3]
% % b = 1.5
% % Reference:
% % "Munkres' Assignment Algorithm, Modified for Rectangular Matrices", 
% % http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
% 
% % version 2.3 by Yi Cao at Cranfield University on 11th September 2011
% 
% assignment = zeros(1,size(costMat,1));
% cost = 0;
% 
% validMat = costMat == costMat & costMat < Inf;
% bigM = 10^(ceil(log10(sum(costMat(validMat))))+1);
% costMat(~validMat) = bigM;
% 
% % costMat(costMat~=costMat)=Inf;
% % validMat = costMat<Inf;
% validCol = any(validMat,1);
% validRow = any(validMat,2);
% 
% nRows = sum(validRow);
% nCols = sum(validCol);
% n = max(nRows,nCols);
% if ~n
%     return
% end
% 
% maxv=10*max(costMat(validMat));
% 
% dMat = zeros(n) + maxv;
% dMat(1:nRows,1:nCols) = costMat(validRow,validCol);
% 
% %*************************************************
% % Munkres' Assignment Algorithm starts here
% %*************************************************
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   STEP 1: Subtract the row minimum from each row.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minR = min(dMat,[],2);
% minC = min(bsxfun(@minus, dMat, minR));
% 
% %**************************************************************************  
% %   STEP 2: Find a zero of dMat. If there are no starred zeros in its
% %           column or row start the zero. Repeat for each zero
% %**************************************************************************
% zP = dMat == bsxfun(@plus, minC, minR);
% 
% starZ = zeros(n,1);
% while any(zP(:))
%     [r,c]=find(zP,1);
%     starZ(r)=c;
%     zP(r,:)=false;
%     zP(:,c)=false;
% end
% 
% while 1
% %**************************************************************************
% %   STEP 3: Cover each column with a starred zero. If all the columns are
% %           covered then the matching is maximum
% %**************************************************************************
%     if all(starZ>0)
%         break
%     end
%     coverColumn = false(1,n);
%     coverColumn(starZ(starZ>0))=true;
%     coverRow = false(n,1);
%     primeZ = zeros(n,1);
%     [rIdx, cIdx] = find(dMat(~coverRow,~coverColumn)==bsxfun(@plus,minR(~coverRow),minC(~coverColumn)));
%     while 1
%         %**************************************************************************
%         %   STEP 4: Find a noncovered zero and prime it.  If there is no starred
%         %           zero in the row containing this primed zero, Go to Step 5.  
%         %           Otherwise, cover this row and uncover the column containing 
%         %           the starred zero. Continue in this manner until there are no 
%         %           uncovered zeros left. Save the smallest uncovered value and 
%         %           Go to Step 6.
%         %**************************************************************************
%         cR = find(~coverRow);
%         cC = find(~coverColumn);
%         rIdx = cR(rIdx);
%         cIdx = cC(cIdx);
%         Step = 6;
%         while ~isempty(cIdx)
%             uZr = rIdx(1);
%             uZc = cIdx(1);
%             primeZ(uZr) = uZc;
%             stz = starZ(uZr);
%             if ~stz
%                 Step = 5;
%                 break;
%             end
%             coverRow(uZr) = true;
%             coverColumn(stz) = false;
%             z = rIdx==uZr;
%             rIdx(z) = [];
%             cIdx(z) = [];
%             cR = find(~coverRow);
%             z = dMat(~coverRow,stz) == minR(~coverRow) + minC(stz);
%             rIdx = [rIdx(:);cR(z)];
%             cIdx = [cIdx(:);stz(ones(sum(z),1))];
%         end
%         if Step == 6
%             % *************************************************************************
%             % STEP 6: Add the minimum uncovered value to every element of each covered
%             %         row, and subtract it from every element of each uncovered column.
%             %         Return to Step 4 without altering any stars, primes, or covered lines.
%             %**************************************************************************
%             [minval,rIdx,cIdx]=outerplus(dMat(~coverRow,~coverColumn),minR(~coverRow),minC(~coverColumn));            
%             minC(~coverColumn) = minC(~coverColumn) + minval;
%             minR(coverRow) = minR(coverRow) - minval;
%         else
%             break
%         end
%     end
%     %**************************************************************************
%     % STEP 5:
%     %  Construct a series of alternating primed and starred zeros as
%     %  follows:
%     %  Let Z0 represent the uncovered primed zero found in Step 4.
%     %  Let Z1 denote the starred zero in the column of Z0 (if any).
%     %  Let Z2 denote the primed zero in the row of Z1 (there will always
%     %  be one).  Continue until the series terminates at a primed zero
%     %  that has no starred zero in its column.  Unstar each starred
%     %  zero of the series, star each primed zero of the series, erase
%     %  all primes and uncover every line in the matrix.  Return to Step 3.
%     %**************************************************************************
%     rowZ1 = find(starZ==uZc);
%     starZ(uZr)=uZc;
%     while rowZ1>0
%         starZ(rowZ1)=0;
%         uZc = primeZ(rowZ1);
%         uZr = rowZ1;
%         rowZ1 = find(starZ==uZc);
%         starZ(uZr)=uZc;
%     end
% end
% 
% % Cost of assignment
% rowIdx = find(validRow);
% colIdx = find(validCol);
% starZ = starZ(1:nRows);
% vIdx = starZ <= nCols;
% assignment(rowIdx(vIdx)) = colIdx(starZ(vIdx));
% pass = assignment(assignment>0);
% pass(~diag(validMat(assignment>0,pass))) = 0;
% assignment(assignment>0) = pass;
% cost = trace(costMat(assignment>0,assignment(assignment>0)));
% 
% end
% 
% 
% %%% Function for munkres
% function [minval,rIdx,cIdx]=outerplus(M,x,y)
% 
% ny=size(M,2);
% minval=inf;
% for c=1:ny
%     M(:,c)=M(:,c)-(x+y(c));
%     minval = min(minval,min(M(:,c)));
% end
% [rIdx,cIdx]=find(M==minval);
% 
% end
%}