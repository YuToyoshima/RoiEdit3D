function roiedit3d_25_16(varargin)
%%% roiedit3d is a GUI tool for editing rois in 3D images.
%%%
%%% This script use Fiji (an ImageJ distribution).
%%% This script use guava (an java library).
%%%     To avoid conflicting, guava should be loaded before matlab's native
%%%     jars. In order to do, this script should have javaclasstxt.txt that
%%%     includes following:
%%%         <before>
%%%         $matlabroot\3rd-party\guava-18.0.jar
%%% This script use reference.mat as reference position for rotation.
%%%     This file should be in the distributed package, but can be
%%%     overwritten by putting a file of the same name in local folder
%%%     (the folder contains the data file (*.mat)).
%%% This script uses large amount of memory (java heap space) and java.opts
%%% should be located in %matlabroot%\bin\$arch or %mcrroot%\bin\$arch.
%%%
%%% Please compile by mcc as following;
%%% mcc('-m','roiedit3d_25_16.m','-a',fileparts(fileparts(which('Miji'))),...
%%% '-a','atlas_20181230_White=1_MaskPha=1_RemoveHYP=1.mat','-a','reference.mat','-a','guava-18.0.jar',...
%%% '-a','java.opts','-a','javaclasspath.txt','-a','peak_detection_14_param_for_*.txt',...
%%% '-a','merge_emtrace_white.mat','-a','visualize_promoter_expression.mat',...
%%% '-a','geneIDs.csv');
%%%

%{
%%%%%% TODO: intensity correction for depth? %%%%%%
%%%%%% TODO: conserve magnification at interpZ %%%%%%
%%%%%% BUG: Cannot set newly loaded image as composite after using
%%%%%%      arbitrarily rotation function. %%%%%%
%%%%%%
%%%%%% TODO: 回転を許容しないトラッキング? -> 角度と軸長を分離して推定すべき
%%%%%% TODO: treetrack、calc_incoherent, check_track, 輝度表示、範囲制約、などの統合
%%%%%% BUG: ROIEDITの旧バージョンで開けない
%%%%%% TODO: 残差のグラフィカルな表示?
%%%%%% BUG: 輝度のLog表示が勝手にリセットされてしまう
%%%%%% TODO: garbage collectionを強制的に入れないとメモリが不足する? 画像キャッシュのクリアでもよいかも。
%%%%%% TODO: 画像読み込みとフィルタの高速化はかれないか。現状では3sec/95枚と遅い…。
%%%%%% TODO: データ量が多いと、ROIのフィルタリング時に待ち時間が長い。。
%%%%%% TODO: tableのUndo、Redo時に生データをフィルタせずそのまま戻す
%%%%%% TODO: 回転などRoiEdit3d自体のUndo実装
%%%%%% TODO: copyanoの実装
%%%%%%
%%%%%% DONE: integration of detection functionality %%%%%%
%%%%%% DONE: integration of roi display functionality to roitable %%%%%%
%%%%%% DONE: dcvやsifの直接読み込み
%%%%%% DONE: blurred imageをつくるための情報をmatファイルに読み書きできるように。
%%%%%% DONE: blurred imageを生成しなくて済むように改変、
%%%%%%
%%%%%% DONE: integrated with Visualize_Promoter_Expression functionality
%%%%%% DONE: introducing ROI maximum size settings & constrained optimization
%%%%%% DONE: Change values of multiple rows in Customized ROI Manager at once
%%%%%% DONE: Change column titles of Customized ROI Manager
%%%%%% DONE: Save and load column titles and width of Customized ROI Manager
%%%%%% DONE: improve annotation functionality with Stephen Wu
%%%%%%
%%%%%% 25_3:
%%%%%% DONE: correct error when the header of hidden column was required
%%%%%%
%%%%%% 25_4:
%%%%%% DONE: fix 2 bugs in the automatic annotation method.
%%%%%%
%%%%%% 25_5:
%%%%%% DONE: add automatic save and load function for automated processing
%%%%%%
%%%%%% 25_6:
%%%%%% DONE: update atlas for annotation
%%%%%%       (atlas_ryo_60-30-30_white_20180515_roiedit.mat)
%%%%%% DONE: fix a bug on approving estimated annotation for the unnamed
%%%%%%       neurons (roiedit3d_ano_v3_mod4.m)
%%%%%%
%%%%%% 25_7:
%%%%%% DONE: fix a bug on saving the result (enbag at 25_5)
%%%%%%
%%%%%% 25_8:
%%%%%% DONE: fix a bug on reading dcv file (for processing unsigned short)
%%%%%% DONE: enable to read splitted tif file (as mif)
%%%%%% DONE: auto update of annotation results
%%%%%% DONE: fix a bug on saving max projected ROIs (for single-channel image)
%%%%%% DONE: enable to scale max projected ROIs on saving
%%%%%%
%%%%%% 25_9:
%%%%%% DONE: fix a bug on detection (interpZ handles correctly after
%%%%%%       peak_detection_14_3)
%%%%%% DONE: optimize parameters of detection for IshiharaLab4DAno
%%%%%%
%%%%%% 25_10:
%%%%%% DONE: update atlas for annotation
%%%%%%       (atlas_ryo_60-30-30_white_20180915.mat)
%%%%%% DONE: update mif file handling
%%%%%%
%%%%%% 25_12:
%%%%%% DONE: enable changing size of image cache
%%%%%%
%%%%%% 25_13:
%%%%%% DONE: fix a bug on right click menu for table without selecting rows
%%%%%% DONE: fix a bug on right click menu when single row selected
%%%%%% DONE: fix a bug on search function for visualize_promoter_expression
%%%%%%       (in visualize_promoter_expression_11.m)
%%%%%% DONE: update atlas for annotation
%%%%%%       (atlas_20181230_White=1_MaskPha=1_RemoveHYP=1.mat)
%%%%%%
%%%%%% 25_14:
%%%%%% DONE: add a functionality for making max projected movie with roi
%%%%%%
%%%%%% 25_15:
%%%%%% DONE: fix a bug in the function for making max projected movie & roi 
%%%%%%
%%%%%% 25_16:
%%%%%% DONE: in visualize_promoter_expression_13, add editor function and 
%%%%%%       median positions of neurons in the Iinolab neuron ID dataset
%%%%%%

%}

%%% copy guava-18.0.jar and java.opts to correct path
staticpath = javaclasspath('-static');
staticpath_guava = staticpath{contains(staticpath,'guava-18.0.jar')};
path_javaopts = fullfile(matlabroot,'bin',computer('arch'));

flag_copy = false;

%%% for debug
disp(staticpath);
disp(staticpath_guava);
disp(path_javaopts);

try
    if ~exist(staticpath_guava,'file')
        flag_copy = true;
        disp(['Copy ',which('guava-18.0.jar')]);
        disp([' to ',staticpath_guava]);
        mkdir(fileparts(staticpath_guava));
        copyfile(which('guava-18.0.jar'),staticpath_guava);
        disp('guava-18.0.jar was copied successfully.');
    end
    
    if ~exist(fullfile(path_javaopts,'java.opts'),'file')
        flag_copy = true;
        disp(['Copy ',which('java.opts')]);
        disp([' to ',fullfile(path_javaopts,'java.opts')]);
        copyfile(which('java.opts'),path_javaopts);
        disp('java.opts was copied successfully.');
    end
    
catch err
    flag_copy = true;
    disp('Copy failed. Please retry as an administrator.');
    disp(getReport(err,'extended'));
end

if flag_copy
    % Component was copied in this session. In order to use the component,
    % it is necessary to restart the program.
    disp('Please restart this application.');
    pause(10000);
    return;
end


%%% run MIJ
if ~exist('MIJ','class') || numel(ij.IJ.getInstance())==0
    strDir = pwd();
    Miji();
    cd(strDir);
end

%%% create GUI
create_GUI();

disp(' ');
disp('Welcome to RoiEdit3D.');
disp(' ');

if nargin>0
    data_load_callback([],[],varargin{1});
end

end


%%% create GUI

function create_GUI()

hf = figure(...
    'Tag','figure_main',...
    'WindowStyle','normal',...
    'Position',[100,100,320,450],...
    'MenuBar','none',...
    'NumberTitle','off',...
    'Name','RoiEdit3D');

%%% tabgroup
htabgroup = uitabgroup('Parent', hf,'TabLocation','left');
htab_main     = uitab(htabgroup,'Tag','tab_main',      'Title','Main');
htab_detect   = uitab(htabgroup,'Tag','tab_detect',    'Title','Detect');
htab_track    = uitab(htabgroup,'Tag','tab_track',     'Title','Track');
htab_trackqp  = uitab(htabgroup,'Tag','tab_trackqp',   'Title','TrackQP');
htab_annotate = uitab(htabgroup,'Tag','tab_annotation','Title','Annotate');
htab_analyze  = uitab(htabgroup,'Tag','tab_analyze',   'Title','Analyze');
htab_vistrack = uitab(htabgroup,'Tag','tab_vistrack',  'Title','VisTrack');
htab_movie    = uitab(htabgroup,'Tag','tab_movie',     'Title','Movie');
htab_option   = uitab(htabgroup,'Tag','tab_option',    'Title','Option');

create_GUI_main(htab_main); % main panel
create_GUI_detect(htab_detect); % detect panel
create_GUI_track(htab_track); % track panel
create_GUI_trackqp(htab_trackqp); % trackqp panel
create_GUI_annotate(htab_annotate); % annotation tab
create_GUI_analyze(htab_analyze); % analyze tab
create_GUI_vistrack(htab_vistrack); % vistrack tab
create_GUI_movie(htab_movie); % movie tab
create_GUI_option(htab_option); % option tab

loadDetectPreset(); % set preset values for detection to UI

enableGUI(false);

end


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


function create_GUI_main(hp)

%%% create panels
hp_main    = fun_panel(hp,'Main',    'panel_main',    [0.05,0.78,0.9,0.20]);
hp_roidisp = fun_panel(hp,'Show ROI & Name',     'panel_roi',     [0.05,0.51,0.9,0.25]);
hp_imdisp  = fun_panel(hp,'Image',   'panel_image',   [0.05,0.27,0.9,0.22]);
% hp_filter = fun_panel(hp,'Filtering',  'panel_condition',  [0.05,0.15,0.9,0.22]);
hp_rotate  = fun_panel(hp,'Rotation','panel_rotation',[0.05,0.02,0.9,0.23]);


%%% main section
fun_pushbutton(hp_main,'Load','load',[0.05,0.67,0.4,0.3],'Callback',@data_load_callback);
fun_pushbutton(hp_main,'Save','save',[0.05,0.35,0.4,0.3],'Callback',@data_save_callback);

fun_pushbutton(hp_main,'Undo','button_main_undo',[0.55,0.67,0.4,0.3],...
    'Callback',@(~,~)undo);
fun_pushbutton(hp_main,'Redo','button_main_redo',[0.55,0.35,0.4,0.3],...
    'Callback',@(~,~)redo);

fun_checkbox(hp_main,'Auto save',    'check_main_autosave',     [0.05,0.02,0.35,0.3],'Value',true);
fun_checkbox(hp_main,'Compatibility','check_main_compatibility',[0.55,0.02,0.4,0.3],'Value',true);


%%% roidisplay section
strcolor = {'red','green','blue','cyan','magenta','yellow'};

fun_checkbox(hp_roidisp,'ROI', 'check_main_showroi', [0.05,0.75,0.3,0.2],...
    'Value',true,'Callback',@checkbox_showRoi_callback);
fun_checkbox(hp_roidisp,'Name','check_main_showname',[0.05,0.55,0.3,0.2],...
    'Value',true,'Callback',@checkbox_showName_callback);
fun_checkbox(hp_roidisp,'Estimated Name', 'check_main_showname_estimated', [0.45,0.75,0.50,0.2],...
    'Value',false,'Callback',@checkbox_showName_callback);
fun_checkbox(hp_roidisp,'Fixed Name','check_main_showname_fixed',[0.45,0.55,0.50,0.2],...
    'Value',true,'Callback',@checkbox_showName_callback);


fun_label(hp_roidisp,   'Stroke Color','label_main_strokecolor',[0.05,0.3,0.45,0.25]);
fun_popupmenu(hp_roidisp,strcolor,     'popup_main_strokecolor',[0.50,0.3,0.45,0.25],...
    'Value',1,'Callback',@popupmenu_roicolor_callback);
fun_label(hp_roidisp,   'Selection Color','label_main_selectioncolor',[0.05,0.05,0.45,0.25]);
fun_popupmenu(hp_roidisp,strcolor,        'popup_main_selectioncolor',[0.50,0.05,0.45,0.25],...
    'Value',4,'Callback',@popupmenu_roicolor_callback);



%%% image display section
fun_checkbox(hp_imdisp,'Mirror Image',       'check_main_rotate_mirror',[0.05,0.76,0.45,0.22],...
    'Callback',{@rotate_callback,'mirror'});
fun_checkbox(hp_imdisp,'Interpolation Z',    'check_main_interpz',      [0.05,0.54,0.45,0.22],...
    'Callback',@checkbox_interpz_callback);
fun_checkbox(hp_imdisp,'Subtract BG',        'check_main_subtract',     [0.05,0.32,0.45,0.22],...
    'Callback',@checkbox_subtract_callback,'TooltipString','Subtract Background');
fun_checkbox(hp_imdisp,'Alignment',          'check_main_align',        [0.50,0.76,0.40,0.22],...
    'Callback',@checkbox_align_callback);
fun_checkbox(hp_imdisp,'Log intensity',      'check_main_logarithm',    [0.50,0.54,0.45,0.22],...
    'Callback',@checkbox_logarithm_callback);
fun_checkbox(hp_imdisp,'MaxProjZ',           'check_main_proj',         [0.50,0.32,0.40,0.22],...
    'Callback',@checkbox_proj_callback);

fun_label(hp_imdisp,'Subtract Background Radius','label_main_subtract',[0.05,0.02,0.75,0.30]);
fun_edit(hp_imdisp, '50',                        'edit_main_subtract', [0.80,0.02,0.15,0.30],...
    'Callback',@edit_subtract_radius_callback,'TooltipString','Radius of rolling-ball for background subtraction');



%%% rotation section
fun_togglebutton(hp_rotate,'X',  'button_main_rotate_x',  [0.1,0.65,0.1,0.30],...
    'Callback',{@rotate_callback,'x'});
fun_togglebutton(hp_rotate,'Y',  'button_main_rotate_y',  [0.3,0.65,0.1,0.30],...
    'Callback',{@rotate_callback,'y'});
fun_togglebutton(hp_rotate,'Z',  'button_main_rotate_z',  [0.5,0.65,0.1,0.30],...
    'Callback',{@rotate_callback,'z'});
% fun_pushbutton(hp_rotate,'Dye','button_main_rotate_dye',[0.7,0.65,0.2,0.30],...
%     'Callback',{@rotate_callback,'dye'}); % not work yet
fun_togglebutton(hp_rotate,'Dye','button_main_rotate_dye',[0.7,0.65,0.2,0.30],...
    'HandleVisibility','off','Enable','off'); % deprecated
fun_text(hp_rotate,'Arbitrary Rotation Along X axis','text_main_rotate_arbit',[0.1,0.35,0.8,0.2]);
fun_edit(hp_rotate,'0','edit_main_rotate_arbit',[0.1,0.05,0.2,0.30],...
    'Callback',{@rotate_callback,'arbit_edit'});
fun_slider(hp_rotate,'slider_main_rotate_arbit',[0.3,0.05,0.6,0.30],...
    0,360,0,'Callback',{@rotate_callback,'arbit_slider'});

end


function create_GUI_detect(hp)

%%% create panels
hp_main   = fun_panel(hp,'Main',     'panel_detect_main',  [0.05,0.85,0.9,0.13]);
hp_align  = fun_panel(hp,'Alignment','panel_detect_align', [0.05,0.61,0.9,0.22]);
hp_detect = fun_panel(hp,'Detection','panel_detect_detect',[0.05,0.19,0.9,0.40]);
hp_fit    = fun_panel(hp,'Fitting'  ,'panel_detect_fit',   [0.05,0.02,0.9,0.15]);

%%% main section
fun_pushbutton(hp_main,'Load Image','button_detect_main_load',[0.05,0.50,0.40,0.48],...
    'Callback',@load_mif_callback);
% fun_pushbutton(hp_main,'Copy ROI info','button_detect_main_detect', [0.50,0.67,0.45,0.3],...
%     'Callback',@(x,y)copyROIInfo());
fun_pushbutton(hp_main,'Detect','button_detect_main_detect', [0.50,0.50,0.45,0.48],...
    'Callback',@(x,y)detectROI());
fun_label(hp_main,'Choose preset','edit_detect_main_preset', [0.05,0.02,0.45,0.48]);
fun_popupmenu(hp_main,'(preset)','popup_detect_main_preset', [0.50,0.02,0.45,0.48],...
    'Callback',@(x,y)loadDetectPreset(),'CreateFcn',@(x,y)initDetectPreset());
% fun_checkbox(hp_main,'debug mode','check_detect_main_debug', [0.05,0.02,0.45,0.3]);
% fun_pushbutton(hp_main,'Detect','button_detect_main_detect', [0.55,0.02,0.40,0.3],...
%     'Callback',@(x,y)detectROI(),'visible','off'); % does not work correctly


%%% alignment option
fun_label(hp_align,'Channel',        'label_detect_align_channel', [0.05,0.78,0.45,0.19]);
fun_edit( hp_align,'1',              'edit_detect_align_channel',  [0.50,0.78,0.45,0.19],...
    'TooltipString','Specify channel for ROI detection');
fun_label(hp_align,'Frame',          'label_detect_align_frame',   [0.05,0.59,0.45,0.19]);
fun_edit( hp_align,'1',              'edit_detect_align_frame',    [0.50,0.59,0.45,0.19],...
    'TooltipString','Specify frame for ROI detection');
fun_label(hp_align,'Filter radius',  'label_detect_align_filter',  [0.05,0.40,0.45,0.19]);
fun_edit(hp_align,'2',               'edit_detect_align_filter',   [0.50,0.40,0.45,0.19],...
    'TooltipString','Specify median and gaussian filter radius for alignment; If set to 0, the filtering step will be ignored.');
fun_label(hp_align,'Max move',       'label_detect_align_max',     [0.05,0.21,0.45,0.19]);
fun_edit(hp_align,'[10, 10]',        'edit_detect_align_max',      [0.50,0.21,0.45,0.19],...
    'TooltipString','[x,y] pixels of maximum difference between adjuscent z-slices; Detected movement larger than this value will be ignored.');
fun_checkbox(hp_align,'Use centroid','check_detect_align_centroid',[0.05,0.02,0.45,0.19],...
    'TooltipString','Whether to use centroid method (faster) than cross-corelation method (precise) in subpixel alignment (only for t-frame)');
fun_checkbox(hp_align,'Subtract BG', 'check_detect_align_subtract',[0.50,0.02,0.45,0.19],...
    'TooltipString','Whether to subtract median value from image in subpixel alignment (sensitive to weak signals)');

%%% detection option
fun_label(hp_detect,'Prefilter: Method', 'label_detect_detect_prefilter_method',  [0.05,0.85,0.45,0.10]);
fun_edit(hp_detect,'Median...',          'edit_detect_detect_prefilter_method',   [0.50,0.85,0.45,0.10],...
    'TooltipString','Method for denoising');
fun_label(hp_detect,'Prefilter: Option', 'label_detect_detect_prefilter_option',  [0.05,0.75,0.45,0.10]);
fun_edit(hp_detect,'radius=2 stack',     'edit_detect_detect_prefilter_option',   [0.50,0.75,0.45,0.10],...
    'TooltipString','Option for denoising');
fun_label(hp_detect,'Blur: Method',      'label_detect_detect_blur_method',       [0.05,0.65,0.45,0.10]);
fun_edit(hp_detect,'Gaussian Blur...',   'edit_detect_detect_blur_method',        [0.50,0.65,0.45,0.10],...
    'TooltipString','Method for blurring');
fun_label(hp_detect,'Blur: Option',      'label_detect_detect_blur_option',       [0.05,0.55,0.45,0.10]);
fun_edit(hp_detect,'sigma=2 stack',      'edit_detect_detect_blur_option',        [0.50,0.55,0.45,0.10],...
    'TooltipString','Option for blurring');
fun_label(hp_detect,'Threshold',         'label_detect_detect_threshold',         [0.05,0.45,0.45,0.10]);
fun_edit(hp_detect,'Triangle',           'edit_detect_detect_threshold',          [0.50,0.45,0.45,0.10],...
    'TooltipString','Method for thresholding; For CCD/CMOS images, try Huang, Li, Triangle, etc.');
fun_label(hp_detect,'Local Maxima',      'label_detect_detect_local',             [0.05,0.35,0.45,0.10]);
fun_edit(hp_detect,'[8,8,8]',            'edit_detect_detect_local',              [0.50,0.35,0.45,0.10],...
    'TooltipString','[x,y,z] (pixels) of dilation radius for searching local maxima');
fun_label(hp_detect,'Min. voxel',        'label_detect_detect_min',               [0.05,0.25,0.45,0.10]);
fun_edit(hp_detect,'1000',               'edit_detect_detect_min',                [0.50,0.25,0.45,0.10],...
    'TooltipString','Minimum voxel size for removing too-small objects');
fun_label(hp_detect,'Curvature: distance','label_detect_detect_curvature_distance',[0.05,0.15,0.45,0.10]);
fun_edit(hp_detect,'[7,7,7]',             'edit_detect_detect_curvature_distance', [0.50,0.15,0.45,0.10],...
    'TooltipString','A negative curvature voxel is ignored when the distance from segmentation border is smaller than this [x,y,z] value');
fun_label(hp_detect,'Curvature: #voxels','label_detect_detect_curvature_#voxels', [0.05,0.05,0.45,0.10]);
fun_edit(hp_detect,'100',                'edit_detect_detect_curvature_#voxels',  [0.50,0.05,0.45,0.10],...
    'TooltipString','A segmented area is marked as under-segmented when the number of negative curvature voxels is larger than this value');

%%% fitting option
fun_label(hp_fit,'RemoveFP: numloop',  'label_detect_fit_removefp_numloop',  [0.05,0.66,0.55,0.32]);
fun_edit(hp_fit,'1',                   'edit_detect_fit_removefp_numloop',   [0.60,0.66,0.35,0.32],...
    'TooltipString','Maximum number of loops for remove false positives.');
fun_label(hp_fit,'RemoveFP: mindist',  'label_detect_fit_removefp_mindist',  [0.05,0.34,0.55,0.32]);
fun_edit(hp_fit,'1.5',                 'edit_detect_fit_removefp_mindist',   [0.60,0.34,0.35,0.32],...
    'TooltipString','[um]; If minimum distance is smaller than this threshold, either roi will be removed');
fun_label(hp_fit,'RemoveFP: minscdist','label_detect_fit_removefp_minscdist',[0.05,0.02,0.55,0.32]);
fun_edit(hp_fit,'1',                   'edit_detect_fit_removefp_minscdist', [0.60,0.02,0.35,0.32],...
    'TooltipString','If minimum relative distance is smaller than this threshold, either roi will be removed');

end


function create_GUI_track(hp)
%%% create panels
hp_source = fun_buttongroup(hp,'Source frame','buttongroup_track_source',[0.05,0.85,0.9,0.13]);
hp_target = fun_buttongroup(hp,'Target frame','buttongroup_track_target',[0.05,0.70,0.9,0.13]);
hp_option = fun_panel(hp,'Option',            'panel_track_option',      [0.05,0.02,0.9,0.66]);


%%% target section
fun_radiobutton(hp_source,'Default (Prev|Curr)','button_track_source_default',[0.05,0.50,0.60,0.45]);
fun_radiobutton(hp_source,'Similar',            'button_track_source_similar',[0.65,0.50,0.30,0.45],...
    'visible','off'); %%% currently not working
fun_radiobutton(hp_source,'Specify',            'button_track_source_specify',[0.05,0.05,0.45,0.45]);
fun_edit(hp_source,'',                          'edit_track_source_specify',  [0.35,0.05,0.60,0.45]);


%%% source section
fun_radiobutton(hp_target,'<<',     'button_track_target_init',    [0.05,0.50,0.15,0.45]);
fun_radiobutton(hp_target,'<',      'button_track_target_previous',[0.20,0.50,0.15,0.45]);
fun_radiobutton(hp_target,'Current','button_track_target_current', [0.35,0.50,0.30,0.45]);
fun_radiobutton(hp_target,'>',      'button_track_target_next',    [0.65,0.50,0.15,0.45]);
fun_radiobutton(hp_target,'>>',     'button_track_target_end',     [0.80,0.50,0.15,0.45]);
fun_radiobutton(hp_target,'Specify','button_track_target_specify', [0.05,0.05,0.40,0.45]);
fun_edit(hp_target,'',              'edit_track_target_specify',   [0.35,0.05,0.60,0.45]);

set(hp_target,'SelectedObject',findobj('Tag','button_track_target_current'));

%%% option section
fun_label(hp_option,'Threshold: Intensity','label_track_option_thrint',[0.05,0.90,0.6,0.08]);
fun_edit( hp_option,'0',                   'edit_track_option_thrint', [0.65,0.90,0.3,0.08]);

fun_label(hp_option,'Threshold: Distance','label_track_option_thrdist',[0.05,0.82,0.6,0.08]);
fun_edit( hp_option,'4',                  'edit_track_option_thrdist', [0.65,0.82,0.3,0.08]);

fun_label(hp_option,'Relative Tolerance','label_track_option_tol',[0.05,0.74,0.6,0.08]);
fun_edit( hp_option,'5e-5',              'edit_track_option_tol', [0.65,0.74,0.3,0.08]);

fun_label(hp_option,'Maximum Iteration','label_track_option_maxiter',[0.05,0.66,0.6,0.08]);
fun_edit( hp_option,'1000',             'edit_track_option_maxiter', [0.65,0.66,0.3,0.08]);

% fun_label(hp_option,'LM: Initial Value','label_track_option_lm_init',[0.05,0.68,0.6,0.06]);
% fun_edit( hp_option,'1',                'edit_track_option_lm_init', [0.65,0.68,0.3,0.06]);
%
% fun_label(hp_option,'LM: Damping Factor','label_track_option_lm_damp',[0.05,0.62,0.6,0.06]);
% fun_edit( hp_option,'2',                 'edit_track_option_lm_damp', [0.65,0.62,0.3,0.06]);
%
% fun_label(hp_option,'LM: Maximum Iteration','label_track_option_lm_maxiter',[0.05,0.56,0.6,0.06]);
% fun_edit( hp_option,'100',                  'edit_track_option_lm_maxiter', [0.65,0.56,0.3,0.06]);
%
% fun_label(hp_option,'Weight: Intensity','label_track_option_lambda_pi',[0.05,0.50,0.6,0.06]);
% fun_edit( hp_option,'0',                'edit_track_option_lambda_pi', [0.65,0.50,0.3,0.06]);
%
% fun_label(hp_option,'Weight: Position','label_track_option_lambda_mu_nondiag',[0.05,0.44,0.6,0.06]);
% fun_edit( hp_option,'0',               'edit_track_option_lambda_mu_nondiag', [0.65,0.44,0.3,0.06]);
%
% fun_label(hp_option,'Weight: Position (scale)','label_track_option_lambda_mu_coeff',[0.05,0.38,0.6,0.06]);
% fun_edit( hp_option,'1',                       'edit_track_option_lambda_mu_coeff', [0.65,0.38,0.3,0.06]);
%
% fun_label(hp_option,'Weight: Shape','label_track_option_lambda_sigma',[0.05,0.32,0.6,0.06]);
% fun_edit( hp_option,'0',            'edit_track_option_lambda_sigma', [0.65,0.32,0.3,0.06]);

fun_label(hp_option,'Size limit','label_track_option_ub_sigma',[0.05,0.58,0.6,0.08]);
fun_edit( hp_option,'0',         'edit_track_option_ub_sigma', [0.65,0.58,0.3,0.08]);

fun_label(    hp_option,'Algorithm',          'label_track_option_algorithm',[0.05,0.50,0.6,0.08]);
fun_popupmenu(hp_option,data_calc_algorithm(),'popup_track_option_algorithm',[0.65,0.50,0.3,0.08]);

fun_checkbox(hp_option,'Fix Intensity',       'check_track_option_fixpik',  [0.05,0.44,0.9,0.06]);
fun_checkbox(hp_option,'Fix Position',        'check_track_option_fixmu',   [0.05,0.38,0.9,0.06]);
fun_checkbox(hp_option,'Fix Shape',           'check_track_option_fixvol',  [0.05,0.32,0.9,0.06]);
fun_checkbox(hp_option,'Selected ROI Only',   'check_track_option_selected',[0.05,0.26,0.9,0.06]);
fun_checkbox(hp_option,'Main Parameters Only','check_track_option_main',    [0.05,0.20,0.9,0.06]);

% fun_pushbutton(hp_option,'L','load2',[0.65,0.18,0.15,0.1],'Callback',@data_load_callback); % temporally
% fun_pushbutton(hp_option,'S','save2',[0.80,0.18,0.15,0.1],'Callback',@data_save_callback); % temporally


% fun_pushbutton(hp_option,'Do Track','button_track_option_track',[0.65,0.03,0.3,0.25]),...
%     'Callback',@data_calc_callback);
fun_pushbutton(hp_option,'Do Track','button_track_option_track',[0.05,0.02,0.9,0.16],...
    'Callback',@data_calc_callback);
end


function create_GUI_trackqp(hp)

fun_label(hp,'Run ID', 'label_trackqp_runid',[0.05,0.90,0.30,0.05]);
fun_edit( hp,'trackqp','edit_trackqp_runid', [0.35,0.90,0.60,0.05]);

fun_label(hp,'Frames',        'label_trackqp_frames',[0.05,0.85,0.30,0.05]);
fun_edit(hp,'{10:-1:1,10:20}','edit_trackqp_frames', [0.35,0.85,0.60,0.05]);

fun_label(hp,'WorkDir','label_trackqp_workdir',[0.05,0.80,0.30,0.05]);
fun_edit(hp,'',        'edit_trackqp_workdir', [0.35,0.80,0.60,0.05]);

fun_pushbutton(hp,'Set WorkDir','button_trackqp_set',[0.05,0.70,0.90,0.05],...
    'Callback',@(x,y)trackqp_setworkdir);

fun_pushbutton(hp,'Export to QPipe','button_trackqp_export',[0.05,0.60,0.90,0.05],...
    'Callback',@(x,y)trackqp_export);

fun_pushbutton(hp,'Import from QPipe (SPF)',...
    'button_trackqp_import_spf',[0.05,0.50,0.90,0.05],...
    'Callback',{@(x,y,z)trackqp_import(z),'spf'});

fun_pushbutton(hp,'Import from QPipe (MODETRACK01)',...
    'button_trackqp_import_modetrack01',[0.05,0.40,0.90,0.05],...
    'Callback',{@(x,y,z)trackqp_import(z),'modetrack01'});

end


function create_GUI_annotate(hp)
%%% create panels
hp_atlas = fun_panel(hp,'Atlas',             'panel_anno_atlas',[0.05,0.83,0.9,0.15]);
%hp_bm    = fun_panel(hp,'Auto-Annotation','panel_anno_BM',     [0.05,0.59,0.9,0.22]);
%hp_mv    = fun_panel(hp,'Majority voting',  'panel_anno_MV',   [0.05,0.35,0.9,0.22]);
%hp_feat  = fun_panel(hp,'Features',         'panel_anno_feat', [0.05,0.02,0.9,0.31]);
hp_AA  = fun_panel(hp,'Auto-Annotation',     'panel_anno_AA',   [0.05,0.50,0.9,0.31]);
hp_MA  = fun_panel(hp,'Manual Annotation',   'panel_anno_MA',   [0.05,0.02,0.9,0.44]);

%%% The Atlas Section
% atlas path
fun_edit(hp_atlas,which('atlas_20181230_White=1_MaskPha=1_RemoveHYP=1.mat'),'edit_anno_path',[0.05,0.5,0.9,0.45]);

% num. of atlas for majority voting
fun_label(hp_atlas,'# Atlas:','text_anno_nA',[0.05,0.05,0.5,0.45]);
fun_edit( hp_atlas,'100',     'edit_anno_nA',[0.55,0.05,0.4,0.45]);


%%% The Auto-Annotation Section
% (35): include eat4p feature and its weight
tmp_height = 0.8; % change to 0.8
fun_checkbox(hp_AA,'Label/promoter',    'check_anno_feat_eat4',[0.05,tmp_height,0.5,0.15]);
%fun_text(hp_feat,'w:',           'text_anno_feat_eat4', [0.65,tmp_height,0.10,0.15]);
%fun_edit(hp_feat,'15',           'edit_anno_feat_eat4', [0.75,tmp_height,0.20,0.15]);
% fun_popupmenu(hp_feat,'(column)','popup_anno_feat_eat4',[0.60,0.20,0.35,0.15]);
%fun_popupmenu(hp_feat,'(column)','popup_anno_feat_eat4',[0.30,tmp_height,0.30,0.15]);
% fun_pushbutton(hp_AA,'Details','button_anno_label_details',[0.65,tmp_height,0.3,0.15],...
%     'Callback',@run_labeldetails_callback);
fun_pushbutton(hp_AA,'Details','button_anno_label_details',[0.65,tmp_height,0.3,0.15]);


% (39): include manual update information and its weight
tmp_height = 0.65; % change to 0.65
fun_checkbox(hp_AA,'Manual Labeling','check_anno_feat_human',[0.05,tmp_height,0.50,0.15],'Value',1);
% fun_text(hp_AA,'w:',               'text_anno_feat_human', [0.65,tmp_height,0.10,0.15]);
% fun_edit(hp_AA,'1',                'edit_anno_feat_human', [0.75,tmp_height,0.20,0.15],...
%     'Enable','off','HandleVisibility','off'); % currently not used

% (14): Max. number of rank to be stored
tmp_height = 0.5; % change to 0.5
fun_text(hp_AA,'# Rank storing:','text_anno_MV_rank',[0.05,tmp_height,0.50,0.15]);
fun_edit(hp_AA,'3',              'edit_anno_MV_rank',[0.55,tmp_height,0.40,0.15]);

% (22): number of CPU running in parallel
tmp_height = 0.35; % change to 0.35
fun_text(hp_AA,'# CPU in par.:','text_anno_MV_nCPU',[0.05,tmp_height,0.5,0.15]);
fun_edit(hp_AA,'1',             'edit_anno_MV_nCPU',[0.55,tmp_height,0.4,0.15]);

% run
% callback: @run_annotation
% fun_pushbutton(hp_AA,'Do Annotation','button_anno_run',[0.05,0.05,0.5,0.25],...
%     'Callback',@run_annotation_callback);
fun_pushbutton(hp_AA,'Do Annotation','button_anno_run',[0.05,0.05,0.45,0.25]);

fun_checkbox(hp_AA,'Auto Update','check_anno_auto',[0.55,0.05,0.40,0.15],'Value',0);

%%% The Manual Annotation Section
% approve auto-annotation
tmp_height = 0.75;
fun_popupmenu(hp_MA,'(column)',    'popup_anno_moveneu',[0.05,tmp_height,0.5,0.2]);
fun_pushbutton(hp_MA,'Edit list',    'button_anno_approvelist',[0.6,tmp_height,0.35,0.2]);
%fun_text(hp_feat,'w:',           'text_anno_feat_eat4', [0.65,tmp_height,0.10,0.15]);
%fun_edit(hp_feat,'15',           'edit_anno_feat_eat4', [0.75,tmp_height,0.20,0.15]);
% fun_popupmenu(hp_feat,'(column)','popup_anno_feat_eat4',[0.60,0.20,0.35,0.15]);
%fun_popupmenu(hp_feat,'(column)','popup_anno_feat_eat4',[0.30,tmp_height,0.30,0.15]);

% edit neuron sets for approval
tmp_height = 0.5;
fun_pushbutton(hp_MA,'Approve Anno.',    'button_anno_approve',[0.05,tmp_height,0.45,0.2]);

% % use CPD or not?
% fun_checkbox(hp_bm,'Use CPD','check_anno_BM_CPD', [0.05,0.72,0.9,0.22]);
%
% % CPD - beta
% fun_text(hp_bm,'CPD - beta:','text_anno_BM_beta',[0.05,0.50,0.5,0.22]);
% fun_edit(hp_bm,'2',          'edit_anno_BM_beta',[0.55,0.50,0.4,0.22]);
%
% % CPD - lambda
% fun_text(hp_bm,'CPD - lambda:','text_anno_BM_lambda',[0.05,0.28,0.5,0.22]);
% fun_edit(hp_bm,'3',            'edit_anno_BM_lambda',[0.55,0.28,0.4,0.22]);
%
% % CPD - outlier
% fun_text(hp_bm,'CPD - outlier:','text_anno_BM_outlier',[0.05,0.06,0.5,0.22]);
% fun_edit(hp_bm,'0.3',           'edit_anno_BM_outlier',[0.55,0.06,0.4,0.22]);
%
%
% %%% The Majority Voting Section
% % (14): Max. number of rank to be stored
% fun_text(hp_mv,'# Rank storing:','text_anno_MV_rank',[0.05,0.72,0.50,0.22]);
% fun_edit(hp_mv,'3',              'edit_anno_MV_rank',[0.55,0.72,0.40,0.22]);
%
% % (16): use Rank-feedback or not and its weight
% fun_checkbox(hp_mv,'Rank-feedback','check_anno_MV_rank', [0.05,0.50,0.60,0.22]);
% fun_text(hp_mv,'w:',               'text_anno_MV_rank_w',[0.65,0.50,0.10,0.22]);
% fun_edit(hp_mv,'0.2',              'edit_anno_MV_rank_w',[0.75,0.50,0.20,0.22]);
%
% % (20): save W for speedup in Rank-feedback or not
% fun_checkbox(hp_mv,'Speedup (save W)','check_anno_MV_saveW', [0.05,0.28,0.9,0.22]);
%
% % (22): number of CPU running in parallel
% fun_text(hp_mv,'# CPU in par.:','text_anno_MV_nCPU',[0.05,0.07,0.5,0.20]);
% fun_edit(hp_mv,'1',             'edit_anno_MV_nCPU',[0.55,0.07,0.4,0.20]);
%
%
% % Feature Section
% % (39): include manual update information and its weight
% fun_checkbox(hp_feat,'Manual Update','check_anno_feat_human',[0.05,0.8,0.50,0.15]);
% fun_text(hp_feat,'w:',               'text_anno_feat_human', [0.65,0.8,0.10,0.15]);
% fun_edit(hp_feat,'1',                'edit_anno_feat_human', [0.75,0.8,0.20,0.15],...
%     'Enable','off','HandleVisibility','off'); % currently not used
%
% % (31): include neuron intensity feature and its weight
% fun_checkbox(hp_feat,'Intensity','check_anno_feat_int',[0.05,0.65,0.60,0.15],...
%     'Enable','off','HandleVisibility','off'); % currently not implemented
% fun_text(hp_feat,'w:',           'text_anno_feat_int', [0.65,0.65,0.10,0.15]);
% fun_edit(hp_feat,'0.1',          'edit_anno_feat_int', [0.75,0.65,0.20,0.15],...
%     'Enable','off','HandleVisibility','off'); % currently not implemented
%
% % (27): include neuron size feature and its weight
% fun_checkbox(hp_feat,'Shape','check_anno_feat_shape',[0.05,0.5,0.60,0.15],...
%     'Enable','off','HandleVisibility','off'); % currently not implemented
% fun_text(hp_feat,'w:',       'text_anno_feat_shape', [0.65,0.5,0.10,0.15]);
% fun_edit(hp_feat,'0.1',      'edit_anno_feat_shape', [0.75,0.5,0.20,0.15],...
%     'Enable','off','HandleVisibility','off'); % currently not implemented
%
% % (35): include eat4p feature and its weight
% fun_checkbox(hp_feat,'eat4p',    'check_anno_feat_eat4',[0.05,0.35,0.25,0.15]);
% fun_text(hp_feat,'w:',           'text_anno_feat_eat4', [0.65,0.35,0.10,0.15]);
% fun_edit(hp_feat,'15',           'edit_anno_feat_eat4', [0.75,0.35,0.20,0.15]);
% % fun_popupmenu(hp_feat,'(column)','popup_anno_feat_eat4',[0.60,0.20,0.35,0.15]);
% fun_popupmenu(hp_feat,'(column)','popup_anno_feat_eat4',[0.30,0.35,0.30,0.15]);
%
% % run
% % callback: @run_annotation
% fun_pushbutton(hp_feat,'Do Annotation','button_anno_run',[0.05,0.05,0.5,0.25],...
%     'Callback',@run_annotation_callback);

roiedit3d_anno_v3_mod7(); % initialize

end


function create_GUI_analyze(hp)

hpn = fun_panel(hp,'Normal Object', 'panel_analyze_normal', [0.05,0.84,0.9,0.14]);
hps = fun_panel(hp,'Special Object','panel_analyze_special',[0.05,0.29,0.9,0.53]);
hpo = fun_panel(hp,'Others',        'panel_analyze_others', [0.05,0.02,0.9,0.25]);

hpt = fun_panel(hps,'Target',    'panel_analyze_target',    [0,0.47,1,0.50]);
hpc = fun_panel(hps,'Clustering','panel_analyze_clustering',[0,0.03,1,0.40]);


%%% plot target
% nomal
fun_label(  hpn,'Plot',  'label_analyze_plot',  [0.05,0.05,0.3,0.9]);
fun_listbox(hpn,'(none)','listbox_analyze_plot',[0.35,0.05,0.6,0.9],...
    'Max',2); % multiple selection enabled

% numerator
fun_label(    hpt,'Special',   'label_analyze_numerator',[0.05,0.73,0.45,0.23]);
fun_popupmenu(hpt,'(none)',    'popup_analyze_numerator',[0.50,0.73,0.45,0.23]);

% denominator
fun_label(    hpt,'divided by','label_analyze_denominator',[0.05,0.50,0.45,0.23]);
fun_popupmenu(hpt,'(none)',    'popup_analyze_denominator',[0.50,0.50,0.45,0.23]);

% smoothing by Savitzky-Golay algorithm
fun_checkbox(hpt, 'Smoothing','check_analyze_smoothing',[0.05,0.27,0.45,0.23]);
fun_edit(    hpt, '31',       'edit_analyze_smoothing', [0.50,0.27,0.45,0.23]);

% normalize to 0-1
fun_checkbox(hpt, 'Normalize','check_analyze_normalize',[0.05,0.04,0.45,0.23]);


%%% clustering
strPopupMetric = {'euclidean','seuclidean','cityblock','minkowski','chebychev',...
    'mahalanobis','cosine','correlation','spearman','hamming','jaccard',...
    'abs_cosine','abs_correlation','abs_spearman','abs_jaccard'};
strPopupMethod = {'ward','average','centroid','complete','median','single','weighted'};

fun_checkbox(hpc,'Clustering','check_analyze_clustering',[0.05,0.65,0.45,0.3]);

fun_label(    hpc,'Metric',      'label_analyze_metric',[0.05,0.35,0.45,0.3]);
fun_popupmenu(hpc,strPopupMetric,'popup_analyze_metric',[0.50,0.35,0.45,0.3]);

fun_label(    hpc,'Method',      'label_analyze_method',[0.05,0.05,0.45,0.3]);
fun_popupmenu(hpc,strPopupMethod,'popup_analyze_method',[0.50,0.05,0.45,0.3]);


%%% other settings
% fps
fun_label(    hpo,'Volume/sec', 'label_analyze_vps',      [0.05,0.72,0.45,0.22]);
fun_edit(     hpo,'4.85',       'edit_analyze_vps',       [0.50,0.72,0.45,0.22]);

% show plot of target movement (alignment amount)
fun_checkbox(hpo, 'Show Move','check_analyze_show_move',  [0.05,0.50,0.45,0.22]);

% show plot of target movement (alignment amount)
fun_checkbox(hpo,'Merge plot','check_analyze_merge',      [0.05,0.28,0.45,0.22]);

% Show heatmap
fun_pushbutton(hpo,'Show heatmap','button_analyze_show',  [0.05,0.06,0.45,0.22],...
    'Callback',@(x,y)createHeatmap);

% plot timecourses
fun_pushbutton(hpo,'Plot target', 'button_analyze_plot',  [0.50,0.06,0.45,0.22],...
    'Callback',@(x,y)plotTcrs);

% set(findobj('-regexp','Tag','label_analyze_*','-regexp'),'enable','inactive');

end


function create_GUI_vistrack(hp)
fun_pushbutton(hp,'Initialize / Update','button_vistrack_init',[0.05,0.94,0.45,0.04],...
    'Callback',@(x,y)visTrackInit());
end


function create_GUI_movie(hp)

fun_label(hp,'Scale', 'label_movie_zoom', [0.10,0.90,0.45,0.05]);
fun_edit( hp,'1',      'edit_movie_zoom', [0.55,0.90,0.35,0.05],...
    'TooltipString','Scale of ROIs for drawing');
fun_label(hp,'Start','label_movie_start', [0.10,0.85,0.45,0.05]);
fun_edit( hp,'1',     'edit_movie_start', [0.55,0.85,0.35,0.05],...
    'TooltipString','Start frame of max projected movie');
fun_label(hp,'End',    'label_movie_end', [0.10,0.80,0.45,0.05]);
fun_edit( hp,'100',     'edit_movie_end', [0.55,0.80,0.35,0.05],...
    'TooltipString','End frame of max projected movie');

fun_checkbox(hp,'Max Projected ROIs','check_movie_roi',                        [0.10,0.75,0.80,0.05]);
fun_checkbox(hp,'Max Projected Movie','check_movie_movie',                     [0.10,0.70,0.80,0.05]);
fun_checkbox(hp,'Max Projected Movie with ROIs Overlaid','check_movie_overlay',[0.10,0.65,0.80,0.05]);

fun_pushbutton(hp,'Save Max Projected ROIs & Movies','button_movie_save',[0.10,0.60,0.80,0.05],...
    'Callback',@(x,y)saveMaxProj);

end


function create_GUI_option(hp)

fun_label(hp,'Visualize radius','label_option_visualize',[0.10,0.90,0.45,0.05]);
fun_edit( hp,'1',               'edit_option_visualize', [0.55,0.90,0.35,0.05],...
    'Callback',@(~,~)setNsigma);

fun_label(hp,'Dist. threshold (um)','label_option_check',[0.10,0.85,0.45,0.05]);
fun_edit( hp,'10',                 'edit_option_check',  [0.55,0.85,0.35,0.05],...
    'Callback',@checkbox_mindist_callback);

fun_label(hp,'Scaling (um/pix)', 'label_option_scaling',[0.10,0.80,0.45,0.05]);
fun_edit( hp,'[0.24,0.24,0.252]','edit_option_scaling', [0.55,0.80,0.35,0.05],...
    'Callback',@(x,y)setScale);

fun_label(hp,'Default Sigma','label_option_sigma',[0.10,0.75,0.45,0.05]);
fun_edit( hp,'[10,10,10]',   'edit_option_sigma', [0.55,0.75,0.35,0.05],...
    'TooltipString','Default value of [x,y,z] (pixels) of sigma of elipsoid (multivariate gaussian distribution)');


fun_label(hp,'MaxProj ROI scale','label_option_maxproj',[0.10,0.65,0.45,0.05]);
fun_edit( hp,'1',                'edit_option_maxproj', [0.55,0.65,0.35,0.05],...
    'TooltipString','Scaling of max projected rois for save');

fun_pushbutton(hp,'Save MaxProj ROI','button_option_save',[0.10,0.60,0.80,0.05],...
    'Callback',@(x,y)saveMaxProj);


% fun_label(hp,'Detection Channel','label_option_detection',[0.10,0.70,0.45,0.05]);
% fun_edit( hp,'1',                'edit_option_detection', [0.55,0.70,0.35,0.05],...
%     'TooltipString','When there is no blurred image (*_blurred.tif), this channel will be used for ROI detection');

% fun_pushbutton(hp,'Incoherency','button_option_incoherency',[0.05,0.05,0.40,0.05],...
%     'Callback',@(x,y)updateIncoherency);
%
% fun_pushbutton(hp,'Plot 4D','button_option_4D',[0.55,0.05,0.40,0.05],...
%     'Callback',@(x,y)visualizeTrack);
%
% fun_pushbutton(hp,'Check track','button_option_check',[0.05,0.05,0.40,0.05],...
%     'Callback',@(x,y)visTrackInit);


fun_label(hp,'Cache size','label_option_cache',[0.10,0.50,0.45,0.05]);
fun_edit( hp,'5000',      'edit_option_cache', [0.55,0.50,0.35,0.05],...
    'Callback',@(~,~)setCacheSize,'TooltipString','Size of image cache');


fun_pushbutton(hp,'visualize promoter expression','button_option_vis_prom_exp',[0.05,0.05,0.40,0.05],...
    'Callback',@(x,y)visualize_promoter_expression_13);





end


%%% data callback functions


function data_load_callback(~,~,strin,flagAuto)

%%% button 1 ：
%%% open mat file and obtain image & roi data
%%% display image in ImageJ
%%% display roi data over the image
%%% set data to customized roi manager and display it
%%% reset listeners
%%% enable orthogonal view
%%% obtain image plus and canvas, then set canvas_callback
%%%

userdata = getUserData();

%%% clear last session
strPathName = [];
if isfield(userdata,'frame') || isfield(userdata,'himp')
    if ~exist('flagAuto','var') || isempty(flagAuto) || ~flagAuto
        answer = questdlg('Do you close current session?',...
            'Closing session','Yes','No','No');
    else
        answer = 'Yes';
    end
    drawnow();
    if strcmpi(answer,'Yes')
        strPathName = userdata.PathName;
        if isfield(userdata,'frame'); userdata.frame.dispose(); end
        if isfield(userdata,'himp'); userdata.himp.changes = false; end
        if isfield(userdata,'himp_proj'); userdata.himp_proj.changes = false; end
        ij.IJ.run('Close All');
        close(gcf);
        close(findobj('-regexp','Tag','figure_*'));
        clear userdata;
        
        %%% create new GUI
        create_GUI();
    else
        return;
    end
end

enableGUI(false);

%%% get data file path
if exist('strin','var') && ~isempty(strin)
    [PathName,FileName,ext] = fileparts(strin);
    FileName = [FileName,ext];
else
    [FileName,PathName] = uigetfile('*.mat','Select detection result',strPathName);
    if  isequal(FileName,0) || isequal(PathName,0) % canceled
        enableGUI(false);
        return;
    end
end
disp(['loading started: ',fullfile(PathName,FileName)]);

%%% loading metadata
% clear userdata;
userdata = getUserData();
userdata.FileName = FileName;
userdata.PathName = PathName;
setUserData(userdata);
data_load_metadata();

%%% load images
data_load_aligned_image(); % load aligned image
data_load_blurred_image(); % load blurred image

%%% load table
data_load_table();
table_createGUI();
initAnalyze();

%%% initialize annotation part
roiedit3d_anno_v3_mod7();

%%% display image
userdata = getUserData();
himp = userdata.himp;

himp.show();
himp.resetStack();


%%% enable orthogonal views
if himp.getNSlices()>1 % multi z-slice image
    while (isempty(himp.getWindow()))
        pause(0.1);
    end
    ij.WindowManager.setCurrentWindow(himp.getWindow());
    hmov = javaObjectEDT(MyOrthogonal_Views());
    hmov.run('');
    while (true)
        himp_xz = hmov.getXZImage();
        himp_yz = hmov.getYZImage();
        if     ~isempty(himp_xz) && ~isempty(himp_xz.getWindow()) ...
                && ~isempty(himp_yz) && ~isempty(himp_yz.getWindow())
            break;
        end
        pause(0.1);
    end
end

userdata.himp = himp;
setUserData(userdata);


%%% visualize roi
init_ellipsoid_roi();


%%% setup Canvas
userdata = getUserData();
userdata.hcanvas_callback_semaphore = userdata.htable_callback_semaphore;
userdata.hcanvas = userdata.himp.getCanvas();
setUserData(userdata);

addCanvasListener();

himp.updateAndDraw();

enableGUI(true);

disp('loading finished.');
disp(' ');

end


function data_load_metadata()
%%% inner function of data_load_callback

userdata = getUserData();
FileName = userdata.FileName;
PathName = userdata.PathName;

%%% load data
warning('off','MATLAB:load:variableNotFound');
data_loaded = load(fullfile(PathName,FileName),...
    'scaling','sigma_default','diameter_check','strStrokeColor','strSelectionColor',...
    'imname','emopt','cum2','flag_virtual_stack','rotinfo','idxColDisp','myluts',...
    'flagAlign','flagRotateX','flagRotateY','flagRotateZ','flagLn','flagInterpZ',...
    'flagSubBG','rotation','flagMirror','imext','detopt','cellColumnHeader',...
    'cellColumnWidth','para_A','sizCache');
warning('on','MATLAB:load:variableNotFound');

% set scaling if data_loaded contains it
if isfield(data_loaded,'scaling')
    set(findobj('Tag','edit_option_scaling'),...
        'String',sprintf('[%g,%g,%g]',data_loaded.scaling));
end

% set sigma_default if data_loaded contains it
if isfield(data_loaded,'sigma_default')
    set(findobj('Tag','edit_option_sigma'),...
        'String',sprintf('[%g,%g,%g]',data_loaded.sigma_default));
end

% update diameter_check
scaling       = str2num(get(findobj('Tag','edit_option_scaling'),'String')); %#ok<ST2NM>
sigma_default = str2num(get(findobj('Tag','edit_option_sigma'),  'String')); %#ok<ST2NM>
set(findobj('Tag','edit_option_check'),'String',...
    num2str(sqrt(mean(sigma_default.*(scaling.^2)))*3,'%.1f'));

% set diameter_check if data_loaded contains it
if isfield(data_loaded,'diameter_check')
    set(findobj('Tag','edit_option_check'),'String',...
        num2str(data_loaded.diameter_check));
end

% set cache size if data_loaded containt it
if isfield(data_loaded,'sizCache')
    set(findobj('Tag','edit_option_cache'),'String',...
        num2str(data_loaded.sizCache));
end

% % set channel_detection if data_loaded contains it
% if isfield(data_loaded,'channel_detection')
%     set(findobj('Tag','edit_option_detection'),'String',...
%         num2str(data_loaded.channel_detection));
% end

% set strokecolor if data_loaded contains it
if isfield(data_loaded,'strStrokeColor')
    huic_strokecolor = findobj('Tag','popup_main_strokecolor');
    strColors = get(huic_strokecolor,'String');
    tmpidxc = find(strcmp(data_loaded.strStrokeColor,strColors));
    if ~isempty(tmpidxc)
        set(huic_strokecolor,'Value',tmpidxc);
    end
end

% set selectioncolor if data_loaded contains it
if isfield(data_loaded,'strSelectionColor')
    huic_selectioncolor = findobj('Tag','popup_main_selectioncolor');
    strColors = get(huic_selectioncolor,'String');
    tmpidxc = find(strcmp(data_loaded.strSelectionColor,strColors));
    if ~isempty(tmpidxc)
        set(huic_selectioncolor,'Value',tmpidxc);
    end
end

% set tracking option with loaded data
if isfield(data_loaded,'emopt')
    setEmopt(data_loaded.emopt);
end

% set checkbox for logarithm intensity
if isfield(data_loaded,'flagLn')
    set(findobj('Tag','check_main_logarithm'),...
        'Value',data_loaded.flagLn);
end

% set checkbox for Alignment
if isfield(data_loaded,'flagAlign')
    set(findobj('Tag','check_main_align'),...
        'Value',data_loaded.flagAlign);
end

% set checkbox for Z interpolarion
if isfield(data_loaded,'flagInterpZ')
    set(findobj('Tag','check_main_interpz'),...
        'Value',data_loaded.flagInterpZ);
end

% set checkbox for subtract background
if isfield(data_loaded,'flagSubBG')
    set(findobj('Tag','check_main_subtract'),...
        'Value',data_loaded.flagSubBG);
end

% set checkbox for Mirror Image
if isfield(data_loaded,'rotation')
    set(findobj('Tag','edit_main_rotate_arbit'),...
        'String',num2str(data_loaded.rotation/pi*180));
    set(findobj('Tag','slider_main_rotate_arbit'),...
        'Value',data_loaded.rotation/pi*180);
end

% set checkbox for Mirror Image
if isfield(data_loaded,'flagMirror')
    set(findobj('Tag','check_main_rotate_mirror'),...
        'Value',data_loaded.flagMirror);
end

if isfield(data_loaded,'flagRotateX')
    set(findobj('Tag','button_main_rotate_x'),...
        'Value',data_loaded.flagRotateX);
end
if isfield(data_loaded,'flagRotateY')
    set(findobj('Tag','button_main_rotate_y'),...
        'Value',data_loaded.flagRotateY);
end
if isfield(data_loaded,'flagRotateZ')
    set(findobj('Tag','button_main_rotate_z'),...
        'Value',data_loaded.flagRotateZ);
end

%%% set up userdata
userdata.imname = data_loaded.imname;
userdata.imext = '.tif';
userdata.cum2 = [];
userdata.flag_virtual_stack = true;
userdata.rotinfo = struct();
userdata.myluts = [];
userdata.detopt = getDetectOpt();

if isfield(data_loaded,'imext')
    userdata.imext = data_loaded.imext;
end
if isfield(data_loaded,'cum2')
    userdata.cum2 = data_loaded.cum2;
end
if isfield(data_loaded,'flag_virtual_stack')
    userdata.flag_virtual_stack = data_loaded.flag_virtual_stack;
end
if isfield(data_loaded,'rotinfo')
    userdata.rotinfo = data_loaded.rotinfo;
end
if isfield(data_loaded,'myluts')
    userdata.myluts = data_loaded.myluts;
end
if isfield(data_loaded,'flagAlign')
    userdata.flagAlign = data_loaded.flagAlign;
end
if isfield(data_loaded,'detopt')
    userdata.detopt = data_loaded.detopt;
    setDetectOpt(userdata.detopt);
end
if isfield(data_loaded,'idxColDisp')
    userdata.idxColDisp = data_loaded.idxColDisp;
end
if isfield(data_loaded,'cellColumnHeader')
    userdata.cellColumnHeader = data_loaded.cellColumnHeader;
end
if isfield(data_loaded,'cellColumnWidth')
    userdata.cellColumnWidth = data_loaded.cellColumnWidth;
end
if isfield(data_loaded,'para_A')
    userdata.para_A = data_loaded.para_A;
end

setUserData(userdata);

end


function data_load_aligned_image()
userdata = getUserData();
PathName = userdata.PathName;
FileName = userdata.FileName;
path_aligned = fullfile(PathName,[userdata.imname,'_aligned.tif']);
disp(['loading aligned image: ',path_aligned]);

if ~exist(path_aligned,'file')
    warning('off','MATLAB:load:variableNotFound');
    data_loaded = load(fullfile(PathName,FileName),'im_aligned');
    warning('on','MATLAB:load:variableNotFound');
    if isfield(data_loaded,'im_aligned')
        himp = setImageMiji(data_loaded.im_aligned,[userdata.imname,'_aligned']);
        ij.IJ.saveAsTiff(himp,path_aligned);
    else
        %         path_aligned = fullfile(PathName,[userdata.imname,'.tif']);
        path_aligned = fullfile(PathName,[userdata.imname,userdata.imext]);
        if ~exist(path_aligned,'file')
            error('Aligned image not found');
        end
    end
end

[~,himp,mfivs] = getImageMiji(path_aligned,1,1,1);

userdata.mfivs = mfivs;
userdata.himp = himp;
setUserData(userdata);

if himp.getNChannels>1
    himp = ij.CompositeImage(himp);
    himp.setMode(ij.CompositeImage.COMPOSITE);
end
setScale(); % set scaling parameter
setLuts(userdata.himp); % recover luts and apply; 'userdata.himp' is required for apply luts
setCacheSize();

if get(findobj('Tag','button_main_rotate_x'),'Value')==1
    mfivs.toggleFlipY(); mfivs.toggleFlipZ();
end
if get(findobj('Tag','button_main_rotate_y'),'Value')==1
    mfivs.toggleFlipX(); mfivs.toggleFlipZ();
end
if get(findobj('Tag','button_main_rotate_z'),'Value')==1
    mfivs.toggleFlipX(); mfivs.toggleFlipY();
end
if get(findobj('Tag','check_main_rotate_mirror'),'Value')==1
    mfivs.toggleFlipZ();
end
if cellfun(@isempty,regexp(cell(userdata.himp.getTitle()),'_aligned\.tif$','once'))
    % original image (not aligned) is shown
    mfivs.setParallelShift(userdata.cum2);
    mfivs.setAlign(get(findobj('Tag','check_main_align'),'Value')==1);
else % aligned image is shown
    if isfield(userdata,'flagAlign')
        set(findobj('Tag','check_main_align'),'Value',userdata.flagAlign);
        mfivs.setParallelShift(-userdata.cum2);
        mfivs.setAlign(~userdata.flagAlign);
        userdata = rmfield(userdata,'flagAlign');
        setUserData(userdata);
    else % for the first time
        set(findobj('Tag','check_main_align'),'Value',1);
        mfivs.setParallelShift(-userdata.cum2);
        mfivs.setAlign(false);
    end
end
mfivs.setLn(get(findobj('Tag','check_main_logarithm'),'Value')==1);
mfivs.setInterpZ(get(findobj('Tag','check_main_interpz'),'Value')==1);
mfivs.setSubBG(get(findobj('Tag','check_main_subtract'),'Value')==1);
rad = str2double(get(findobj('Tag','edit_main_rotate_arbit'),'String'))/180*pi;
mfivs.setRotation(rad);

end


function data_load_blurred_image()
userdata = getUserData();
PathName = userdata.PathName;
FileName = userdata.FileName;
path_blurred  = fullfile(PathName,[userdata.imname,'_blurred.tif']);
path_blurred2 = fullfile(PathName,[userdata.imname,'_thresholded.tif']); % for version 6
disp(['loading blurred image: ',path_blurred]);

if ~exist(path_blurred,'file')
    if exist(path_blurred2,'file')
        path_blurred = path_blurred2;
    else % try recover from im_blurred_nn (for old data)
        warning('off','MATLAB:load:variableNotFound');
        data_loaded = load(fullfile(PathName,FileName),'im_blurred_nn');
        warning('on','MATLAB:load:variableNotFound');
        if isfield(data_loaded,'im_blurred_nn')
            himp_blurred = setImageMiji(...
                data_loaded.im_blurred_nn,[userdata.imname,'_blurred']);
            ij.IJ.saveAsTiff(himp_blurred,path_blurred);
        else
            %             path_blurred = fullfile(PathName,[userdata.imname,'.tif']);
            path_blurred = fullfile(PathName,[userdata.imname,userdata.imext]);
            if ~exist(path_blurred,'file')
                path_blurred = fullfile(PathName,[userdata.imname,'_aligned.tif']);
                if ~exist(path_blurred,'file')
                    error('Blurred image not found');
                end
            end
        end
    end
end

[~,himp_blurred,mfivs_blurred] = getImageMiji(path_blurred,1,1,1);

userdata.mfivs_blurred = mfivs_blurred;
userdata.himp_blurred = himp_blurred;
setUserData(userdata);

end


function data_load_table()

%%% reconstruct table from dat file
userdata = getUserData();
strFileMat = fullfile(userdata.PathName,userdata.FileName);
strFileDat = [strFileMat(1:end-4),'.dat'];
disp(['Trying to load the data table from: ',strFileDat]);
if exist(strFileDat,'file')
    try % current version
        hmtm = MyTableModel4.load(strFileDat);
    catch
        try % old version
            hmtm_old = MyTableModel2.load(strFileDat);
            hmtm = MyTableModel4(hmtm_old.getDataVector(),hmtm_old.getColumnIdentifiers());
        catch % dat file may be created by old verion of this program
            try
                hois = java.io.ObjectInputStream(java.io.FileInputStream(strFileDat));
                hdv = hois.readObject();
                hois.close();
                data_loaded = load(strFileMat,'columnNames');
                hmtm = MyTableModel4(hdv,data_loaded.columnNames);
                
            catch % file is exist but not relevant to this program
                % do nothing
            end
        end
    end
end


%%% reconstruct table from mat file
if ~exist('hmtm','var')
    disp('loading mat file; please wait... ')
    warning('off','MATLAB:load:variableNotFound');
    data_loaded = load(strFileMat,'data_edited','columnNames');
    warning('on','MATLAB:load:variableNotFound');
    
    if ~isfield(data_loaded,'data_edited') || ~isfield(data_loaded,'columnNames')
        %%% data_edited not found; reconstruct table from other variables
        % data_loaded = load(strFileMat,'params_em','flag_oob','flag_repaired','idx_frame');
        % numrows = numel(data_loaded.flag_oob);
        data_loaded = load(strFileMat,'params_em','idx_frame');
        numrows = size(data_loaded.params_em,2);
        uid = num2cell((1:numrows)');
        % data_loaded.columnNames = {...
        %     'Check','uniqueID','Name','frame','OoB','Repaired','min_dist','Name_min_dist',...
        %     'pi_k','xc','yc','zc','S_11','S_12','S_13','S_22','S_23','S_33'};
        data_loaded.columnNames = {...
            'Check','uniqueID','Name','frame','min_dist','Name_min_dist',...
            'pi_k','xc','yc','zc','S_11','S_12','S_13','S_22','S_23','S_33'};
        data_loaded.data_edited = cat(2,...
            num2cell(false(numrows,1)),...
            uid,...
            cellfun(@num2str,uid,'UniformOutput',false),...
            num2cell(data_loaded.idx_frame'),...
            ... % num2cell(data_loaded.flag_oob'),...
            ... % num2cell(data_loaded.flag_repaired'),...
            num2cell(zeros(numrows,1)),...
            repmat({'0'},numrows,1),...
            num2cell(data_loaded.params_em'));
    end
    hmtm = MyTableModel4(data_loaded.data_edited,data_loaded.columnNames);
end


%%% define required data model
numcolor = userdata.himp.getNChannels();
str_numcolor = cellfun(@num2str,num2cell(1:numcolor)','UniformOutput',false);
str_mu = {'xc','yc','zc'};
str_sigma = {'S_11','S_12','S_13','S_22','S_23','S_33'};
str_pi_k_ch = strcat('pi_k_Ch',str_numcolor)';
str_pi_k_ch_normres = strcat('pi_k_Ch',str_numcolor,'_normres')';
str_landmarks = {'Landmark1','Landmark2','Landmark3','Landmark4','Landmark5'};
str_incoherency = {'xd','yd','zd','Incoherency'};

strRequiredColNames = cat(2,{'Unsure','Check','uniqueID','Name','frame','OoB',...
    'Repaired','min_dist','Name_min_dist','pi_k'},str_mu,str_sigma,...
    {'vol','min_scdist','Name_min_scdist'},str_pi_k_ch,str_pi_k_ch_normres,...
    str_landmarks,{'Name_log','Name_estim','Name_estim_unique',...
    'Name_anno_human'},str_incoherency);

strregexp_string = '^Name.*'; % Data in the column is String
strregexp_bool = '^(Unsure|Check|OoB|Repaired|Landmark.*)$'; % Data in the column is boolean

% register keys that indicate multi columns
hmtm.setMap('mu',str_mu);
hmtm.setMap('sigma',str_sigma);
hmtm.setMap('params_em',cat(2,{'pi_k'},str_mu,str_sigma));
hmtm.setMap('pi_k_Ch_all',cat(2,str_pi_k_ch,str_pi_k_ch_normres));
hmtm.setMap('Landmarks',str_landmarks);
hmtm.setMap('Name_popup',{'Name','Name_estim','Name_estim_unique','Name_anno_human'});


%%% establish consistency: whether the table contains required column
for p=1:numel(strRequiredColNames)
    tmpName = strRequiredColNames{p};
    if (hmtm.getColumnIndex(tmpName)==-1) % required column not found
        numrows = hmtm.getRowCount();
        switch tmpName
            case 'frame'
                tmpData = num2cell(ones(1,numrows));
            case 'uniqueID'
                tmpData = num2cell(1:numrows);
            case 'Name'
                tmpData = cellfun(@num2str,num2cell(1:numrows),'UniformOutput',false);
            otherwise
                if ~isempty(regexp(tmpName,strregexp_string, 'once')) % string
                    tmpData = repmat({'0'},[1,numrows]);
                    %                     tmpData = repmat({' '},[1,numrows]);
                elseif ~isempty(regexp(tmpName,strregexp_bool, 'once')) % boolean
                    tmpData = num2cell(false(1,numrows));
                else
                    tmpData = num2cell(zeros(1,numrows));
                end
        end
        hmtm.addColumn(tmpName,tmpData);
    end
end


%%% establish consistency: create default data for matlab
hmtm.addAutoConverterCharToString();
hmtm.convertCharToString();
numcol = hmtm.getColumnCount();
defaultdata = cell(1,numcol);
for p=1:numcol
    elem = hmtm.getValueAt(0,p-1);
    if ischar(elem); defaultdata(p) = {'0'};
    elseif isnumeric(elem); defaultdata(p) = {0};
    elseif islogical(elem); defaultdata(p) = {false};
    end
end

% %%% add choice for column of eat4p
% hci = hmtm.getColumnIdentifiers();
% strCol = cell(hci.toArray());
% set(findobj('Tag','popup_anno_feat_eat4'),'String',strCol(cellfun(@islogical,defaultdata)));


%%% check flags and settings for image filtering
%%% set filtered values and convert to raw values.
idxRowFull = 0:hmtm.getRowCount()-1;
params_em = hmtm.getValuesAt(idxRowFull,'params_em'); % filtered values

mfivs = userdata.mfivs;
hmtm.nxorig = mfivs.nxorig;
hmtm.nyorig = mfivs.nyorig;
hmtm.nzorig = mfivs.nzorig;
hmtm.ncorig = mfivs.ncorig;
hmtm.ntorig = mfivs.ntorig;
hmtm.zscale = mfivs.zscale;

if get(findobj('Tag','button_main_rotate_x'),'Value')==1
    hmtm.toggleFlipY(); hmtm.toggleFlipZ();
end
if get(findobj('Tag','button_main_rotate_y'),'Value')==1
    hmtm.toggleFlipX(); hmtm.toggleFlipZ();
end
if get(findobj('Tag','button_main_rotate_z'),'Value')==1
    hmtm.toggleFlipX(); hmtm.toggleFlipY();
end
if get(findobj('Tag','check_main_rotate_mirror'),'Value')==1
    hmtm.toggleFlipZ();
end

if isempty(userdata.cum2)
    hmtm.setParallelShift(zeros(2,mfivs.nzorig,mfivs.ntorig));
else
    hmtm.setParallelShift(userdata.cum2);
end
hmtm.setAlign(get(findobj('Tag','check_main_align'),'Value')==1);
hmtm.setInterpZ(get(findobj('Tag','check_main_interpz'),'Value')==1);
rad = str2double(get(findobj('Tag','edit_main_rotate_arbit'),'String'))/180*pi;
hmtm.setRotation(rad);
hmtm.setValuesAt(params_em,idxRowFull,'params_em',[false,false,false]);


%%%
hci = hmtm.getColumnIdentifiers();
userdata.columnNames = cell(hci.toArray());
userdata.defaultdata = defaultdata;
userdata.hmtm = hmtm;
setUserData(userdata);

end


function data_save_callback(~,~,flag_auto,pathSave)

%%% get save file path
userdata = getUserData();
[~, tmpname, ~] = fileparts(userdata.FileName);
if exist('flag_auto','var') && ~isempty(flag_auto) && flag_auto % for autosave
    tmpname = regexprep(tmpname,'_tmp$','');
    pathSaveMat = fullfile(userdata.PathName,[tmpname,'_tmp.mat']);
    pathSaveDat = fullfile(userdata.PathName,[tmpname,'_tmp.dat']);
    idx_filter = 2; % fast mode
elseif exist('pathSave','var') && ~isempty(pathSave)
    pathSaveMat = pathSave;
    pathSaveDat = [pathSaveMat(1:end-4),'.dat'];
    idx_filter = 1; % slow mode
else
    tmpname = regexprep(tmpname,'_edited$','');
    defaultSavePathMat = fullfile(userdata.PathName,[tmpname,'_edited.mat']);
    % slow mode is default settings when backward compatibility is checked
    if get(findobj('Tag','check_main_compatibility'),'Value')==1
        [FileName,PathName,idx_filter] = ...
            uiputfile({'*.mat','slow mode';'*.mat','fast mode'},'Save',defaultSavePathMat);
    else
        [FileName,PathName,idx_filter] = ...
            uiputfile({'*.mat','fast mode';'*.mat','slow mode'},'Save',defaultSavePathMat);
        idx_filter = 3 - idx_filter; % invert index
    end
    if  isequal(FileName,0) || isequal(PathName,0)
        return;
    end
    pathSaveMat = fullfile(PathName,FileName);
    pathSaveDat = fullfile(PathName,strrep(FileName,'.mat','.dat'));
end
disp(['saving data as: ',pathSaveMat,' and .dat']);

enableGUI(false);

%%% get ROI colors
huic_strokecolor    = findobj('Tag','popup_main_strokecolor');
huic_selectioncolor = findobj('Tag','popup_main_selectioncolor');
strColors1 = get(huic_strokecolor,   'String');
strColors2 = get(huic_selectioncolor,'String');
strStrokeColor    = strColors1{get(huic_strokecolor,   'Value')}; %#ok<NASGU>
strSelectionColor = strColors2{get(huic_selectioncolor,'Value')}; %#ok<NASGU>

%%% get column display index
idxColDisp = zeros(1,userdata.hmtm.getColumnCount());
for p=1:numel(idxColDisp)
    idxColDisp(p) = userdata.jtable.convertColumnIndexToView(p-1);
end

%%% get column width and header name
cellColumnHeader = getColumnHeader(); %#ok<NASGU>
cellColumnWidth  = getColumnWidth(); %#ok<NASGU>

%%% get parameter for annotation (para_A)
userdataAno = get(findobj('Tag','tab_annotation'),'UserData');
para_A = userdataAno.para_A; %#ok<NASGU>

myluts = getLuts(); %#ok<NASGU>
emopt  = getEmopt(false); %#ok<NASGU>
detopt = getDetectOpt(); %#ok<NASGU>
scaling       = str2num(get(findobj('Tag','edit_option_scaling'),'String')); %#ok<ST2NM,NASGU>
sigma_default = str2num(get(findobj('Tag','edit_option_sigma'),  'String')); %#ok<ST2NM,NASGU>
sizCache      = str2num(get(findobj('Tag','edit_option_cache'),  'String')); %#ok<ST2NM,NASGU>
% channel_detection = str2double(get(findobj('Tag','edit_option_detection'),'String')); %#ok<NASGU>
imname             = userdata.imname; %#ok<NASGU>
imext              = userdata.imext; %#ok<NASGU>
cum2               = userdata.cum2; %#ok<NASGU>
flag_virtual_stack = userdata.flag_virtual_stack; %#ok<NASGU>
rotinfo            = userdata.rotinfo; %#ok<NASGU>
flagAlign   = get(findobj('Tag','check_main_align'),'Value')==1; %#ok<NASGU>
flagRotateX = get(findobj('Tag','button_main_rotate_x'),'Value')==1; %#ok<NASGU>
flagRotateY = get(findobj('Tag','button_main_rotate_y'),'Value')==1; %#ok<NASGU>
flagRotateZ = get(findobj('Tag','button_main_rotate_z'),'Value')==1; %#ok<NASGU>
flagLn      = get(findobj('Tag','check_main_logarithm'),'Value')==1; %#ok<NASGU>
flagInterpZ = get(findobj('Tag','check_main_interpz'),'Value')==1; %#ok<NASGU>
flagSubBG   = get(findobj('Tag','check_main_subtract'),'Value')==1; %#ok<NASGU>

rotation    = userdata.mfivs.getRotation(); %#ok<NASGU>
flagMirror  = get(findobj('Tag','check_main_rotate_mirror'),'Value')==1; %#ok<NASGU>

save(pathSaveMat,'imname','emopt','myluts','scaling','sigma_default','cum2',...
    'flag_virtual_stack','rotinfo','idxColDisp','strStrokeColor','strSelectionColor',...
    'flagAlign','flagRotateX','flagRotateY','flagRotateZ','flagLn','flagInterpZ',...
    'flagSubBG','rotation','flagMirror','detopt','imext','cellColumnHeader',...
    'cellColumnWidth','para_A','sizCache');


%%% slow mode; for backward compatibility
%%% Because old roiedit programs has inadequate load algorithm,
%%% complete backward compatibility is impossible...
%%% Following codes are for tools other than roiedit.
if idx_filter==1
    hmtm = userdata.hmtm;
    idxRowFull = 0:hmtm.getRowCount()-1;
    
    % reconstruct data_edited with resort column order
    numcolor = userdata.himp.getNChannels();
    columnNames = {
        'Check','uniqueID','Name','frame','OoB','Repaired','min_dist','Name_min_dist',...
        'pi_k','xc','yc','zc','S_11','S_12','S_13','S_22','S_23','S_33',...
        'vol','min_scdist','Name_min_scdist'};
    for p=1:numcolor; columnNames(21+p) = {['pi_k_Ch',num2str(p)]}; end
    for p=1:numcolor; columnNames(21+numcolor+p) = {['pi_k_Ch',num2str(p),'_normres']}; end
    
    %     data_edited = cell(hmtm.getRowCount(),numel(columnNames));
    %     for p=1:numel(columnNames)
    %         tmpdata = hmtm.getValuesAt(idxRowFull,columnNames{p});
    %         if isa(tmpdata(1,1),'java.lang.String')
    %             data_edited(:,p) = ...
    %                 cellfun(@java.lang.String,cell(tmpdata),'UniformOutput',false);
    %         else
    %             data_edited(:,p) = num2cell(tmpdata);
    %         end
    %     end
    
    data_edited = cell(hmtm.getRowCount(),hmtm.getColumnCount());
    idxLast = numel(columnNames);
    for p=1:hmtm.getColumnCount()
        tmpColName = cell(hmtm.getColumnName(p-1));
        tmpIdx = find(strcmpi(tmpColName,columnNames));
        if isempty(tmpIdx)
            idxLast = idxLast+1;
            tmpIdx = idxLast;
            columnNames(tmpIdx) = tmpColName;
        end
        if ischar(hmtm.getValueAt(0,p-1))
            data_edited(:,tmpIdx) = ...
                cellfun(@java.lang.String,cell(hmtm.getValuesAt(idxRowFull,p-1)),...
                'UniformOutput',false);
        else
            data_edited(:,tmpIdx) = num2cell(hmtm.getValuesAt(idxRowFull,p-1));
        end
    end
    
    huicontrol_logarithm = findobj('Tag','logarithm');
    flag_logarithm = get(huicontrol_logarithm,'Value')==1; %#ok<NASGU>
    idx_frame     = hmtm.getValuesAt(idxRowFull,'frame')'; %#ok<NASGU>
    flag_oob      = hmtm.getValuesAt(idxRowFull,'OoB')'; %#ok<NASGU>
    flag_repaired = hmtm.getValuesAt(idxRowFull,'Repaired')'; %#ok<NASGU>
    params_em     = hmtm.getValuesAt(idxRowFull,'params_em')'; %#ok<NASGU>
    
    % save columnNames and data_edited with append mode
    save(pathSaveMat,'-append','columnNames','data_edited','params_em',...
        'flag_oob','flag_repaired','idx_frame','flag_logarithm');
end


%%% save table data
hmtm = userdata.hmtm;

% remove listeners before save
hl_tmp = hmtm.getTableModelListeners();
for p=1:numel(hl_tmp)
    hmtm.removeTableModelListener(hl_tmp(p));
end

% remove undo/redo stacks before save
hus = hmtm.undoStack.clone();
hrs = hmtm.redoStack.clone();
hmtm.undoStack.clear();
hmtm.redoStack.clear();

% remove cache before save
hMuCached = hmtm.muCached;
hSigmaCached = hmtm.sigmaCached;
hRotmat = hmtm.rotmat;
hmtm.muCached = [];
hmtm.sigmaCached = [];
hmtm.rotmat = [];

% save table data
hmtm.save(pathSaveDat);

% recover listeners
for p=1:numel(hl_tmp)
    hmtm.addTableModelListener(hl_tmp(p));
end

% recover undo/redo stacks
hmtm.undoStack = hus;
hmtm.redoStack = hrs;

% recover caches
hmtm.muCached = hMuCached;
hmtm.sigmaCached = hSigmaCached;
hmtm.rotmat = hRotmat;

disp(['saved to ',pathSaveMat]);

enableGUI(true);

end


function data_calc_callback(~,~)

t_total = tic;
userdata = getUserData();
hmtm = userdata.hmtm;
himp = userdata.himp;
emopt = getEmopt(true);
[~,optimfun] = data_calc_algorithm();
if emopt.ub_sigma>0 && emopt.selected>0
    optimfun = @optimze_constrained;
end
enableGUI(false);

%%% obtain target frame
frameCurrent = himp.getT();
diff_frame = -1;
hctrl_target = findobj('Tag','buttongroup_track_target');
switch hctrl_target.SelectedObject.String
    case '<<'
        diff_frame = 1;
        frameTarget = frameCurrent-1:-1:1;
    case '<'
        diff_frame = 1;
        frameTarget = frameCurrent-1;
    case 'Current'
        diff_frame = 0;
        frameTarget = frameCurrent;
    case '>'
        frameTarget = frameCurrent+1;
    case '>>'
        frameTarget = frameCurrent+1:himp.getNFrames();
    case 'Specify'
        frameTarget ...
            = str2num(get(findobj('Tag','edit_track_target_specify'),'String')); %#ok<ST2NM>
end
if ~all(isfinite(frameTarget))
    disp('Error: Target frame is not numeric.');
    enableGUI(true);
    return;
end


%%% obtain source frame
hctrl_origin = findobj('Tag','buttongroup_track_source');
switch hctrl_origin.SelectedObject.String
    case 'Default (Prev|Curr)'
        frameOrigin = frameTarget+diff_frame;
    case 'Similar'
        disp('Error: ''Similar'' does not work currently');
        enableGUI(true);
        return;
    case 'Specify'
        frameOrigin ...
            = str2num(get(findobj('Tag','edit_track_source_specify'),'String')); %#ok<ST2NM>
end
if ~all(isfinite(frameOrigin))
    disp('Error: Source frame is not numeric.');
    enableGUI(true);
    return;
end
if numel(frameTarget)~=numel(frameOrigin)
    disp('Error: Number of source frame is different from that of target frame.');
    enableGUI(true);
    return;
end

setupBlurredImage(); % apply settings of aligned image to blurred image

%%% optimize selected ROI only
numrow = hmtm.getRowCount();
flagSelected = false(numrow,1);
if emopt.selected==1
    flagSelected(getSelectedRowIndexInFigureAsModel()+1) = true;
    if ~any(flagSelected)
        disp('"Selected ROI only" option is valid when some ROIs are selected.');
        enableGUI(true);
        return;
    end
end


%%% main loop
for p=1:numel(frameTarget)
    t_each = tic;
    disp(['Optimize rois in frame #',num2str(frameTarget(p)),...
        ' based on frame #',num2str(frameOrigin(p)),' ...']);
    
    idxFull = 0:hmtm.getRowCount()-1;
    idxFrame = hmtm.getValuesAt(idxFull,'frame');
    flagTarget = idxFrame==frameTarget(p);
    flagOrigin = idxFrame==frameOrigin(p);
    
    %     im_target = double(getImageMiji(himp_blurred,[],1,frameTarget(p)));
    im_target = getBlurredImage([],frameTarget(p));
    im_target(im_target<emopt.thrint)=0;
    emopt.thrint = 0;
    
    %%% optimize selected ROI only
    if emopt.selected==1
        hmrm = userdata.hmrm;
        names = cell(hmrm.getSelectedRowNames());
        flagSelected = false(size(idxFrame));
        for q=1:numel(names)
            flagSelected = flagSelected | hmtm.strcmp(names{q},idxFull,'Name');
        end
        params_em_in_tmp = hmtm.getValuesAt(flagOrigin&~flagSelected,'params_em')';
        emopt_tmp = emopt;
        emopt_tmp.maxiter = 0;
        emopt_tmp.fixpik  = 1;
        emopt_tmp.fixmu   = 1;
        emopt_tmp.fixvol  = 1;
        [~,~,im_synth] = optimfun(params_em_in_tmp,im_target,emopt_tmp);
        im_target = im_target - im_synth;
        flagOrigin = flagOrigin & flagSelected;
        flagTarget = flagTarget & flagSelected;
    end
    
    %%% if no rois were found in source frame, return
    params_em_in = hmtm.getValuesAt(flagOrigin,'params_em')';
    if isempty(params_em_in)
        disp('ROIs not found in source frame.');
        enableGUI(true);
        return;
    end
    
    %%% do track
    params_em_out = optimfun(params_em_in,im_target,emopt);
    
    
    %%% do update
    huic = findobj('Tag','check_main_autosave');
    huic.Value = 0; % temporally disable autosave
    if frameTarget(p)==frameOrigin(p)
        % update
        hmtm.setValuesAt(params_em_out',flagTarget,'params_em',false);
    else
        % prepare data first
        numpos2 = size(params_em_out,2);
        maxID = max(hmtm.getValuesAt(~flagTarget,'uniqueID'));
        if isempty(maxID)
            maxID = 0;
        end
        uniqueID = (1:numpos2)' + maxID;
        
        data_tmp = userdata.defaultdata(ones(1,numpos2),:);
        data_tmp(:,cat(1,...
            hmtm.getColumnIndex('Check'),...
            hmtm.getColumnIndex('uniqueID'),...
            hmtm.getColumnIndex('Name'),...
            hmtm.getColumnIndex('frame'),...
            hmtm.getColumnIndex('params_em'),...
            hmtm.getColumnIndex('Landmarks') )+1 ...
            ) = cat(2,...
            num2cell(hmtm.getValuesAt(flagOrigin,'Check')),... % Check
            num2cell(uniqueID),... % uniqueID
            cell(hmtm.getValuesAt(flagOrigin,'Name')),... % Name
            num2cell(frameTarget(p)*ones(numpos2,1)),... % frame
            num2cell(params_em_out'),... % params_em
            num2cell(hmtm.getValuesAt(flagOrigin,'Landmarks'))); % landmarks
        
        % remove old data
        hmtm.removeRowsAt(flagTarget,[false,true]); % unfire and add undo
        
        % insert new data
        idx_insert = find(...
            hmtm.getValuesAt((1:hmtm.getRowCount())-1,'frame')<frameTarget(p),1,'last');
        if isempty(idx_insert)
            idx_insert = 0;
        end
        hmtm.insertRowsAt(data_tmp,idx_insert,[false,true]); %  unfire and add undo
        hmtm.convertCharToString();
    end
    
    if get(findobj('Tag','check_track_option_main'),'Value')~=1
        data_calc_auxiliarly(frameTarget(p));
    end
    
    
    userdata.hmrm.invalidateAll(); % force update roi display
    displayCurrentSelection();
    
    fprintf('Optimize #%d->#%d: finished in %g (sec) \n',...
        frameOrigin(p),frameTarget(p),toc(t_each));
    
end

fprintf('All tracking process finished in %g (sec) \n\n',toc(t_total));
enableGUI(true);

end


function [strlist,fun] = data_calc_algorithm()
%%% this function manages tracking algorithm

% list of algorithm
algorithm = {
    'em_22',@elipsoid_em_22;
    'em_13',@elipsoid_em_13;
    'em_12',@elipsoid_em_12;
    'em_9', @elipsoid_em_9;
    'em_8', @elipsoid_em_8;
    'em_7', @elipsoid_em_7;
    };

% return strlist for GUI
strlist = algorithm(:,1);

% return function for tracking
if nargout==2
    tmpidx = get(findobj('Tag','popup_track_option_algorithm'),'value');
    fun = algorithm{tmpidx,2};
end

end


function data_calc_auxiliarly(frame_target)
% This function calculate values for auxiliarly columns, in silent

% setup
userdata = getUserData();
hmtm = userdata.hmtm;
idx_frame = hmtm.getValuesAt((1:hmtm.getRowCount())-1,'frame');
flag_target = idx_frame==frame_target;
names = cell(hmtm.getValuesAt(flag_target,'Name'));
params_em = hmtm.getValuesAt(flag_target,'params_em')';
flag_silent = [false,false,false];
scaling = str2num(get(findobj('Tag','edit_option_scaling'),'String')); %#ok<ST2NM>

% minimum scaled distance and the partner's name
numpos2 = sum(flag_target);
scaled_distance = calc_scaled_distance(params_em);
scaled_distance(find(eye(numpos2))) = inf; %#ok<FNDSB>
[min_scdist_val,min_scdist_idx] = min(scaled_distance,[],2);
hmtm.setValuesAt(min_scdist_val,flag_target,'min_scdist',flag_silent);
hmtm.setValuesAt(names(min_scdist_idx),flag_target,'Name_min_scdist',flag_silent);

% volume
vol = zeros(size(params_em,2),1);
for p=1:size(params_em,2)
    vol(p) = det(params_em([5,6,7;6,8,9;7,9,10]+(p-1)*10));
end
hmtm.setValuesAt(vol,flag_target,'vol',flag_silent);

% minimum distance and the partner's name
pos_fit = params_em(2:4,:);
disttable_fit = squareform(pdist(bsxfun(@times,pos_fit,scaling')'));
disttable_fit(find(eye(numpos2))) = inf; %#ok<FNDSB>
[mindistval,mindistidx] = min(disttable_fit,[],2);
hmtm.setValuesAt(mindistval,flag_target,'min_dist',flag_silent);
hmtm.setValuesAt(names(mindistidx),flag_target,'Name_min_dist',flag_silent);

% obtain fitted intensity for each channel of aligned image
disp('obtain intensities in each channel ...');
emopt_tmp = getEmopt(true);
emopt_tmp.maxiter = 1;
emopt_tmp.fixpik  = 0;
emopt_tmp.fixmu   = 1;
emopt_tmp.fixvol  = 1;
emopt_tmp.thrint  = -inf;
emopt_tmp.lambda_pi = 0;

[~,optimfun] = data_calc_algorithm();
numcolor = userdata.himp.getNChannels();
pi_k = hmtm.getValuesAt(flag_target,'pi_k');
pi_k_ch = zeros(numpos2,numcolor);
pi_k_ch_normres = zeros(numpos2,numcolor);
for p=1:numcolor
    im_target = double(getImageMiji(userdata.himp,[],p,frame_target));
    params_em_out = optimfun(params_em,im_target,emopt_tmp);
    pi_k_ch(:,p) = params_em_out(1,:);
    mdl = fitlm(pi_k,pi_k_ch(:,p),'RobustOpts','on');
    pi_k_ch_normres(:,p) = mdl.Residuals.Standardized;
end
hmtm.setValuesAt(cat(2,pi_k_ch,pi_k_ch_normres),flag_target,'pi_k_Ch_all',flag_silent);

end



function [params_em_out,ofv,im_synth] = ...
    optimze_constrained(params_em_in,im_target,emopt)
%%% Called instead of em_methods when the upper bound of sigma was specified
%%% Have same inputs and outputs as em_methods
%%% [params_em_out,ofv,im_synth,benchmark] = em_methods(params_em_in,im_target,emopt);
%%% this function handles following fileds of emopt:
%%%    maxiter, tol, fixvol, fixmu, fixpik, ub_sigma
[~,optimfun] = data_calc_algorithm();
tmpemopt = emopt;
tmpemopt.maxiter=0;

optimopt = optimset('maxiter',emopt.maxiter,'TolFun',emopt.tol,'Display','iter');
if emopt.maxiter==1
    optimopt = optimset(optimopt,'Display','off');
end
numrois = size(params_em_in,2);
lb_orig = zeros(10,1);
ub_orig = inf(4,1);
ub_orig(5:10) = emopt.ub_sigma;

% handle fixed parameters
flagFixed = false(10,1);
if emopt.fixpik==1
    flagFixed(1) = true;
end
if emopt.fixmu==1
    flagFixed(2:4) = true;
end
if emopt.fixvol==1
    flagFixed(5:10) = true;
end

x0 = reshape(params_em_in(~flagFixed,:),sum(~flagFixed)*numrois,1);
lb = lb_orig(~flagFixed);
ub = ub_orig(~flagFixed);

params_em_out = params_em_in;
ofv = inf;
if any(~flagFixed)
    [x,ofv] = fmincon(@(x)getSecondRet(optimfun,x,im_target,tmpemopt,params_em_in,flagFixed),...
        x0,[],[],[],[],lb,ub,[],optimopt);
    params_em_out(~flagFixed,:) = reshape(x,sum(~flagFixed),numrois);
end

if nargout==3
    [~,~,im_synth] = optimfun(params_em_out,im_target,tmpemopt);
end

end


function ret = getSecondRet(hfun,x,im_target,tmpemopt,params_em,flagFixed)
numrois = size(params_em,2);
params_em_in = params_em;
params_em_in(~flagFixed,:) = reshape(x,sum(~flagFixed),numrois);
[~,ret] = hfun(params_em_in,im_target,tmpemopt);
end



%%% table callback functions


function table_createGUI()

userdata = getUserData();

%%% make TableModel and override getColumnClass function
%%% myTableModel.class should be placed in current dir or javaclasspath
hmtm = javaObjectEDT(userdata.hmtm);

%%% create SortableTable of JIDE from TableModel and enable sorting & filtering
hqtff = javaObjectEDT(com.jidesoft.grid.QuickTableFilterField(hmtm));
hftm = javaObjectEDT(com.jidesoft.grid.FilterableTableModel(hqtff.getDisplayTableModel()));
hstm = javaObjectEDT(com.jidesoft.grid.SortableTableModel(hftm));
jtable = javaObjectEDT(com.jidesoft.grid.SortableTable(hstm));
tableHeader = javaObjectEDT(com.jidesoft.grid.AutoFilterTableHeader(jtable));

% hstm = javaObjectEDT(com.jidesoft.grid.SortableTableModel(hmtm));
% jtable = javaObjectEDT(com.jidesoft.grid.SortableTable(hstm));
% tableHeader = javaObjectEDT(com.jidesoft.grid.AutoFilterTableHeader(jtable));


tableHeader.setAutoFilterEnabled(true);
tableHeader.setShowFilterName(true);
tableHeader.setShowFilterIcon(true);
jtable.setTableHeader(tableHeader);

% % create rowsorter for narrowing
% htrs = javax.swing.table.TableRowSorter(jtable.getModel());
% jtable.setRowSorter(htrs);

% hqtff.setColumnIndices(hmtm.getColumnIndex('Name_popup'));
hqtff.setSearchingColumnIndices(hmtm.getColumnIndex('Name_popup'));

%%% display boolean check box
idxlogical = find(cellfun(@islogical,userdata.defaultdata));
for p=1:numel(idxlogical)
    hcolumn = jtable.getColumnModel.getColumn(idxlogical(p)-1);
    hcolumn.setCellEditor(com.jidesoft.grid.BooleanCheckBoxCellEditor);
    hcolumn.setCellRenderer(com.jidesoft.grid.BooleanCheckBoxCellRenderer);
end

%%% In order to enable sliders, disable AutoResize
jtable.setAutoResizeMode(jtable.AUTO_RESIZE_OFF);

%%% If clicked twice, enable editing
jtable.setClickCountToStart(2);

%%% display
frame = javax.swing.JFrame('Customized ROI Manager');
frame.add(javax.swing.JScrollPane(jtable));
frame.pack();
frame.setVisible(true);
frame.setFocusable(true);

%%% enabling convenient features
installer = com.jidesoft.grid.TableHeaderPopupMenuInstaller(jtable);
pmCustomizer1=com.jidesoft.grid.AutoResizePopupMenuCustomizer;
installer.addTableHeaderPopupMenuCustomizer(pmCustomizer1);
pmCustomizer2=com.jidesoft.grid.TableColumnChooserPopupMenuCustomizer;
installer.addTableHeaderPopupMenuCustomizer(pmCustomizer2);

%%% In order to hold listneres, hold handles in matlab variable
hlistener.jtable(1) = addlistener(jtable,'MouseReleased',@table_mouse_callback);
hlistener.jtable(2) = addlistener(jtable,'KeyReleased',@table_key_callback);
hlistener.myTableModel = addlistener(hmtm,'tableChanged',@table_data_callback);
userdata.htable_callback_semaphore = java.util.concurrent.Semaphore(1,true);
userdata.htable_changed_semaphore = java.util.concurrent.Semaphore(1,true);

%%% hold handles
userdata.hmtm   = hmtm;
userdata.hqtff  = hqtff;
userdata.hftm   = hftm;
userdata.hstm   = hstm;
userdata.jtable = jtable;
userdata.frame  = frame;
userdata.hlistener = hlistener;

%%% sorting by frame number -> SortableTableModel.setMasterSortColumns
hstm.setMasterSortColumns(hmtm.getColumnIndex('frame'));

%%% enable or disable auto resort (enable is the default setting)
jtable.setAutoResort(false);


setUserData(userdata);

%%% set column width and header name
if isfield(userdata,'cellColumnHeader')
    renameColumnHeader(userdata.cellColumnHeader);
end
if isfield(userdata,'cellColumnWidth')
    setColumnWidth(userdata.cellColumnWidth);
end

userdata = getUserData();


%%% update column order
if isfield(userdata,'idxColDisp')
    idxColDisp = userdata.idxColDisp;
    difference = hmtm.getColumnCount() - numel(idxColDisp);
    idxColDisp = cat(2,idxColDisp,(1:difference) + max(idxColDisp));
else
    idxColDisp = 0:hmtm.getColumnCount()-1;
end

[colDispSorted,sortedidx] = sort(idxColDisp);
[~,~,idxColDispSorted] = unique(colDispSorted); % avoiding non-continous colDispSorted
if colDispSorted(1)==-1
    idxColDispSorted = idxColDispSorted-2;
else
    idxColDispSorted = idxColDispSorted-1;
end
hcm = javaObjectEDT(jtable.getColumnModel);
hci = userdata.hmtm.getColumnIdentifiers();
columnNames = cell(hci.toArray());
for p=1:numel(colDispSorted)
    if colDispSorted(p)==-1 % remove column
        hcm.removeColumn(jtable.getColumn(columnNames{sortedidx(p)}));
    else % move column
        %         hcm.moveColumn(hcm.getColumnIndex(columnNames{sortedidx(p)}),...
        %             colDispSorted(p));
        hcm.moveColumn(hcm.getColumnIndex(columnNames{sortedidx(p)}),...
            idxColDispSorted(p));
        
    end
end

setUserData(userdata);
%
% %%% set column width and header name
% if isfield(userdata,'cellColumnHeader')
%     renameColumnHeader(userdata.cellColumnHeader);
% end
% if isfield(userdata,'cellColumnWidth')
%     setColumnWidth(userdata.cellColumnWidth);
% end

end


function table_mouse_callback(hObj,ev)
userdata = getUserData();

if userdata.htable_callback_semaphore.availablePermits() == 0
    return;
end
userdata.htable_callback_semaphore.acquire();

if ev.getButton()==ev.BUTTON1 % left click
    table_selection_update();
elseif ev.getButton()==ev.BUTTON2 || ev.getButton()==ev.BUTTON3 % right click
    show_popup_menu_callback(hObj,ev);
end

if userdata.htable_callback_semaphore.availablePermits() == 0
    userdata.htable_callback_semaphore.release();
end

end


function table_key_callback(~,ev)
userdata = getUserData();

if userdata.htable_callback_semaphore.availablePermits() == 0
    return;
end
userdata.htable_callback_semaphore.acquire();

if         ev.getKeyCode() == ev.VK_UP...
        || ev.getKeyCode() == ev.VK_DOWN...
        || ev.getKeyCode() == ev.VK_PAGE_UP...
        || ev.getKeyCode() == ev.VK_PAGE_DOWN...
        || ev.getKeyCode() == ev.VK_HOME...
        || ev.getKeyCode() == ev.VK_END...
        || ev.getKeyCode() == ev.VK_SPACE
    % || ev.getKeyCode() == ev.VK_CONTROL...
    % || ev.getKeyCode() == ev.VK_SHIFT
    
    table_selection_update();
    
elseif ev.getKeyCode() == ev.VK_DELETE % delete table rows
    removeROIs();
    
elseif bitand(ev.getModifiersEx(),ev.CTRL_DOWN_MASK)==ev.CTRL_DOWN_MASK % ctrl pressed
    if ev.getKeyCode() == ev.VK_A % select all
        table_selection_update();
        
    elseif ev.getKeyCode() == ev.VK_Z % undo
        undo();
        
    elseif ev.getKeyCode() == ev.VK_Y % redo
        redo();
    end
    
end

if userdata.htable_callback_semaphore.availablePermits() == 0
    userdata.htable_callback_semaphore.release();
end

end


function table_selection_update()
%%% Modify selection status for MyRoiManager.
%%% Assuming selection status were changed through the table GUI.
%%% Selection status for the cells which is visible in the table sholud be
%%% matched whether the cells are selected in the table.
%%% Selection status for the cells which is invisible in the table should be held.

userdata = getUserData();
hmtm = userdata.hmtm;
hmrm = userdata.hmrm;

idxRowsSelectedInTable  = getSelectedRowIndexInTableAsModel();
flag_valid = idxRowsSelectedInTable~=-1;
idxRowsInTable = idxRowsSelectedInTable(flag_valid);
idxRowsSelectedInFigure = getSelectedRowIndexInFigureAsModel();
flag_visible = convertRowIndexToView(idxRowsSelectedInFigure)>=0;
idxRowsInFigure = idxRowsSelectedInFigure(flag_visible);
idxRowsToBeDeselected = setdiff(idxRowsInFigure,idxRowsInTable);
idxRowsToBeSelected   = setdiff(idxRowsInTable, idxRowsInFigure);

if ~isempty(idxRowsToBeDeselected)
    names  = hmtm.getValuesAt(idxRowsToBeDeselected,'Name');
    frames = hmtm.getValuesAt(idxRowsToBeDeselected,'frame');
    hmrm.removeSelection(frames,names(:,1));
end

if ~isempty(idxRowsToBeSelected)
    names  = hmtm.getValuesAt(idxRowsToBeSelected,'Name');
    frames = hmtm.getValuesAt(idxRowsToBeSelected,'frame');
    hmrm.addSelection(frames,names(:,1),idxRowsToBeSelected);
end

displayCurrentSelection();

end


function table_switch_roiname(input1,input2)
%%% change or switch roiname
%%% If name_new does not exist in the table yet, name_old was changed to
%%% name_new. If not, name_old and name_new were switched. These two
%%% functionality can be implemented as single code using boolean mask.
%%%
%%% If called from table_data_callback or show_popup_menu_callback, input1
%%% and input2 should be scalar row indices in matlab style (1-based index),
%%% indicating rows switching their names.
%%% If called from approveEstimatedName, input1 is a scalar row index in
%%% Java style (0-based index) and input2 is new name as an char array.
%%%

userdata = getUserData();
hmtm = userdata.hmtm;
hmrm = userdata.hmrm;

%%% branch based on caller type (direct callback or indirect)
if isnumeric(input2) % direct callback
    row_old = input1;
    row_new = input2;
    if row_old == row_new %%% roi name was changed in the table
        % Because roi name was already changed, corrupt the change at first.
        name_old = char(hmtm.oldobj);
        name_new = char(hmtm.getValuesAt(row_new-1,'Name'));
        hmtm.undoStack.pop(); % remove the most recent element from undoStack
        if strcmp(name_old,name_new); return; end
        hmtm.setValuesAt({name_old},row_old-1,'Name',[false,false,false]);
        semaphore = userdata.htable_changed_semaphore;
    else %%% two rois were specified in the canvas
        name_old = char(hmtm.getValuesAt(row_old-1,'Name'));
        name_new = char(hmtm.getValuesAt(row_new-1,'Name'));
        semaphore = userdata.hcanvas_callback_semaphore;
    end
else % indirect callback
    row_new = input1+1;
    name_old = char(hmtm.getValuesAt(input1,'Name'));
    name_new = input2;
    semaphore = userdata.htable_changed_semaphore;
end

%%% detect overlaps
idx_full = (1:hmtm.getRowCount())-1;
frame = hmtm.getValuesAt(idx_full,'frame');
flag_frame = (frame>=frame(row_new));
flag_overlap_old           = hmtm.strcmp(name_old,idx_full,'Name')            & flag_frame;
flag_overlap_new           = hmtm.strcmp(name_new,idx_full,'Name')            & flag_frame;
flag_overlap_old_mindist   = hmtm.strcmp(name_old,idx_full,'Name_min_dist')   & flag_frame;
flag_overlap_new_mindist   = hmtm.strcmp(name_new,idx_full,'Name_min_dist')   & flag_frame;
flag_overlap_old_minscdist = hmtm.strcmp(name_old,idx_full,'Name_min_scdist') & flag_frame;
flag_overlap_new_minscdist = hmtm.strcmp(name_new,idx_full,'Name_min_scdist') & flag_frame;


%%% Confirmation dialog:
%%% If name_new is already exist in the table, confirmation dialog for
%%% switching name_new and name_old will be displayed.
%%%
%%% Temporally release the semaphore because 'drawnow' in
%%% questdlg permits interruption and deadlock may occur
flag_switch = any(flag_overlap_new);
if flag_switch
    if semaphore.availablePermits() == 0
        semaphore.release();
    end
    answer = questdlg(...
        ['Do you want to switch ROI names with ',name_old,...
        ' and ',name_new,' in the current & subsequent frames?'],...
        'switch ROI names','Yes','No','No');
    semaphore.acquire();
    
    %%% if cancel, display warnings, and return (do nothing)
    if ~strcmpi(answer,'Yes')
        if semaphore.availablePermits() == 0
            semaphore.release();
        end
        enableGUI(true);
        warndlg(['Changing to new name is cancelled because the name ''',...
            name_new,''' is already exist in frame #',...
            num2str(frame(find(flag_overlap_new,1))),' or after.']);
        return;
    end
end

%%% update table
if semaphore.availablePermits() == 0
    semaphore.release();
end

flag_change_all = flag_overlap_old|flag_overlap_new ...
    | flag_overlap_old_mindist | flag_overlap_new_mindist ...
    | flag_overlap_old_minscdist | flag_overlap_new_minscdist ;
hmtm.setValuesAt(frame(flag_change_all),flag_change_all,'frame',[false,true,false]); % dummy for undo
tmpfun = @(name,flag,colname) ...
    hmtm.setValuesAt(repmat({name},[sum(flag),1]),flag,colname,false(1,3));
tmpfun(name_new,flag_overlap_old,'Name');
tmpfun(name_old,flag_overlap_new,'Name');
tmpfun(name_new,flag_overlap_old_mindist,'Name_min_dist');
tmpfun(name_old,flag_overlap_new_mindist,'Name_min_dist');
tmpfun(name_new,flag_overlap_old_minscdist,'Name_min_scdist');
tmpfun(name_old,flag_overlap_new_minscdist,'Name_min_scdist');

% update name_log
if any(flag_overlap_old)
    name_log_old = strcat(name_old,',',cell(hmtm.getValuesAt(flag_overlap_old,'Name_log')));
    hmtm.setValuesAt(name_log_old,flag_overlap_old,'Name_log',false(1,3));
end
if any(flag_overlap_new)
    name_log_new = strcat(name_new,',',cell(hmtm.getValuesAt(flag_overlap_new,'Name_log')));
    hmtm.setValuesAt(name_log_new,flag_overlap_new,'Name_log',false(1,3));
end

hmtm.convertCharToString(); % call directly
hmrm.invalidateAll(); % force update roi display % call directly
autosave(); % call directly
clearSelection();


%%% set selection
jtable = userdata.jtable;
% row_at = com.jidesoft.grid.TableModelWrapperUtils.getRowAt(...
%     jtable.getModel,row_new-1);
% if row_at>=0
%     jtable.addRowSelectionInterval(row_at,row_at);
% end
row_at = convertRowIndexToView(row_new-1);
if row_at>=0
    jtable.scrollRowToVisible(row_at);
end
jtable.repaint();
if flag_switch % switch two rois
    disp(['ROI names were switched with ',name_old,' and ',name_new,...
        ' in frame #',num2str(frame(find(flag_overlap_old,1))),' or after.']);
else % change rois
    disp(['ROI names were changed from ',name_old,' to ',name_new,...
        ' in frame #',num2str(frame(find(flag_overlap_old,1))),' or after.']);
end

end


function table_data_callback(~,ev)
%%% callback for changing roi data only thorugh the table GUI.

userdata = getUserData();
hmtm = userdata.hmtm;
hmrm = userdata.hmrm;

if exist('ev','var') ...
        && ev.getType() == ev.UPDATE ...
        && ev.getColumn() == hmtm.getColumnIndex('Name') ...
        && ev.getFirstRow()==ev.getLastRow()
    % name was changed
    row = ev.getFirstRow()+1;
    table_switch_roiname(row,row);
    return;
    % elseif exist('ev','var') ...
    %         && ev.getType() == ev.UPDATE ...
    %         && ev.getColumn() == hmtm.getColumnIndex('Name_anno_human') ...
    %         && ev.getFirstRow()==ev.getLastRow()
    %     % name_anno_human was changed
    %     row = ev.getFirstRow()+1;
    %     edit_name_anno_human(row);
    %     return;
    %
end

hmrm.invalidateAll(); % force update roi display
autosave();
table_selection_update();
end


%%% Canvas Callback


function canvas_mouse_callback(hObj,ev)
%%% Left Click With Ctrl: generate and register ROI
% %%% Left Click With Ctrl+Shift: remove ROI
%%% Left Click With Ctrl+Alt: Move ROI position
%%% Left Click With Shift: Multi-Select ROI
%%% Left Click With Shift+Alt: deselect ROI

userdata = getUserData();

if userdata.hcanvas_callback_semaphore.availablePermits() == 0
    return;
end
userdata.hcanvas_callback_semaphore.acquire();

if ev.getButton == ev.BUTTON3 ... % right click
        || ev.getButton == ev.BUTTON2 % middle click
    show_popup_menu_callback(hObj,ev); % show popup
    
else % assuming left click
    
    modifier = ev.getModifiersEx();
    mask_alt   = ev.ALT_DOWN_MASK;
    mask_ctrl  = ev.CTRL_DOWN_MASK;
    mask_shift = ev.SHIFT_DOWN_MASK;
    flag_alt   = bitand(modifier,mask_alt  )==mask_alt;
    flag_ctrl  = bitand(modifier,mask_ctrl )==mask_ctrl;
    flag_shift = bitand(modifier,mask_shift)==mask_shift;
    
    if flag_ctrl && ~flag_shift &&  flag_alt % Ctrl+Alt: Move ROI position
        canvas_mouse_move(hObj,ev);
        
    elseif flag_ctrl && ~flag_shift && ~flag_alt % Ctrl: generate and register ROI
        canvas_mouse_generate(hObj,ev);
        
    elseif ~flag_ctrl && flag_shift &&  flag_alt % Shift+Alt: deselect ROI
        clearSelection();
        
    elseif flag_shift && ~flag_alt % Multi-Select
        canvas_mouse_select(hObj,ev);
        
    end
end

if userdata.hcanvas_callback_semaphore.availablePermits() == 0
    userdata.hcanvas_callback_semaphore.release();
end

end


function canvas_mouse_move(hObj,ev)
%%%Move ROI position
userdata = getUserData();
hmtm = userdata.hmtm;
pos = canvas_mouse_setup(hObj,ev);
% modelRows = getSelectedRowIndexInTableAsModel();
modelRows = getSelectedRowIndexInFigureAsModel();

if numel(modelRows)==1
    name_select = char(hmtm.getValuesAt(modelRows,'Name'));
    oldpos = hmtm.getValuesAt(modelRows,'mu');
    hmtm.setValuesAt(pos',modelRows,'mu',[true,true]);
    disp(['ROI (',name_select,') was moved from (',...
        num2str(oldpos(:)'),') to (',num2str(pos(:)'),').']);
end

end


function canvas_mouse_generate(hObj,ev)
%%% generate and register ROI
userdata = getUserData();
hmtm = userdata.hmtm;
jtable = userdata.jtable;
[pos,flag_near,~,frame_target,idx_insert] = canvas_mouse_setup(hObj,ev);

flag_add_ROI = true;
if flag_near
    flag_add_ROI = false;
    answer = questdlg('Do you add ROI near the existing ROI(s)?',...
        'adding ROI','Yes','No','No');
    drawnow();
    if strcmpi(answer,'Yes')
        flag_add_ROI = true;
    end
end

if flag_add_ROI
    uid = max(hmtm.getValuesAt((1:hmtm.getRowCount())-1,'uniqueID'))+1;
    if isempty(uid); uid=1; end
    name_uid = num2str(uid,'%-d');
    sigma_default = str2num(get(findobj('Tag','edit_option_sigma'),  'String')); %#ok<ST2NM>
    sigma = zeros(1,6);
    sigma([1,4,6]) = sigma_default;
    tabledata = userdata.defaultdata;
    tabledata(hmtm.getColumnIndex('uniqueID')+1) = {uid};
    tabledata(hmtm.getColumnIndex('Name')+1)     = {name_uid};
    tabledata(hmtm.getColumnIndex('frame')+1)    = {frame_target};
    tabledata(hmtm.getColumnIndex('mu')+1)       = num2cell(pos);
    tabledata(hmtm.getColumnIndex('sigma')+1)    = num2cell(sigma);
    hmtm.insertRowsAt(tabledata,idx_insert,[true,true]);
    disp(['Point ( ',num2str(pos'),' ) was added as a ROI (',name_uid,')']);
    
    %     row_at = com.jidesoft.grid.TableModelWrapperUtils.getRowAt(...
    %         jtable.getModel(),idx_insert);
    row_at = convertRowIndexToView(idx_insert);
    if row_at>=0
        jtable.addRowSelectionInterval(row_at,row_at);
    end
    jtable.scrollRowToVisible(row_at);
    jtable.repaint();
end

end


function canvas_mouse_select(hObj,ev)
%%% Multi-Select ROI
userdata = getUserData();
hmtm = userdata.hmtm;
hmrm = userdata.hmrm;
jtable = userdata.jtable;
[~,flag_near,minidx] = canvas_mouse_setup(hObj,ev);

if flag_near
    %     row_at = com.jidesoft.grid.TableModelWrapperUtils.getRowAt(...
    %         jtable.getModel,minidx-1);
    row_at = convertRowIndexToView(minidx-1);
    frame_select = hmtm.getValuesAt(minidx-1,'frame');
    name_select  = char(hmtm.getValuesAt(minidx-1,'Name'));
    if hmrm.isSelect(frame_select,name_select) % remove ROI from the selection list
        if row_at>=0
            jtable.removeRowSelectionInterval(row_at,row_at);
        end
        hmrm.removeSelection(frame_select,name_select);
        disp(['ROI (',name_select,') was removed from selection.']);
    else % add ROI to the selection list
        if row_at>=0
            jtable.addRowSelectionInterval(row_at,row_at);
            jtable.scrollRowToVisible(row_at);
        end
        jtable.repaint();
        hmrm.addSelection(frame_select,name_select,minidx-1);
        disp(['ROI (',name_select,') was added to selection.']);
    end
end
displayCurrentSelection();
end


function [pos,flag_near,minidx,frame_target,idx_insert] = canvas_mouse_setup(hObj,ev)

%%% helper function for canvas_mouse_callback

userdata = getUserData();
hmtm = userdata.hmtm;
himp = userdata.himp;

idx_frame = hmtm.getValuesAt((1:hmtm.getRowCount())-1,'frame')';
frame_target = himp.getT();
flag_target = idx_frame==frame_target;

diameter_check = str2double(get(findobj('Tag','edit_option_check'),'String'));
scaling = str2num(get(findobj('Tag','edit_option_scaling'),'String')); %#ok<ST2NM>
if get(findobj('Tag','check_main_interpz'),'Value')==1; scaling(3) = scaling(2); end
z_ratio = scaling(3)/scaling(2);

% cursorLoc = get(hObj,'cursorLoc');
% cursorLoc = ev.getPoint();
hc = ev.getComponent();
hp = ev.getPoint();
clx = hc.offScreenX(hp.getX);
cly = hc.offScreenY(hp.getY);

if himp.getNSlices()==1 % no z-slice
    pos = [clx+1; cly+1; 1];
else % z-sliced
    hmov = javaObjectEDT(MyOrthogonal_Views.getInstance());
    if isempty(hmov)
        pos = [clx+1; cly+1; himp.getZ()];
    else
        himp_tmp = hObj.getImage();
        if himp_tmp == hmov.getImage()
            pos = [clx+1; cly+1; himp.getZ()];
        elseif himp_tmp == hmov.getXZImage()
            crossLoc = double(hmov.getCrossLoc());
            pos = [clx+1; crossLoc(2); round(cly/z_ratio)+1];
        elseif himp_tmp == hmov.getYZImage()
            crossLoc = double(hmov.getCrossLoc());
            pos = [crossLoc(1); cly+1; round(clx/z_ratio)+1];
        end
    end
end

if sum(flag_target)>0
    pos_orig = hmtm.getValuesAt(flag_target,'mu')';
    numpos = size(pos_orig,2);
    disttable = sum((  (pos_orig - pos(:,ones(1,numpos))) ...
        .* scaling(ones(1,numpos),:)'  ).^2);
    [mindist,minidx_in_targets] = min(disttable);
    idx_target = find(flag_target);
    idx_insert = find(flag_target,1,'last');
    minidx = idx_target(minidx_in_targets);
else %%
    minidx = inf;
    mindist = inf;
    idx_insert = find(idx_frame==himp.getT()-1,1,'last');
    if isempty(idx_insert)
        idx_insert = 0;
    end
end

flag_near = mindist<=diameter_check^2;

end


function canvas_key_callback(~,ev)

userdata = getUserData();
if userdata.hcanvas_callback_semaphore.availablePermits() == 0
    return;
end
userdata.hcanvas_callback_semaphore.acquire();

modifier = ev.getModifiersEx();
mask_ctrl  = ev.CTRL_DOWN_MASK;
mask_shift = ev.SHIFT_DOWN_MASK;
flag_ctrl  = bitand(modifier,mask_ctrl )==mask_ctrl;
flag_shift = bitand(modifier,mask_shift)==mask_shift;

if ev.getKeyCode() == ev.VK_SPACE % toggle ROI names or optimization
    if ~flag_ctrl && ~flag_shift % toggle display roi
        huicontrol_roi = findobj('Tag','check_main_showroi');
        if get(huicontrol_roi,'Value')==1
            set(huicontrol_roi,'Value',false);
        else
            set(huicontrol_roi,'Value',true);
        end
        checkbox_showRoi_callback();
        
    elseif  flag_ctrl && ~flag_shift % optimize current frame without warnings
        set(findobj('Tag','buttongroup_track_target'),'SelectedObject',...
            findobj('Tag','button_track_target_current'));
        set(findobj('Tag','buttongroup_track_source'),'SelectedObject',...
            findobj('Tag','button_track_source_default'));
        data_calc_callback();
        
    elseif ~flag_ctrl &&  flag_shift % toggle display name
        huicontrol_roi = findobj('Tag','check_main_showroi');
        if get(huicontrol_roi,'Value')==1
            huicontrol_name = findobj('Tag','check_main_showname');
            if get(huicontrol_name,'Value')==1
                set(huicontrol_name,'Value',false);
            else
                set(huicontrol_name,'Value',true);
            end
            checkbox_showName_callback();
        end
    end
    
elseif ev.getKeyCode() == ev.VK_Z  && flag_ctrl % undo
    undo();
    
elseif ev.getKeyCode() == ev.VK_Y  && flag_ctrl % redo
    redo();
    
elseif ev.getKeyCode() == ev.VK_DELETE % remove selected roi
    removeROIs();
end

if userdata.hcanvas_callback_semaphore.availablePermits() == 0
    userdata.hcanvas_callback_semaphore.release();
end

end


%%% callback functions for other GUIs

function checkbox_interpz_callback(~,~)
% obtain parameters
userdata = getUserData();
hcanvas = userdata.hcanvas;
himp = userdata.himp;
mfivs = userdata.mfivs;
hmtm = userdata.hmtm;

% pre-processing
disp('Preparing for Z interpolation...');
ht = tic;
magnification = hcanvas.getMagnification();
origstate = MyOrthogonal_Views.isOrthoViewsImage(himp);
enableGUI(false);

% % calc zscale
hobj_scaling = findobj('Tag','edit_option_scaling');

% calc interpolated parameters
flag_interpz = get(findobj('Tag','check_main_interpz'),'Value')==1;
if flag_interpz % enabling interpolation
    set(hobj_scaling,'Enable','inactive'); % disable input for scaling
else % disabling interpolation
    set(hobj_scaling,'Enable','on'); % enable input for scaling
end

% do update
if MyOrthogonal_Views.isOrthoViewsImage(himp)
    hmov = javaObjectEDT(MyOrthogonal_Views());
    hmov.run(''); % force dispose
end
mfivs.setInterpZ(flag_interpz);
hmtm.setInterpZ(flag_interpz);
userdata.frame.repaint();
userdata.hmrm.invalidateAll(); % force update roi display


% post-processing
hcanvas.setMagnification(magnification);
if origstate ~= MyOrthogonal_Views.isOrthoViewsImage(himp)
    ij.WindowManager.setCurrentWindow(himp.getWindow());
    hmov = javaObjectEDT(MyOrthogonal_Views());
    hmov.run('');
end
addCanvasListener();
checkbox_showRoi_callback();

% update image
if himp.isComposite()
    tmpc = himp.getC();
    tmpluts = himp.getLuts();
    himp.setChannelsUpdated();
    himp.setLuts(tmpluts);
    himp.updateChannelAndDraw();
    himp.setC(tmpc);
end

enableGUI(true);

disp(['finished in ',num2str(toc(ht)),'(sec)']);

end


function checkbox_logarithm_callback(~,~)
userdata = getUserData();
flag_log = get(findobj('Tag','check_main_logarithm'),'Value')==1;
userdata.mfivs.setLn(flag_log);
end


function checkbox_proj_callback(~,~)

userdata = getUserData();
flag_showProj = get(findobj('Tag','check_main_proj'),'Value')==1;
himp = userdata.himp;

origstate = MyOrthogonal_Views.isOrthoViewsImage(himp);

enableGUI(false);

if flag_showProj
    fprintf('calculating Max Projection for Z... ');
    tic;
    if     isfield(userdata,'himp_proj')...
            && ~isempty(userdata.himp_proj) ...
            && isempty(userdata.himp_proj.getWindow()) ...
            && ~isempty(userdata.himp_proj.getImage())  % proj is deleted manually
        userdata.himp_proj.show();
    else
        userdata.himp_proj = []; % delete handle; initialize
        userdata.hmrm.imp_proj = []; % delete handle; initialize
        mmpfivs = MyMaxProjFileInfoVirtualStack(userdata.mfivs); % max projection
        %         saveMaxProjROI(); % save max-projected ROIs
        himp_proj = mmpfivs.getImage();
        userdata.himp_proj = himp_proj;
        setUserData(userdata);
        userdata.himp_proj.show();
    end
    userdata.hmrm.setProj(userdata.himp_proj);
else
    fprintf('processing ... ');
    tic;
    if isfield(userdata,'himp_proj') && ~isempty(userdata.himp_proj)
        userdata.himp_proj.hide();
    end
end

userdata.himp_proj.updateAndDraw();

if origstate ~= MyOrthogonal_Views.isOrthoViewsImage(himp)
    ij.WindowManager.setCurrentWindow(himp.getWindow());
    hmov = javaObjectEDT(MyOrthogonal_Views());
    hmov.run('');
end

addCanvasListener();

enableGUI(true);

fprintf('finished in %g (sec) \n\n',toc);

end


function checkbox_showRoi_callback(~,~)

userdata = getUserData();
flag_showRoi = get(findobj('Tag','check_main_showroi'),'Value')==1;

userdata.hmrm.showAll(flag_showRoi);

huicontrol_tmp(1) = findobj('Tag','check_main_showname');
huicontrol_tmp(2) = findobj('Tag','check_main_showname_estimated');
huicontrol_tmp(3) = findobj('Tag','check_main_showname_fixed');
huicontrol_tmp(4) = findobj('Tag','popup_main_strokecolor');
huicontrol_tmp(5) = findobj('Tag','popup_main_selectioncolor');

if flag_showRoi
    set(huicontrol_tmp,'Enable','on');
else
    set(huicontrol_tmp,'Enable','off');
end

end


function checkbox_showName_callback(~,~)

userdata = getUserData();
flag_showName = get(findobj('Tag','check_main_showname'),'Value')==1;
flag_showNameEstimated = get(findobj('Tag','check_main_showname_estimated'),'Value')==1;
flag_showNameFixed = get(findobj('Tag','check_main_showname_fixed'),'Value')==1;
userdata.hmrm.showName(flag_showName);
userdata.hmrm.setFlagUseNameDisp(flag_showNameEstimated);
userdata.hmrm.setFlagUseFixedName(flag_showNameFixed);

huicontrol_tmp(1) = findobj('Tag','check_main_showname_estimated');
huicontrol_tmp(2) = findobj('Tag','check_main_showname_fixed');

if flag_showName
    set(huicontrol_tmp,'Enable','on');
    if ~flag_showNameEstimated
        set(huicontrol_tmp(2),'Enable','off');
    end
else
    set(huicontrol_tmp,'Enable','off');
end

end


function popupmenu_roicolor_callback(~,~)

huicontrol_strokecolor    = findobj('Tag','popup_main_strokecolor');
huicontrol_selectioncolor = findobj('Tag','popup_main_selectioncolor');

strColors1 = get(huicontrol_strokecolor,   'String');
strColors2 = get(huicontrol_selectioncolor,'String');
strStrokeColor    = strColors1{get(huicontrol_strokecolor,   'Value')};
strSelectionColor = strColors2{get(huicontrol_selectioncolor,'Value')};

strokeColor = java.awt.Color.(strStrokeColor);
selectionColor = java.awt.Color.(strSelectionColor);

userdata = getUserData();
userdata.hmrm.setColor(strokeColor,selectionColor);

end


function checkbox_mindist_callback(~,~)
userdata = getUserData();

if ~isfield(userdata,'jtable') || isempty(userdata.jtable)
    return;
end

% remove existing filter
ftm = userdata.jtable.getTableHeader.getFilterableTableModel();
target_column = userdata.hmtm.getColumnIndex('min_dist');
tmpfilters = ftm.getFilters(target_column);
for p = 1:numel(tmpfilters)
    if isa(tmpfilters(p),'com.jidesoft.filter.LessOrEqualFilter')
        ftm.removeFilter(target_column,tmpfilters(p));
    end
end

% add filter
if get(findobj('Tag','check_main_mindist'),'Value')==1 % filter on
    thr_mindist = str2double(get(findobj('Tag','edit_option_check'),'String'));
    target_filter = com.jidesoft.filter.LessOrEqualFilter(...
        ['<=',num2str(thr_mindist)],thr_mindist);
    ftm.addFilter(target_column,target_filter);
end

% do filtering
ftm.setFiltersApplied(true);
ftm.refresh();

end


function rotate_callback(~,~,rotaxis)
userdata = getUserData();
if ~isfield(userdata,'himp') ...
        || isempty(userdata.himp) ...
        || isempty(userdata.himp.getWindow())
    return;
end
himp = userdata.himp;

enableGUI(false);

tic;
switch rotaxis
    case {'x','X'}
        fprintf('rotate 180 degree around X-axis... ');
        flipXYZ(~[true,false,false]);
        
    case {'y','Y'}
        fprintf('rotate 180 degree around Y-axis... ');
        flipXYZ(~[false,true,false]);
        
    case {'z','Z'}
        fprintf('rotate 180 degree around Z-axis... ');
        flipXYZ(~[false,false,true]);
        
    case {'mirror'}
        fprintf('Mirroring (flip around Z-axis) ... ');
        mirrorZ();
        
    case {'dye'}
        fprintf('rotate around X-axis based on named rois with dye ... ');
        
        % load reference
        flag_findref = false;
        strfile_ref = 'reference.mat';
        strfile_ref_ext = fullfile(userdata.PathName,strfile_ref);
        if exist(strfile_ref_ext,'file')
            fprintf('\n loading: %s',strfile_ref_ext);
            data_ref = load(strfile_ref_ext);
            flag_findref = true;
        elseif exist(strfile_ref,'file')
            fprintf('\n loading: %s',strfile_ref);
            data_ref = load(strfile_ref);
            flag_findref = true;
        else
            fprintf('\n reference data (%s) was not found \n',strfile_ref);
        end
        
        if flag_findref
            pos_ref = data_ref.posref;
            weight  = data_ref.weight;
            
            % collect positions of required neurons
            hmtm = userdata.hmtm;
            idx_frame = hmtm.getValuesAt(0:hmtm.getRowCount()-1,'frame')';
            flag_target = idx_frame==himp.getT();
            roinames = cell(hmtm.getValuesAt(flag_target,'Name'));
            roiname_ref = {'ADLL','ADLR','ASHL','ASHR','ASIL','ASIR',...
                'ASJL','ASJR','ASKL','ASKR','AWBL','AWBR'};
            numref = numel(roiname_ref);
            pos_target = hmtm.getValuesAt(flag_target,'pos')';
            pos_sub = nan(3,numref);
            for p=1:numref
                tmpflag = strcmpi(roinames,roiname_ref(p));
                if any(tmpflag)
                    pos_sub(:,p) = pos_target(:,tmpflag);
                end
            end
            
            %%% if a part of rois is not found, do not use it
            flag_roiexist = ~any(isnan(pos_sub),1);
            pos_sub = pos_sub(:,flag_roiexist);
            pos_ref = pos_ref(:,flag_roiexist);
            numref = sum(flag_roiexist);
            
            % obtain rotation parameter
            x0 = [   ones(1,3),    zeros(1,4)];
            lb = [  zeros(1,3),-pi, -inf(1,3)];
            ub = [10*ones(1,3), pi,  inf(1,3)];
            pos_sub_mean = mean(pos_sub,2);
            pos_ref_mean = mean(pos_ref,2);
            pos_sub_0 = pos_sub - pos_sub_mean(:,ones(1,numref));
            pos_ref_0 = pos_ref - pos_ref_mean(:,ones(1,numref));
            minx = lsqnonlin(...
                @(x) partial_affine(x,pos_sub_0,pos_ref_0,weight(:,flag_roiexist)),...
                x0,lb,ub,optimset('display','none'));
            
            %%% rotation; only rotation for x is used.
            set(findobj('Tag','edit_main_rotate_arbit'),'String',...
                num2str(minx(4)*180/pi));
            rotateX();
        end
        
    case {'arbit_edit'}
        fprintf('rotate arbitrarily around X-axis... ');
        deg = str2double(get(findobj('Tag','edit_main_rotate_arbit'),'String'));
        set(findobj('Tag','slider_main_rotate_arbit'),'Value',deg);
        rotateX();
        
    case {'arbit_slider'}
        fprintf('rotate arbitrarily around X-axis... ');
        deg = get(findobj('Tag','slider_main_rotate_arbit'),'Value');
        set(findobj('Tag','edit_main_rotate_arbit'),'String',num2str(deg));
        rotateX();
end



% update image
if himp.isComposite()
    tmpc = himp.getC();
    tmpluts = himp.getLuts();
    himp.setChannelsUpdated();
    himp.setLuts(tmpluts);
    himp.updateChannelAndDraw();
    himp.setC(tmpc);
end

enableGUI(true);

fprintf('finished in %g (sec) \n\n',toc);

end


function checkbox_align_callback(~,~)

userdata = getUserData();
if ~isfield(userdata,'cum2') || isempty(userdata.cum2)  % calc cum2
    calc_alignment();
    userdata = getUserData();
end
cum2 = userdata.cum2;
flag_align = get(findobj('Tag','check_main_align'),'Value')==1;

enableGUI(false);


hmtm = userdata.hmtm;
mfivs = userdata.mfivs;

%%% apply change
if cellfun(@isempty,regexp(cell(userdata.himp.getTitle()),'_aligned\.tif$','once'))
    % original image (not aligned) is shown
    mfivs.setParallelShift(cum2);
    mfivs.setAlign(flag_align);
else % aligned image is shown
    mfivs.setParallelShift(-cum2);
    mfivs.setAlign(~flag_align);
end

hmtm.setAlign(flag_align);
userdata.frame.repaint();
userdata.hmrm.invalidateAll(); % force update roi display

enableGUI(true);

end


function checkbox_subtract_callback(~,~)
userdata = getUserData();
flag_subtract = get(findobj('Tag','check_main_subtract'),'Value')==1;
userdata.mfivs.setSubBG(flag_subtract);
end


function edit_subtract_radius_callback(~,~)
userdata = getUserData();
radius = str2double(get(findobj('Tag','edit_main_subtract'),'String'));
userdata.mfivs.setRadiusSubBG(radius);
end


%%% ROI detection functions

function load_mif_callback(~,~)
userdata = getUserData();

%%% clear last session
strPathName = [];
if isfield(userdata,'frame') || isfield(userdata,'himp')
    answer = questdlg('Do you close current session?',...
        'Closing session','Yes','No','No');
    drawnow();
    if strcmpi(answer,'Yes')
        strPathName = userdata.PathName;
        if isfield(userdata,'frame'); userdata.frame.dispose(); end
        if isfield(userdata,'himp'); userdata.himp.changes = false; end
        if isfield(userdata,'himp_proj'); userdata.himp_proj.changes = false; end
        ij.IJ.run('Close All');
        close(gcf);
        clear userdata;
        
        %%% create new GUI
        create_GUI();
    else
        return;
    end
end
enableGUI(false);

%%% create and save temporary data
[FileName,PathName] = uigetfile('*.tif;*.mif;*.dcv;*.sif;*.pif','Select image file',strPathName);
if  isequal(FileName,0) || isequal(PathName,0) % canceled
    return;
end
% imname = strrep(FileName,'.tif','');
[~,imname,imext] = fileparts(FileName); %#ok<ASGLU>
strSaveTmp = fullfile(PathName,[imname,'.mat']);
params_em = [1,1,1,1,1,0,0,1,0,1]'; %#ok<NASGU>
idx_frame = 1; %#ok<NASGU>
save(strSaveTmp,'imname','imext','params_em','idx_frame');

%%% load temporary data as usual
data_load_callback([],[],strSaveTmp);

end


function initDetectPreset()
userdata = getUserData();
hobj = findobj('Tag','popup_detect_main_preset');
flist = dir(fullfile(fileparts(mfilename('fullpath')),'peak_detection_14_param_for_*.txt'));
tmpstr = regexp({flist.name},'peak_detection_14_param_for_(.*)\.txt','tokens');
hobj.String = [cellfun(@char,[tmpstr{:}],'UniformOutput',false),{'Default','(Load from file)'}];
userdata.detoptlist = flist;
userdata.detopt = [];
userdata.detoptrel = { % relationship between variable names and tags, and its species.
    %     'flag_debug','check_detect_main_debug','logical';
    'align_filter_radius','edit_detect_align_filter','scalar';
    'align_max_move','edit_detect_align_max','vector2';
    'flag_align_centroid','check_detect_align_centroid','logical';
    'flag_align_subtract','check_detect_align_subtract','logical';
    'channel_use','edit_detect_align_channel','scalar';
    'frame_use','edit_detect_align_frame','scalar';
    'prefilter_method','edit_detect_detect_prefilter_method','string';
    'prefilter_option','edit_detect_detect_prefilter_option','string';
    'radius_background','edit_main_subtract','scalar';
    'blur_method','edit_detect_detect_blur_method','string';
    'blur_option','edit_detect_detect_blur_option','string';
    'threshold_method','edit_detect_detect_threshold','string';
    'radius_dilation','edit_detect_detect_local','vector3';
    'min_voxel_remove','edit_detect_detect_min','scalar';
    'curvature_distance','edit_detect_detect_curvature_distance','vector3';
    'curvature_numvoxels','edit_detect_detect_curvature_#voxels','scalar';
    'removefp_numloop','edit_detect_fit_removefp_numloop','scalar';
    'removefp_min_dist','edit_detect_fit_removefp_mindist','scalar';
    'removefp_min_scdist','edit_detect_fit_removefp_minscdist','scalar';
    'scaling','edit_option_scaling','vector3';
    'sigma_default','edit_option_sigma','vector3';
    'thrdist','edit_track_option_thrdist','scalar';
    'reltol','edit_track_option_tol','scalar';
    };
setUserData(userdata);

% hobj.Value = find(~cellfun(@isempty,regexp({flist.name},'confocal')),1,'first');
hobj.Value = find(strcmp(hobj.String,'Default'),1,'first');

end


function loadDetectPreset()
userdata = getUserData();
hobj = findobj('Tag','popup_detect_main_preset');
if hobj.Value<=numel(userdata.detoptlist)
    stropttxt = userdata.detoptlist(hobj.Value).name;
elseif hobj.Value == numel(userdata.detoptlist)+1 % Default
    if isfield(userdata,'detopt') && ~isempty(userdata.detopt)
        setDetectOpt(userdata.detopt); % set parameters from userdata
    end
    return;
else % load new parameter file
    [FileName,PathName] = uigetfile('*.txt','Select parameter file',userdata.PathName);
    stropttxt = fullfile(PathName,FileName);
end

% load option from file
fid = fopen(stropttxt,'rt');
eval(fread(fid,inf,'*char'));
fclose(fid);

% set option value to UI
detoptrel = userdata.detoptrel;
for p=1:size(detoptrel,1)
    hobj2 = findobj('Tag',detoptrel{p,2});
    switch detoptrel{p,3}
        case 'logical'
            hobj2.Value = eval(detoptrel{p,1});
        case 'scalar'
            hobj2.String = num2str(eval(detoptrel{p,1}));
        case 'vector2'
            hobj2.String = sprintf('[%g,%g]',eval(detoptrel{p,1}));
        case 'vector3'
            hobj2.String = sprintf('[%g,%g,%g]',eval(detoptrel{p,1}));
        case 'string'
            hobj2.String = eval(detoptrel{p,1});
    end
end

if isfield(userdata,'himp')
    setScale();
    edit_subtract_radius_callback();
end

end


function setDetectOpt(detopt)
userdata = getUserData();

% set option value to UI
detoptrel = userdata.detoptrel;
for p=1:size(detoptrel,1)
    hobj2 = findobj('Tag',detoptrel{p,2});
    switch detoptrel{p,3}
        case 'logical'
            hobj2.Value = detopt.(detoptrel{p,1});
        case 'scalar'
            hobj2.String = num2str(detopt.(detoptrel{p,1}));
        case 'vector2'
            hobj2.String = sprintf('[%g,%g]',detopt.(detoptrel{p,1}));
        case 'vector3'
            hobj2.String = sprintf('[%g,%g,%g]',detopt.(detoptrel{p,1}));
        case 'string'
            hobj2.String = detopt.(detoptrel{p,1});
    end
end

if isfield(userdata,'himp')
    setScale();
    edit_subtract_radius_callback();
end

end


function detopt = getDetectOpt()
userdata = getUserData();

%%% gather option values from UI
detopt = userdata.detopt;
detoptrel = userdata.detoptrel;
for p=1:size(detoptrel,1)
    hobj2 = findobj('Tag',detoptrel{p,2});
    switch detoptrel{p,3}
        case 'logical'
            detopt.(detoptrel{p,1}) = hobj2.Value==1;
        case 'scalar'
            detopt.(detoptrel{p,1}) = str2double(hobj2.String);
        case {'vector2','vector3'}
            detopt.(detoptrel{p,1}) = str2num(hobj2.String); %#ok<ST2NM>
        case 'string'
            detopt.(detoptrel{p,1}) = hobj2.String;
    end
end

detopt.flag_debug = false;
detopt.flag_virtual_stack = true;
detopt.flag_reslice_z = get(findobj('Tag','check_main_interpz'),'Value')==1;

userdata.detopt = detopt;
setUserData(userdata);

end


function detectROI()
userdata = getUserData();
enableGUI(false);
pathImage = fullfile(userdata.PathName,[userdata.imname,userdata.imext]);
detopt = getDetectOpt();
peak_detection_14_3(pathImage,detopt);
strSaveTmp = fullfile(userdata.PathName,[userdata.imname,'.mat']);
data_load_callback([],[],strSaveTmp);
end


function setupBlurredImage()
%%% apply settings of aligned image to blurred image
userdata = getUserData();
mfivs = userdata.mfivs;
mfivs_blurred = userdata.mfivs_blurred;
mfivs_blurred.setParallelShift(mfivs.parallelShift);
mfivs_blurred.setAlign(mfivs.flagAlign);
mfivs_blurred.setFlipX(mfivs.flagFlipX);
mfivs_blurred.setFlipY(mfivs.flagFlipY);
mfivs_blurred.setFlipZ(mfivs.flagFlipZ);
mfivs_blurred.zscale = mfivs.zscale;
mfivs_blurred.setInterpZ(mfivs.flagInterpZ);
mfivs_blurred.setRotation(mfivs.getRotation);
mfivs_blurred.clearAndResizeCache(mfivs.zscale*mfivs.nzorig); % to avoid memory-consuming problem
setUserData(userdata);
end


function im_blurred = getBlurredImage(z,t)
%%% return blurred image when the aligned image is used as blurred mage
userdata = getUserData();
himp_blurred = userdata.himp_blurred;

flagBlurredImage = ~cellfun(@isempty,...
    regexp(cell(himp_blurred.getTitle()),'_blurred\.tif$','once'));
if ~flagBlurredImage
    % original image (not blurred) is shown
    strPrefilterMethod = ...
        get(findobj('Tag','edit_detect_detect_prefilter_method'),'String');
    strPrefilterOption = ...
        get(findobj('Tag','edit_detect_detect_prefilter_option'),'String');
    strRadiusBackground = ...
        get(findobj('Tag','edit_main_subtract'),'String');
    strBlurMethod = ...
        get(findobj('Tag','edit_detect_detect_blur_method'),'String');
    strBlurOption = ...
        get(findobj('Tag','edit_detect_detect_blur_option'),'String');
    %     strThresholdMethod = ...
    %         get(findobj('Tag','edit_detect_detect_threshold'),'String');
    idxc = str2double(get(findobj('Tag','edit_detect_align_channel'),'String'));
    himp_tmp = setImageMiji(getImageMiji(himp_blurred,[],idxc,t),'blurred');
    ij.IJ.run(himp_tmp,strPrefilterMethod,strPrefilterOption);
    ij.IJ.run(himp_tmp,'Subtract Background...',['rolling=',strRadiusBackground,' stack']);
    ij.IJ.run(himp_tmp,strBlurMethod,strBlurOption);
    im_blurred = double(getImageMiji(himp_tmp,z,[],[]));
else
    im_blurred = double(getImageMiji(himp_blurred,z,[],t));
end

end



%%% inner functions


function init_ellipsoid_roi()
%%% initialize function for update_elipsoid_roi().
userdata = getUserData();

hmov = MyOrthogonal_Views.getInstance();
if isempty(hmov) % for xyt image
    userdata.hmrm = javaObjectEDT(MyRoiManager2(userdata.himp,userdata.hmtm));
else % for xyz and xyzt image
    userdata.hmrm = javaObjectEDT(MyRoiManager2(javaObjectEDT(hmov),userdata.hmtm));
end

setUserData(userdata);

popupmenu_roicolor_callback();
checkbox_showRoi_callback();
checkbox_showName_callback();
setNsigma();

end


function saveMaxProj()
userdata = getUserData();
disp('saveMaxProj');

zoom   = str2double(get(findobj('Tag','edit_movie_zoom'),'String'));
tStart = str2double(get(findobj('Tag','edit_movie_start'),'String'));
tEnd   = str2double(get(findobj('Tag','edit_movie_end'),'String'));
flagSaveRoi       = get(findobj('Tag','check_movie_roi'),'Value')==1;
flagSaveMovie     = get(findobj('Tag','check_movie_movie'),'Value')==1;
flagSaveOverlay   = get(findobj('Tag','check_movie_overlay'),'Value')==1;

if ~flagSaveRoi && ~flagSaveMovie && ~ flagSaveOverlay
    disp('All checkbox are empty; Please specify items for saving and retry. Aborting...');
    return
end

strFileBase = sprintf('%s_MaxProj_zoom=%d_t=%d-%d',userdata.FileName(1:end-4),zoom,tStart,tEnd);
pathSaveRoi     = fullfile(userdata.PathName,[strFileBase,'_RoiSet.zip']);
pathSaveMovie   = fullfile(userdata.PathName,[strFileBase,'_RGB.tif']);
pathSaveOverlay = fullfile(userdata.PathName,[strFileBase,'_RGB_WithROI.tif']);

dims = userdata.mfivs.getImage().getDimensions(); % [x,y,c,z,t]
tEnd = min(tEnd,dims(5));
numT = tEnd-tStart+1;
dimsNew = [dims(1)*zoom,dims(2)*zoom,1,1,numT];

hrm = ij.plugin.frame.RoiManager.getInstance();
if isempty(hrm)
    hrm = ij.plugin.frame.RoiManager();
end
if (hrm.getCount()>0)
    hrm.runCommand('Deselect');
    hrm.runCommand('Delete');
end
hrm.hide();

if flagSaveRoi
zos = java.util.zip.ZipOutputStream(java.io.BufferedOutputStream(...
    java.io.FileOutputStream(pathSaveRoi)));
out = java.io.DataOutputStream(java.io.BufferedOutputStream(zos));
re = ij.io.RoiEncoder(out);
end

if flagSaveMovie
    mfsMovie = MyFileSaver();
    mfsMovie.saveAsTiffStack_open(pathSaveMovie,dimsNew,ij.io.FileInfo.RGB);
end

if flagSaveOverlay
    mfsOverlay = MyFileSaver();
    mfsOverlay.saveAsTiffStack_open(pathSaveOverlay,dimsNew,ij.io.FileInfo.RGB);
end

if flagSaveMovie || flagSaveOverlay
%     himp = userdata.mfivs.getImage();
    mmpfivs = MyMaxProjFileInfoVirtualStack(userdata.mfivs);
    himp_proj = mmpfivs.getImage();
    flagChannelActive = himp_proj.getActiveChannels;
    idxChannelActive = find(flagChannelActive(1:dims(3)));
%     himp_proj.show();
%     userdata.hmrm.setProj(userdata.himp_proj);
    getLuts(); % myluts saved in userdata
    vWidth  = ( 0:double(dimsNew(1))-1) / zoom+1;
    vHeight = ( 0:double(dimsNew(2))-1) / zoom+1;
    vCol = 1:double(dims(3));
end

ht = tic;

for ct = tStart:tEnd
    hme = userdata.hmrm.get(ct).values().toArray();
    for p=1:numel(hme)
        hroi = hme(p).getMyEllipseAt(1,3); % position=1, mode=3
        strRoiName = hroi.getName();
        strFileName = [cell2mat(cell(strRoiName)),'-',num2str(ct),'.roi'];        
        hroi = ij.plugin.RoiScaler.scale(hroi,zoom,zoom,false);
        hroi.setName(strRoiName);
        hroi.setPosition(ct); % set Slice;        
        if flagSaveRoi
            zos.putNextEntry(java.util.zip.ZipEntry(strFileName));
            re.write(hroi);
            out.flush();
        end
        hrm.addRoi(hroi)
    end
    
    if flagSaveMovie || flagSaveOverlay
%         imMax = max(getImageMiji(himp,[],[],ct),[],3);
        imMax = zeros(dims(1),dims(2),dims(3));
        imMax(:,:,idxChannelActive) = getImageMiji(himp_proj,1,idxChannelActive,ct);
        g = griddedInterpolant(imMax);
        impMax = setImageMiji(permute(g({vWidth,vHeight,vCol}),[1,2,4,3]),'imMax');
        setLuts(impMax);
        impMaxRGB = ij.ImagePlus('imMaxRGB',ij.process.ColorProcessor(impMax.getImage()));
        ij.WindowManager.setTempCurrentImage(impMaxRGB);
    end
    
    if flagSaveMovie
        mfsMovie.saveAsTiffStack_append(impMaxRGB.getProcessor().getPixels(),1);
    end
    
    if flagSaveOverlay        
        hrm.runCommand('Deselect');
        hrm.runCommand('Draw');
        mfsOverlay.saveAsTiffStack_append(impMaxRGB.getProcessor().getPixels(),1);
    end
    
    if flagSaveMovie || flagSaveOverlay
        impMax.close();
        impMaxRGB.close();
    end
    
    hrm.runCommand('Deselect');
    hrm.runCommand('Delete');
    
    fprintf('Frame #%d processed.\n',ct);
end

hrm.close();

if flagSaveRoi
out.close();
disp(['Max Projected ROIs were saved to: ',pathSaveRoi]);
end

if flagSaveMovie || flagSaveOverlay
    himp_proj.close();
end

if flagSaveMovie
mfsMovie.saveAsTiffStack_close();
disp(['Max Projected Movie was saved to: ',pathSaveMovie]);
end

if flagSaveOverlay
mfsOverlay.saveAsTiffStack_close();
disp(['Max Projected Movie with ROIs Overlaied was saved to: ',pathSaveOverlay]);
end

fprintf('elapsed time: %f sec\n',toc(ht));
end


function enableGUI(flag)

userdata = getUserData();

if flag
    set(findobj('Enable','off'),'Enable','on');
else
    set(findobj('Enable','on'),'Enable','off');
end

set(findobj('-regexp','Tag','^(load|save)$'),'Enable','on');
set(findobj('-regexp','Tag','button_detect_main_load'),'Enable','on');

if isfield(userdata,'himp')
    if userdata.himp.getNSlices()==1
        %%% single z-slice image
        huicontrol_single_z = findobj('Tag','check_main_proj');
        set(huicontrol_single_z,'Enable','off');
    end
    
    if userdata.himp.getNFrames()==1
        %%% single t-frame image
        set(findobj('-regexp',...
            'Tag','button_track_(source|target)_(?!(default|current))'),...
            'Enable','off');
    end
    
    %     if userdata.flag_virtual_stack
    %         huicontrol_virtual = findobj('Tag','check_main_logarithm');
    %         huicontrol_virtual = cat(1,huicontrol_virtual,...
    %             findobj('-regexp','Tag','button_main_rotate_'));
    %         set(huicontrol_virtual,'Enable','off');
    %     end
end

drawnow();

end


function myluts = getLuts()
%%% save luts for display
userdata = getUserData();
himp = userdata.himp;
luts = himp.getLuts();
if isempty(luts) % single channel image
    luts = himp.getProcessor().getLut();
end
for p=1:numel(luts)
    tmplut = luts(p).getBytes();
    myluts(p).r = typecast(tmplut(  1:256),'uint8');  %#ok<AGROW>
    myluts(p).g = typecast(tmplut(257:512),'uint8');  %#ok<AGROW>
    myluts(p).b = typecast(tmplut(513:768),'uint8');  %#ok<AGROW>
    myluts(p).min = luts(p).min;  %#ok<AGROW>
    myluts(p).max = luts(p).max;  %#ok<AGROW>
end

% if called with output argument, do not copy myluts to userdata, and
% return it directly.
if nargout==0
    userdata.myluts = myluts;
    setUserData(userdata);
end

end


function setLuts(himp)
%%% recover luts and apply
userdata = getUserData();
myluts = userdata.myluts;
if isempty(myluts); return; end % if empty, do nothing

luts_r = [];
for p=1:numel(myluts)
    luts_r = [luts_r,ij.process.LUT(myluts(p).r,myluts(p).g,myluts(p).b)]; %#ok<AGROW>
    luts_r(p).min = myluts(p).min;
    luts_r(p).max = myluts(p).max;
end
% luts_r = [luts_r,ij.process.LUT(zeros(1,256,'uint8'),uint8(0:255),zeros(1,256,'uint8'))];
% luts_r(end).min = 0;
% luts_r(end).max = 1;
if himp.isComposite()
    himp.setLuts(luts_r);
else
    himp.getProcessor().setLut(luts_r(1));
end

end


function setScale()
%%% set scaling parameter
userdata = getUserData();
himp = userdata.himp;
scaling = str2num(get(findobj('Tag','edit_option_scaling'),'String')); %#ok<ST2NM>
hcal = himp.getCalibration();
hcal.pixelWidth  = scaling(1);
hcal.pixelHeight = scaling(2);
hcal.pixelDepth  = scaling(3);
hcal.setUnit('micron');
himp.setCalibration(hcal);

userdata.mfivs.zscale = scaling(3)/scaling(2);
userdata.hmtm.zscale = scaling(3)/scaling(2);

if MyOrthogonal_Views.isOrthoViewsImage(himp)
    hmov = javaObjectEDT(MyOrthogonal_Views.getInstance());
    hmov.imageClosed(himp);
    ij.WindowManager.setCurrentWindow(himp.getWindow());
    hmov = javaObjectEDT(MyOrthogonal_Views());
    hmov.run('');
end

end


function autosave()
%%% save data automatically

hbutton = findobj('Tag','button_anno_run');

if get(findobj('Tag','check_anno_auto'),'Value')==1 % autoupdate
    if hbutton.Value==0 % keep avoiding from infinite loop of autoupdate
        hbutton.Value = 1;
        hbutton.Callback(); % this update induce the autosave callback
    else
        hbutton.Value = 0;
    end
else
    hbutton.Value = 0;
end

if get(findobj('Tag','check_main_autosave'),'Value')==1 % autosave is valid
    data_save_callback([],[],true);
end

end


function res = partial_affine(x,pos_sub,pos_ref,weight)
% return residuals after transformation of scaling, shift, x-rotation.

%%% Let positions of reference animal as y, and that of sample
%%% animal as x. Then registration of x to y by Affine
%%% transformation matrix A is
%%%   A*x = y,
%%% and matrix A can be obtained by
%%%   A = y/x.
%%% Matlab function rdivide ('/') returns a solution that
%%% minimize norm(A*x-y). Let introduce weight w (diagonal
%%% matrix) for residuals. w means distribution of positions in
%%% reference animal(s). Weighted sum of square of residuals is
%%%   norm( (A*x-y)*w ) = norm( A*x*w - y*w ),
%%% and we can obtain affine matrix A for the weighted sum by
%%%   A = (y*w) / (x*w).
%%%

scale = diag(x(1:3));
rotx = [1,0,0;0,cos(x(4)),-sin(x(4));0,sin(x(4)),cos(x(4))];
shift = x(5:7);
affinemat = [rotx*scale,shift(:);0,0,0,1];
tmpones = ones(1,size(pos_ref,2));
tmpres = affinemat*[pos_sub;tmpones]-[pos_ref;tmpones];
res = reshape(tmpres(1:3,:).*weight,1,[]);
end


function show_popup_menu_callback(hObj,ev)
%%% obtain names from Name_estim and Name_anno_human, and create popup menu

userdata = getUserData();
hmtm = userdata.hmtm;
% htrs = userdata.jtable.getRowSorter();
hqtff = userdata.hqtff;
jtable = userdata.jtable;

hjpm = javax.swing.JPopupMenu(); % new JPopupMenu
hlistener = [];
% modelRows = getSelectedRowIndexInTableAsModel();
modelRows = getSelectedRowIndexInFigureAsModel();
% col_name_popup = hmtm.getColumnIndex('Name_popup');

% undo & redo
hjmi = javax.swing.JMenuItem('UNDO (Ctrl+Z)');
hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',@(~,~)undo));
hjpm.add(hjmi);
hjmi = javax.swing.JMenuItem('REDO (Ctrl+Y)');
hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',@(~,~)redo));
hjpm.add(hjmi);
hjpm.addSeparator();

% rename column header
if isa(hObj,'com.jidesoft.grid.SortableTable')
    th = hObj.getTableHeader();
    cm = th.getColumnModel();
    ci = cm.getColumnIndexAtX(ev.getX);
    tc = cm.getColumn(ci);
    strValue = tc.getHeaderValue();
    strId = tc.getIdentifier();
    hjmi = javax.swing.JMenuItem(['Rename Column Header of ',strValue]);
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',@(~,~)renameColumnHeader(strId)));
    hjpm.add(hjmi);
    hjmi = javax.swing.JMenuItem('Reset All Column Header');
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',@(~,~)renameColumnHeader));
    hjpm.add(hjmi);
    hjpm.addSeparator();
end

% add ROI
if isa(hObj,'ij.gui.ImageCanvas')
    hjmi = javax.swing.JMenuItem('Add New ROI Here (Ctrl+Click)');
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
        {@(~,~,x,y)canvas_mouse_generate(x,y),hObj,ev}));
    hjpm.add(hjmi);
end

% move ROI
if numel(modelRows)==1 && isa(hObj,'ij.gui.ImageCanvas')
    name = char(hmtm.getValuesAt(modelRows,'Name'));
    hjmi = javax.swing.JMenuItem(['Move ',name,' Here (Ctrl+Alt+Click)']);
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
        {@(~,~,x,y)canvas_mouse_move(x,y),hObj,ev}));
    hjpm.add(hjmi);
end

% remove ROI
if numel(modelRows)>=1
    hjmi = javax.swing.JMenuItem('Remove Selected ROIs (Del)');
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
        @(~,~)removeROIs));
    hjpm.add(hjmi);
end

% switch ROI name
if numel(modelRows)==2
    names = cell(hmtm.getValuesAt(modelRows,'Name'));
    hjmi = javax.swing.JMenuItem(['Switch Name: ',names{1},' <-> ',names{2}]);
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
        {@(~,~,x,y)table_switch_roiname(x,y),modelRows(1)+1,modelRows(2)+1}));
    hjpm.add(hjmi);
end

% add selection
if isa(hObj,'ij.gui.ImageCanvas')
    [~,flag_near,minidx] = canvas_mouse_setup(hObj,ev);
    if flag_near
        tmpstr = char(hmtm.getValuesAt(minidx-1,'Name'));
        hjmi = javax.swing.JMenuItem(['Add/Clear Selection of ',tmpstr,' (Shift+Click)']);
        hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
            {@(~,~,x,y)canvas_mouse_select(x,y),hObj,ev}));
        hjpm.add(hjmi);
    end
end

% clear selection
if numel(modelRows) > 0
    hjmi = javax.swing.JMenuItem('Clear All Selection (Shift+Alt+Click)');
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',@(~,~)clearSelection));
    hjpm.add(hjmi);
end

hjpm.addSeparator();

if numel(modelRows) > 0
    hjm = javax.swing.JMenu('Narrowing');
    hjpm.add(hjm);
    
    % obtain names
    name = strrep(cell(hmtm.getValuesAt(modelRows,'Name')),'_LIKE','');
    name_estim = regexp(cell(hmtm.getValuesAt(modelRows,'Name_estim')),...
        '([^, ]+)','tokens');
    name_human = regexp(cell(hmtm.getValuesAt(modelRows,'Name_anno_human')),...
        '([^, ]+)','tokens');
    names = unique(cat(2,name(:)',cellfun(@char,cat(2,name_estim{:},name_human{:}),...
        'UniformOutput',false)));
    names = names(~strcmp(names,'0') & ~strcmp(names,'N/A') & ~strcmp(names,'Null'));
    
    % add narrowing
    for p=1:numel(names)
        hjmi = javax.swing.JMenuItem(['Narrowing by ',names{p}]); % filter
        %         hrf = javax.swing.RowFilter.regexFilter(names(p),col_name_popup);
        %         hlistener = cat(1,hlistener,...
        %             addlistener(hjmi,'ActionPerformed',{@(~,~,x)htrs.setRowFilter(x),hrf}));
        hlistener = cat(1,hlistener,...
            addlistener(hjmi,'ActionPerformed',{@(~,~,x)hqtff.setSearchingText(x),names{p}}));
        hjm.add(hjmi);
    end
end

% reset narrowing
hjmi = javax.swing.JMenuItem('Reset Narrowing');
% hlistener = cat(1,hlistener,...
%     addlistener(hjmi,'ActionPerformed',{@(~,~,x)htrs.setRowFilter(x),[]}));
hlistener = cat(1,hlistener,...
    addlistener(hjmi,'ActionPerformed',{@(~,~,x)hqtff.setSearchingText(x),''}));
hjpm.add(hjmi);
hjpm.addSeparator();

% set name as
if numel(modelRows>0)
    hjm = javax.swing.JMenu(['Set ',getColumnHeaderValue(jtable,hmtm,'Name'),' as...']);
    hjpm.add(hjm);
    
    % Approve top rank candidate
    name_estim_top = cellfun(@char,regexp(cell(hmtm.getValuesAt(modelRows,...
        'Name_estim')),'([^, ]+)','tokens','once'),'UniformOutput',false);
    hjmi = javax.swing.JMenuItem('Approve top rank estimation'); % set name as
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
        {@(~,~,x,y)approveEstimatedName(x,y),modelRows,name_estim_top}));
    hjm.add(hjmi);
    
    % if single ROI is selected; display for each candidate
    if numel(modelRows) == 1
        name_estim = cellfun(@char,regexp(char(hmtm.getValuesAt(modelRows,'Name_estim')),...
            '([^, ]+)','tokens'),'UniformOutput',false);
        name_estim_unique = cellfun(@char,regexp(char(hmtm.getValuesAt(modelRows,'Name_estim_unique')),...
            '([^, ]+)','tokens'),'UniformOutput',false);
        name_human = cellfun(@char,regexp(char(hmtm.getValuesAt(modelRows,'Name_anno_human')),...
            '([^, ]+)','tokens'),'UniformOutput',false);
        names = unique(cat(2,name_estim,name_human,name_estim_unique));
        names = names(~strcmp(names,'0') & ~strcmp(names,'N/A') & ~strcmp(names,'Null'));
        %     name_anno_human_old = cell(hmtm.getValuesAt(modelRows,'Name_anno_human'));
        
        % set name as, for each candidate
        %     idxcol_name = hmtm.getColumnIndex('Name');
        for p=1:numel(names)
            % hjmi = javax.swing.JMenuItem(['Set ',strValue,' as ',names{p}]); % set name as
            hjmi = javax.swing.JMenuItem(names{p}); % set name as
            hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
                @(~,~)approveEstimatedName(modelRows,names(p))));
            hjm.add(hjmi);
        end
    end
end

% send anno info to name_anno_human
if numel(modelRows>0)
    hjpm.addSeparator();
    hjm = javax.swing.JMenu(['Send anno info to ',getColumnHeaderValue(jtable,hmtm,'Name_anno_human'),'...']);
    hjpm.add(hjm);
    
    name_old = cell(hmtm.getValuesAt(modelRows,'Name_anno_human'));
    tmpcbfun = @(x)hmtm.setValuesAt(x,modelRows,'Name_anno_human',[true,true,false]);
    
    % send Name
    name_new = strcat(name_old,{', '},cell(hmtm.getValuesAt(modelRows,'Name')));
    hjmi = javax.swing.JMenuItem(getColumnHeaderValue(jtable,hmtm,'Name'));
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',@(~,~)tmpcbfun(name_new)));
    hjm.add(hjmi);
    
    % send Name_estim_unique
    name_new = strcat(name_old,{', '},cell(hmtm.getValuesAt(modelRows,'Name_estim_unique')));
    hjmi = javax.swing.JMenuItem(getColumnHeaderValue(jtable,hmtm,'Name_estim_unique'));
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',@(~,~)tmpcbfun(name_new)));
    hjm.add(hjmi);
    
    % send Name_estim
    name_new = strcat(name_old,{', '},cell(hmtm.getValuesAt(modelRows,'Name_estim')));
    hjmi = javax.swing.JMenuItem(['All of ',getColumnHeaderValue(jtable,hmtm,'Name_estim')]);
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',@(~,~)tmpcbfun(name_new)));
    hjm.add(hjmi);
    
    % send top rank of Name_estim
    name_estim_top = cellfun(@char,regexp(cell(hmtm.getValuesAt(modelRows,...
        'Name_estim')),'([^, ]+)','tokens','once'),'UniformOutput',false);
    name_new = strcat(name_old,{', '},name_estim_top);
    hjmi = javax.swing.JMenuItem(['Top of ',getColumnHeaderValue(jtable,hmtm,'Name_estim')]);
    hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',@(~,~)tmpcbfun(name_new)));
    hjm.add(hjmi);
    
    % if single ROI is selected; display for each candidate
    if numel(modelRows) == 1
        name_estim = cellfun(@char,regexp(char(hmtm.getValuesAt(modelRows,'Name_estim')),...
            '([^, ]+)','tokens'),'UniformOutput',false);
        name_estim_unique = cellfun(@char,regexp(char(hmtm.getValuesAt(modelRows,'Name_estim_unique')),...
            '([^, ]+)','tokens'),'UniformOutput',false);
        name_human = cellfun(@char,regexp(char(hmtm.getValuesAt(modelRows,'Name_anno_human')),...
            '([^, ]+)','tokens'),'UniformOutput',false);
        names = unique(cat(2,name_estim,name_human,name_estim_unique));
        names = names(~strcmp(names,'0') & ~strcmp(names,'N/A') & ~strcmp(names,'Null'));
        
        for p=1:numel(names)
            name_new = strcat(name_old,{', '},names{p});
            hjmi = javax.swing.JMenuItem(['''',names{p},'''']); % set name as
            hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',@(~,~)tmpcbfun(name_new)));
            hjm.add(hjmi);
        end
    end
end

% set multi-selected cell values
if numel(modelRows>0) && isa(hObj,'com.jidesoft.grid.SortableTable')
    th = hObj.getTableHeader();
    cm = th.getColumnModel();
    ci = cm.getColumnIndexAtX(ev.getX);
    tc = cm.getColumn(ci);
    hmtm = userdata.hmtm;
    flagUpdate = [true,true,false];
    strValue = tc.getHeaderValue();
    strId = tc.getIdentifier();
    hjpm.addSeparator();
    hjm = javax.swing.JMenu(['Set ',strValue,'...']);
    hjpm.add(hjm);
    if islogical(hmtm.getValuesAt(0,tc.getModelIndex))
        hjmi = javax.swing.JMenuItem('Toggle');
        hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
            @(~,~)hmtm.setValuesAt(~(hmtm.getValuesAt(modelRows,strId)),modelRows,strId,flagUpdate)));
        hjm.add(hjmi);
        hjmi = javax.swing.JMenuItem('Set True');
        hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
            @(~,~)hmtm.setValuesAt(true(numel(modelRows),1),modelRows,strId,flagUpdate)));
        hjm.add(hjmi);
        hjmi = javax.swing.JMenuItem('Set False');
        hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
            @(~,~)hmtm.setValuesAt(false(numel(modelRows),1),modelRows,strId,flagUpdate)));
        hjm.add(hjmi);
    elseif isa(hmtm.getValuesAt(0,tc.getModelIndex),'java.lang.String[][]')
        if ~strcmp(strId,'Name')
            hjmi = javax.swing.JMenuItem('Set 0');
            hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
                @(~,~)hmtm.setValuesAt(repmat({'0'},numel(modelRows),1),modelRows,strId,flagUpdate)));
            hjm.add(hjmi);
        end
        strUid = cellfun(@num2str,num2cell(hmtm.getValuesAt(modelRows,'uniqueID')),'UniformOutput',false');
        hjmi = javax.swing.JMenuItem('Set UID');
        hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
            @(~,~)hmtm.setValuesAt(strUid,modelRows,strId,flagUpdate)));
        hjm.add(hjmi);
    else % numeric
        if ~strcmp(strId,'uniqueID')
            hjmi = javax.swing.JMenuItem('Set 0');
            hlistener = cat(1,hlistener,addlistener(hjmi,'ActionPerformed',...
                @(~,~)hmtm.setValuesAt(zeros(numel(modelRows),1),modelRows,strId,flagUpdate)));
            hjm.add(hjmi);
        end
    end
end


% show popup menu
hjpm.show(ev.getComponent,ev.getX(),ev.getY());
userdata.hlistener.popup = hlistener;
setUserData(userdata);

end


function strOut = getColumnHeaderValue(jtable,hmtm,strValue)
strOut = strValue;
if jtable.convertColumnIndexToView(hmtm.getColumnIndex(strValue))~=-1 % displayed
    cm = jtable.getColumnModel();
    tc = cm.getColumn(cm.getColumnIndex(strValue));
    strOut = tc.getHeaderValue();
end
end


function approveEstimatedName(modelRows,names_new)

% check overlap between names_new
[names_unique,~,idx_unique] = unique(names_new);
flag_overlap = histcounts(idx_unique)~=1;
if any(flag_overlap)
    names_overlap = strcat(names_unique(flag_overlap),', ');
    strdisp = cat(2,names_overlap{:});
    disp(['Overlapped Name_estim: ',strdisp(1:end-2)]);
    disp('Operation is canceled.');
    return;
end

% temporally disable autosave
huic = findobj('Tag','check_main_autosave');
flag_autosave_orig = huic.Value;
huic.Value = 0;

% set top estimaed name to name for each roi
for p=1:numel(names_new)
    table_switch_roiname(modelRows(p),names_new{p});
end

% recover autosave state and do autosave
huic.Value = flag_autosave_orig;
autosave();

end


function clearSelection(flag_disp)
userdata = getUserData();
userdata.hmrm.clearSelection();
userdata.jtable.clearSelection();
if ~exist('flag_disp','var') || flag_disp
    disp('ROIs were deselected');
end
end


function idxRowModel = getSelectedRowIndexInTableAsModel()
%%% get indices of selected rows in myTableModel
userdata = getUserData();
jtable = userdata.jtable;
selections = jtable.getSelectedRows(); %% terrible bug; existing but wrong line number returned;;;

idxRowModel = -ones(size(selections));
for p=1:numel(selections)
    idxRowModel(p) = com.jidesoft.grid.TableModelWrapperUtils.getActualRowAt(...
        jtable.getModel(),selections(p));
end

idxRowModel(idxRowModel>=userdata.hmtm.getRowCount()) = -1; % avoid error

end


function idxRowModel = getSelectedRowIndexInFigureAsModel()
%%% get indices of selected rows in Figure (= MyRoiManager)
%%% TableModel does not hold selection of invisible rois but myRoiManager do
userdata = getUserData();
hmrm = userdata.hmrm;
idxRowModel = hmrm.getSelectedRowIndices();
end


function idxRowModel = convertRowIndexToModel(idxRowView) %#ok<DEFNU>
userdata = getUserData();
jtable = userdata.jtable;
idxRowModel  = zeros(size(idxRowView));
for p=1:numel(idxRowView)
    idxRowModel(p) = com.jidesoft.grid.TableModelWrapperUtils.getActualRowAt(...
        jtable.getModel(),jtable.convertRowIndexToModel(idxRowView));
end
end


function idxRowView = convertRowIndexToView(idxRowModel)
userdata = getUserData();
jtable = userdata.jtable;
idxRowView  = -ones(size(idxRowModel));
for p=1:numel(idxRowModel)
    idxRowView(p) = com.jidesoft.grid.TableModelWrapperUtils.getRowAt(...
        jtable.getModel(),idxRowModel(p));
end
end


function addCanvasListener()

userdata = getUserData();
himp = userdata.himp;
hcanvas = himp.getCanvas();
hcanvas.disablePopupMenu(true); % disable imagej popup menu when right-clicked
hlistener(1) = addlistener(hcanvas,'MouseReleased',@canvas_mouse_callback);
hlistener(2) = addlistener(hcanvas,  'KeyReleased',@canvas_key_callback);

if himp.getNSlices()>1 % multi z-slice image
    hmov = javaObjectEDT(MyOrthogonal_Views.getInstance());
    while (isempty(hmov.getXZImage()) ||  isempty(hmov.getXZImage()))
        pause(0.1);
    end
    hcanvas_xz = hmov.getXZImage().getCanvas();
    hcanvas_yz = hmov.getYZImage().getCanvas();
    hcanvas_xz.disablePopupMenu(true); % disable imagej popup menu when right-clicked
    hcanvas_yz.disablePopupMenu(true); % disable imagej popup menu when right-clicked
    hlistener(3) = addlistener(hcanvas_xz,'MouseReleased',@canvas_mouse_callback);
    hlistener(4) = addlistener(hcanvas_yz,'MouseReleased',@canvas_mouse_callback);
    hlistener(5) = addlistener(hcanvas_xz,'KeyReleased',@canvas_key_callback);
    hlistener(6) = addlistener(hcanvas_yz,'KeyReleased',@canvas_key_callback);
end

userdata.hlistener.hcanvas = hlistener;

setUserData(userdata);

end


function undo()
userdata = getUserData();
hmtm = userdata.hmtm;
if hmtm.undoStack.isEmpty()
    disp('Undo history is empty.');
else
    hmtm.undo(true);
    disp('Undo');
end
end


function redo()
userdata = getUserData();
hmtm = userdata.hmtm;
if hmtm.redoStack.isEmpty()
    disp('Redo history is empty.');
else
    hmtm.redo(true);
    disp('Redo');
end
end


function renameColumnHeader(strIn)
userdata = getUserData();
jtable = userdata.jtable;
th = jtable.getTableHeader();
cm = th.getColumnModel();
if ~exist('strIn','var') || isempty(strIn) % reset all if strIn is empty
    for p=1:cm.getColumnCount()
        tc = cm.getColumn(p-1);
        tc.setHeaderValue(tc.getIdentifier());
    end
elseif ~iscell(strIn) % rename single column if strIn is char
    ci = cm.getColumnIndex(strIn);
    tc = cm.getColumn(ci);
    strHeaderValue = inputdlg(['Enter new header name of ',strIn],...
        'Rename column header',1,{tc.getHeaderValue()});
    tc.setHeaderValue(strHeaderValue{1});
else % rename all column if strIn is cell array of labels and identifiers
    for p=1:size(strIn,1)
        ci = cm.getColumnIndex(strIn{p,1});
        tc = cm.getColumn(ci);
        tc.setHeaderValue(strIn{p,2});
    end
end

end


function cellColumnHeader = getColumnHeader()
userdata = getUserData();
jtable = userdata.jtable;
th = jtable.getTableHeader();
cm = th.getColumnModel();
numcol = cm.getColumnCount();
cellColumnHeader = cell(numcol,2);
for p=1:numcol
    tc = cm.getColumn(p-1);
    cellColumnHeader(p,1) = {tc.getIdentifier()};
    cellColumnHeader(p,2) = {tc.getHeaderValue};
end
end


function setColumnWidth(cellColumnWidth)
userdata = getUserData();
jtable = userdata.jtable;
th = jtable.getTableHeader();
cm = th.getColumnModel();
for p=1:size(cellColumnWidth,1)
    ci = cm.getColumnIndex(cellColumnWidth{p,1});
    tc = cm.getColumn(ci);
    tc.setPreferredWidth(cellColumnWidth{p,2});
    tc.setWidth(cellColumnWidth{p,2});
end

end


function cellColumnWidth = getColumnWidth()
userdata = getUserData();
jtable = userdata.jtable;
th = jtable.getTableHeader();
cm = th.getColumnModel();
numcol = cm.getColumnCount();
cellColumnWidth = cell(numcol,2);
for p=1:numcol
    tc = cm.getColumn(p-1);
    cellColumnWidth(p,1) = {tc.getIdentifier()};
    cellColumnWidth(p,2) = {tc.getWidth};
end
end


function setEmopt(emopt)
hctrl = findobj('-regexp','Tag','^edit_track_option_*');
strfield = strrep(get(hctrl,'Tag'),'edit_track_option_','');
for p=1:numel(strfield)
    if isfield(emopt,strfield{p})
        set(hctrl(p),'String',num2str(emopt.(strfield{p})));
    end
end
end


function emopt = getEmopt(flag_checkbox)
% for backward compatibility
emopt.lambda_sigma      = 0;
emopt.lambda_mu_coeff   = 1;
emopt.lambda_mu_nondiag = 0;
emopt.lambda_pi         = 0;
emopt.lm_maxiter        = 10;
emopt.lm_damp           = 2;
emopt.lm_init           = 1;
emopt.m1                = 3;
emopt.m2                = 7;

% editbox
hctrl = findobj('-regexp','Tag','^edit_track_option_*');
strfield = strrep(get(hctrl,'Tag'),'edit_track_option_','');
for p=1:numel(strfield)
    emopt.(strfield{p}) = str2double(get(hctrl(p),'String'));
end

if flag_checkbox
    % checkbox
    hctrl = findobj('-regexp','Tag','^check_track_option_*');
    strfield = strrep(get(hctrl,'Tag'),'check_track_option_','');
    for p=1:numel(strfield)
        emopt.(strfield{p}) = get(hctrl(p),'Value');
    end
end

end


function userdata = getUserData()
userdata = get(findobj('Tag','figure_main'),'UserData');
end


function setUserData(userdata)
set(findobj('Tag','figure_main'),'UserData',userdata);
end


function displayCurrentSelection()
userdata = getUserData();
hmrm = userdata.hmrm;
names = cell(hmrm.getSelectedRowNames());
if isempty(names)
    strrois = '(none)';
else
    tmpstrrois = strcat(unique(names),{', '});
    strrois = cat(2,tmpstrrois{:});
    strrois = strrois(1:end-2);
end
disp(['Current Selection: ',strrois]);
end


function flipXYZ(flagFlips)
%%% inner function for flip X, Y, and Z

userdata = getUserData();
mfivs = userdata.mfivs;
hmtm = userdata.hmtm;

if flagFlips(1) % FlipX
    mfivs.toggleFlipX();
    hmtm.toggleFlipX();
end
if flagFlips(2) % FlipY
    mfivs.toggleFlipY();
    hmtm.toggleFlipY();
end
if flagFlips(3) % FlipZ
    mfivs.toggleFlipZ();
    hmtm.toggleFlipZ();
end

userdata.frame.repaint();
userdata.hmrm.invalidateAll(); % force update roi display

end


function rotateX()
%%% inner function for arbitrary rotation along X axis
userdata = getUserData();
mfivs = userdata.mfivs;
hmtm = userdata.hmtm;
rad = str2double(get(findobj('Tag','edit_main_rotate_arbit'),'String'))/180*pi;
mfivs.setRotation(rad);
hmtm.setRotation(rad);
userdata.frame.repaint();
userdata.hmrm.invalidateAll(); % force update roi display
end


function mirrorZ
%%% inner function for flip X, Y, and Z
userdata = getUserData();
mfivs = userdata.mfivs;
hmtm = userdata.hmtm;
mfivs.toggleFlipZ();
hmtm.toggleFlipZ();
userdata.frame.repaint();
userdata.hmrm.invalidateAll(); % force update roi display
end


function calc_alignment()

disp('Calculation of parallel shift amount (alignment)...');
ht2 = tic;

enableGUI(false);

userdata = getUserData();
himp = userdata.himp;
% path_im_orig = fullfile(userdata.PathName,[userdata.imname,'.tif']);
path_im_orig = fullfile(userdata.PathName,[userdata.imname,userdata.imext]);
[~,himp_orig,mfivs_orig] = getImageMiji(path_im_orig,[],[],1);

flag_centroid = get(findobj('Tag','check_detect_align_centroid'),'Value')==1;
flag_subtract = get(findobj('Tag','check_detect_align_subtract'),'Value')==1;
filter_radius = str2double(get(findobj('Tag','edit_detect_align_filter'),'String'));
max_move = str2num(get(findobj('Tag','edit_detect_align_max'),'String')); %#ok<ST2NM>

channel_use = himp.getC();
mini_scale = 4;
flag_use_mini_z = false; % if true, fast but less precise alignment for z

width  = himp_orig.getWidth();
height = himp_orig.getHeight();
depth  = himp_orig.getNSlices();
% numcolor = himp_orig.getNChannels();
numframe = himp_orig.getNFrames();

xm = 1:(width/mini_scale/2);
ym = 1:(height/mini_scale/2);
numxm = numel(xm);
numym = numel(ym);
cmat = false(width,height);
cmat(xm            ,[ym,ym-numym+height]) = true;
cmat(xm-numxm+width,[ym,ym-numym+height]) = true;


% calc alignment with z-slice
output = zeros(4,depth-1,numframe);
filtermatrix = fspecial('average');
fiji_directory = fileparts(fileparts(which('Miji')));
if depth>2
    % align with z-slice
    % get image data with padding
    pad_width  = round( width*0.1);
    pad_height = round(height*0.1);
    % im_pad = zeros(width+pad_width*2,height+pad_height*2);
    range_width  =  (1:width) + pad_width;
    range_height = (1:height) + pad_height;
    
    % fast fourier transform and subpixel registration
    fprintf('\n');
    parfor r=1:numframe
        %     for r=1:numframe
        ht = tic;
        fprintf('calc movement in frame # %d ... ',r);
        
        % for noisy image
        Miji_mod(false,fiji_directory); % setup Miji for parfor environment
        % himp2 = setImageMiji(getImageMiji(path_im_orig,[],channel_use,r),'im2');
        [~,himp2,mfivs2] = getImageMiji(path_im_orig,[],[],1);
        mfivs2.setMappingSize(mfivs2.sliceBytesWithGap*100); % set MappedBuffer Size as 100 image
        mfivs2.clearAndResizeCache(100); % Cache 100 image
        
        if filter_radius~=0
            % ij.IJ.run(himp2,'Median...', ['radius=',num2str(filter_radius),' stack']);
            % ij.IJ.run(himp2,'Gaussian Blur...',['sigma=',num2str(filter_radius),' stack']);
            mfivs2.setMedian(true);
            mfivs2.setBlur(true);
        end
        
        tmpim = double(getImageMiji(himp2,1,channel_use,r));
        im_pad = zeros(width+pad_width*2,height+pad_height*2);
        im_pad(range_width,range_height) ...
            = tmpim - flag_subtract*mode(tmpim(:));
        if flag_subtract
            im_fft_old = fft2(imfilter(im_pad,filtermatrix));
        else
            im_fft_old = fft2(im_pad);
        end
        if flag_use_mini_z
            im_fft_old = reshape(im_fft_old(cmat),[numxm*2,numym*2]); %#ok<UNRCH>
        end
        for p=1:depth-1
            tmpim = double(getImageMiji(himp2,p+1,channel_use,r));
            im_pad(range_width,range_height) ...
                = tmpim - flag_subtract*mode(tmpim(:));
            if flag_subtract
                im_fft_new = fft2(imfilter(im_pad,filtermatrix));
            else
                im_fft_new = fft2(im_pad);
            end
            if flag_use_mini_z
                im_fft_new = reshape(im_fft_new(cmat),[numxm*2,numym*2]); %#ok<UNRCH>
            end
            [output(:,p,r)] = dftregistration(im_fft_old,im_fft_new,100);
            im_fft_old = im_fft_new;
        end
        fprintf('%g (sec) \n',toc(ht));
    end
    clear('im_pad','im_fft_old','im_fft_new');
end

% alignment
if flag_use_mini_z
    output(3:4,:,:) = output(3:4,:,:)*4; %#ok<UNRCH>
end

% for noisy image
flag_max_move = any(bsxfun(@ge,abs(output(3:4,:,:)),max_move(:)));
output(flag_max_move(ones(1,4),:,:))=0;

% temporary aligned in z but currently not for t
cum = cat(2,zeros(2,1,numframe),cumsum(output(3:4,:,:),2));
mfivs_orig.setParallelShift(cum);
mfivs_orig.setAlign(true);

if flag_centroid % image alignment by centroid
    im = getImageMiji(himp_orig,[],channel_use,1);
    im_proj = double(max(im,[],3));
    sum_proj = sum(im_proj(:));
    [width_proj,height_proj] = size(im_proj);
    gx = (1:width_proj);
    gy = (1:height_proj)';
    centroid = zeros(2,numframe);
    centroid(1,1) = sum(gx*im_proj) / sum_proj;
    centroid(2,1) = sum(im_proj*gy) / sum_proj;
    for r=2:numframe
        ht = tic;
        fprintf('calc movement of frame # %d from former frame ... ',r);
        im = getImageMiji(himp_orig,[],channel_use,1);
        im_proj = double(max(im,[],3));
        sum_proj = sum(im_proj(:));
        centroid(1,r) = sum(gx*im_proj) / sum_proj;
        centroid(2,r) = sum(im_proj*gy) / sum_proj;
        fprintf('%g (sec) \n',toc(ht));
    end
    output2([3,4],:,:) = -diff(centroid,[],2);
    cumshift = bsxfun(@plus,cum,cat(3,zeros(2,1,1),cumsum(output2(3:4,:,:),3)));
    
else % image alignment by minimum spanning tree
    
    %%% make mini fft image
    disp('FFT and make low resolution images...');
    tic;
    im_fft_mini = zeros(numxm*2*numym*2,numframe);
    for ct=1:numframe
        im = getImageMiji(himp_orig,[],channel_use,ct);
        im_proj = double(max(im,[],3)); % assume z-stack is already aligned
        im_fft_tmp = fft2(im_proj);
        im_fft_mini(:,ct) = im_fft_tmp(cmat);
    end
    im_fft_mini = reshape(im_fft_mini,[numxm*2,numym*2,numframe]);
    toc;
    
    %%% dft registration by mini images
    disp('Obtain image similarity by dftregistration using low resolution images...')
    tic;
    dftres_mini = zeros(numframe,numframe,4);
    parfor ct1 = 2:numframe
        % for ct1 = 2:numframe
        ht = tic;
        dftres_tmp = zeros(numframe,4);
        for ct2 = 1:ct1-1
            dftres_tmp(ct2,:) = ...
                dftregistration(im_fft_mini(:,:,ct1),im_fft_mini(:,:,ct2)); %#ok<PFBNS>
        end
        dftres_mini(ct1,:,:) = dftres_tmp;
        disp(['Asymmetric loop of dftregistration for frame #',num2str(ct1),...
            ' ...  ',num2str(toc(ht)),' (sec)']);
    end
    disttable = dftres_mini(:,:,1) + dftres_mini(:,:,1)' + diag(inf*ones(1,numframe));
    toc;
    
    %%% calc minimum spanning tree (MST)
    %%% MST minimizes sum of edge weights when the tree covers all the nodes
    disp('Building minimum spanning tree based on similarity...');
    tic;
    g = graph(disttable);
    g.Nodes.Name = cellfun(@num2str,num2cell(1:numframe)','UniformOutput',false);
    t = minspantree(g);
    idx_ecc = 1;
    edgetonew = t.bfsearch(idx_ecc,'edgetonew'); % edges in table format
    numedges = size(edgetonew,1);
    % figure;
    % plot(t,'NodeLabel',t.Nodes.Name,'Layout','layered','Source',num2str(idx_ecc));
    toc;
    
    %%% determine shift amount by 2nd dft registration
    disp('Determine shift amount by 2nd dftregistration using high resolution images...');
    tic;
    dftres = zeros(numedges,4);
    parfor c=1:numedges
        dftres(c,:) = ...
            dftfun(path_im_orig,channel_use,edgetonew(c,:),cum,fiji_directory);
    end
    cumshift = zeros(numframe,2);
    for c=1:numedges
        cumshift(edgetonew(c,2),:) = cumshift(edgetonew(c,1),:)+dftres(c,3:4);
    end
    toc;
end

cum2 = bsxfun(@plus,permute(cumshift,[2,3,1]),cum);
userdata.cum2 = cum2;
setUserData(userdata);

enableGUI(true);

fprintf('Done: Calculation alignment; ');
toc(ht2);

end


function dftres = dftfun(path_im_orig,channel_use,tvec,cum,fiji_directory)
disp(['treefit (fitfunc): ',num2str(tvec(1)),' to ',num2str(tvec(2))]);
Miji_mod(false,fiji_directory); % setup Miji for parfor environment

[~,himp_orig,mfivs_orig] = getImageMiji(path_im_orig,[],[],1);
mfivs_orig.setParallelShift(cum);
mfivs_orig.setAlign(true);
im = getImageMiji(himp_orig,[],channel_use,tvec);

im_proj = double(max(im,[],3)); % assume aligned
im_fft = fft2(im_proj);
dftres = dftregistration(im_fft(:,:,:,:,1),im_fft(:,:,:,:,2),100);
end


function removeROIs()
userdata = getUserData();
hmtm = userdata.hmtm;
strdisp = 'Selected ROIs';
modelRows = getSelectedRowIndexInFigureAsModel();
if numel(modelRows)>0
    if numel(modelRows)==1
        strdisp = ['ROI (',char(hmtm.getValuesAt(modelRows(1),'Name')),')'];
    elseif numel(modelRows)>1
        answer = questdlg('Do you want to remove multiple ROIs?',...
            'remove multiple ROIs?','Yes','No','No');
        drawnow();
        if strcmpi(answer,'No'); return; end % do nothing
    end
    hmtm.removeRowsAt(modelRows,[true,true]);
    userdata.jtable.resort(); % for repaint
    disp([strdisp,' were removed.']);
    clearSelection();
end
end


function setNsigma()
userdata = getUserData();
nsigma = str2double(get(findobj('Tag','edit_option_visualize'),'String'));
userdata.hmrm.setNsigma(nsigma); % including invalidateAll implicitly
end


function setCacheSize()
userdata = getUserData();
sizCache = str2double(get(findobj('Tag','edit_option_cache'),'String'));
userdata.mfivs.clearAndResizeCache(sizCache); % including invalidateAll implicitly
end


%%%%%% analyze function %%%%%%
function [tcrs,frameUnique,namesUnique,tvec,strDesc,leaforder] = getTcrs()
userdata = getUserData();
hmtm = userdata.hmtm;

% convert to matrix format
idxFull = 0:hmtm.getRowCount()-1;
frame = hmtm.getValuesAt(idxFull,'frame');
names = cell(hmtm.getValuesAt(idxFull,'Name'));
[frameUnique,~,idxFrameUnique] = unique(frame);
[namesUnique,~,idxNamesUnique] = unique(names);
numFrameUnique = numel(frameUnique);
numRoiUnique = numel(namesUnique);
idxInMatrix = sub2ind([numFrameUnique,numRoiUnique],idxFrameUnique,idxNamesUnique);

tcrs = [];
strDesc = {};
leaforder = [];

% filtering
namesSelected = userdata.hmrm.getSelectedRowNames();
flagFilter = false(numRoiUnique,1);
if isempty(namesSelected)
    flagFilter = ~flagFilter; % use all rois
else
    for p=1:numel(namesSelected)
        flagFilter = flagFilter | strcmp(namesSelected(p),namesUnique);
    end
end
namesUnique = namesUnique(flagFilter); % filtering


% obtain nomral target
hObjNorm = findobj('Tag','listbox_analyze_plot');
strColNorm = hObjNorm.String(hObjNorm.Value());
if ~isempty(strColNorm)
    idxColNorm = zeros(size(strColNorm));
    strIdentifier = cell(hmtm.getColumnIdentifiers().toArray());
    numCol = numel(strColNorm);
    for p=1:numCol
        idxColNorm(p) = find(strcmp(strColNorm{p},strIdentifier));
    end
    dataNorm = hmtm.getValuesAt(idxFull,idxColNorm-1);
    tcrs = nan(numFrameUnique,numRoiUnique,numCol);
    for p=1:numCol
        tmpmat = nan(numFrameUnique,numRoiUnique);
        tmpmat(idxInMatrix) = dataNorm(:,p);
        tcrs(:,:,p) = tmpmat;
    end
    strDesc = strColNorm;
    tcrs = tcrs(:,flagFilter,:); % filering
end



% obtain numerator
hObjNum = findobj('Tag','popup_analyze_numerator');
strColNum = hObjNum.String{hObjNum.Value};
if ~strcmp(strColNum,'(none)')
    dataNum = hmtm.getValuesAt(idxFull,strColNum);
    strDescNum = strColNum;
    
    % obtain denominator
    hObjDen = findobj('Tag','popup_analyze_denominator');
    strColDen = hObjDen.String{hObjDen.Value};
    if ~strcmp(strColDen,'(none)')
        dataNum = dataNum./hmtm.getValuesAt(idxFull,strColDen);
        strDescNum = [strDescNum,'/',strColDen];
    end
    
    % convert to matrix
    tcrsNum = nan(numFrameUnique,numRoiUnique);
    tcrsNum(idxInMatrix) = dataNum;
    tcrsNum = tcrsNum(:,flagFilter,:); % filering
    
    % Savitzky-Golay filtering
    if get(findobj('Tag','check_analyze_smoothing'),'Value')==1
        sgwindow = str2double(get(findobj('Tag','edit_analyze_smoothing'),'String'));
        %         tcrsNum = sgolayfilt(tcrsNum,2,sgwindow);
        tcrsNum = sgolayfilt(medfilt1(tcrsNum,11),2,sgwindow);
        strDescNum = [strDescNum,', smoothed'];
    end
    
    % normalize
    if get(findobj('Tag','check_analyze_normalize'),'Value')==1
        normalize = @(x)bsxfun(@rdivide,bsxfun(@minus,x,min(x)),max(x)-min(x));
        tcrsNum = normalize(tcrsNum);
        strDescNum = [strDescNum,', normalized'];
    end
    
    % clustering
    if size(tcrsNum,2)>1 && get(findobj('Tag','check_analyze_clustering'),'Value')==1
        hObjMetric = findobj('Tag','popup_analyze_metric');
        strClustMetric = hObjMetric.String{hObjMetric.Value};
        flag_abs = strncmp('abs_',strClustMetric,4);
        hObjMethod = findobj('Tag','popup_analyze_method');
        strClustMethod = hObjMethod.String{hObjMethod.Value};
        strClustCriteria = 'adjacent';
        strClustTransform = 'linear';
        tcrsNoNan = tcrsNum;
        tcrsNoNan(isnan(tcrsNoNan)) = 0;
        if flag_abs
            distance = 1-abs(1-pdist(tcrsNoNan',strClustMetric(5:end)));
        else
            distance = pdist(tcrsNoNan',strClustMetric);
        end
        distance(~isfinite(distance)) = inf;
        linkmat = linkage(distance,strClustMethod);
        leaforder = optimalleaforder(linkmat,distance,'criteria',strClustCriteria,...
            'transformation',strClustTransform);
        strDescNum = [strDescNum,', ',strClustMetric,', ',strClustMethod];
    end
    
    % merge
    tcrs = cat(3,tcrs,tcrsNum);
    strDesc = cat(1,strDesc,{strDescNum});
end

% time vector
vps = str2double(get(findobj('Tag','edit_analyze_vps'),'String'));
tvec = (frameUnique-1)/vps;

% leaforder
if isempty(leaforder)
    leaforder = 1:size(tcrs,2);
end

end


function createHeatmap()
userdata = getUserData();
[tcrs,frameUnique,namesUnique,tvec,strDesc,leaforder] = getTcrs();
numRoiUnique = numel(namesUnique);

% create figure
hf = figure(...
    'paperorientation','portrait',...
    'paperUnits','normalized',...
    'Paperposition',[0,0,1,1],...
    'Windowstyle','docked',...
    'ColorMap',jet(256));

strTitle = ['Heatmap of ',strDesc{end}];

% whether show target movement
if get(findobj('Tag','check_analyze_show_move'),'Value')==1
    % plot movement
    ax1 = subplot(10,1,1);
    plot(tvec,permute(userdata.cum2(1,1,frameUnique),[3,2,1]));
    axis tight;
    set(ax1,'XTick','');
    ylabel('Movement');
    title(strTitle,'interpreter','none');
    
    % show heatmap
    ax2 = subplot(10,1,2:10);
    imagesc(tvec([1,end]),[1,size(tcrs,2)],tcrs(:,leaforder,end)');
    colorbar('Location','southoutside');
    drawnow();
    
    % align two axes
    pos1 = get(ax1,'Position');
    pos2 = get(ax2,'Position');
    pos1([1,3]) = pos2([1,3]);
    pos1(4) = pos1(2)+pos1(4) - (pos2(2)+pos2(4));
    pos1(2) = pos2(2)+pos2(4);
    set(ax1,'position',pos1);
    linkaxes([ax1,ax2],'x');
else
    % show heatmap
    imagesc(tvec([1,end]),[1,size(tcrs,2)],tcrs(:,leaforder,end)');
    colorbar;
    title(strTitle,'interpreter','none');
end

strYTickLabel = namesUnique(leaforder);
for p=1:2:numel(strYTickLabel)
    strYTickLabel(p) = {[strYTickLabel{p},'--------']};
end
xlabel('Time (sec)');
set(gca,'YTick',1:numRoiUnique,'YTickLabel',strYTickLabel,'TickLabelInterpreter','none');

strPathFigSave = fullfile(userdata.PathName,['heatmap_',userdata.FileName(1:end-4),'.pdf']);
print(hf,'-dpdf',strPathFigSave);

strPathTableSave = fullfile(userdata.PathName,['heatmap_',userdata.FileName]);
clear userdata hmtm hf ax* hObj*;
save(strPathTableSave);

end


function plotTcrs()
userdata = getUserData();
strPathName = userdata.PathName;
strFileName = userdata.FileName(1:end-4);
maxCol = 5;

%%% figure file save path
strPathMatSave = fullfile(strPathName,['plot_',strFileName,'.mat']);
strPathFigSave = fullfile(strPathName,['plot_',strFileName,'.ps']);
if exist(strPathFigSave,'file')
    delete(strPathFigSave);
end

%%% obtain data
[tcrs,frameUnique,namesUnique,tvec,strDesc,leaforder] = getTcrs();
numRoiUnique = numel(namesUnique);

%%% wheter show move
flagShowMove = get(findobj('Tag','check_analyze_show_move'),'Value')==1;
if flagShowMove
    movement = permute(userdata.cum2(1,1,frameUnique),[3,2,1]);
    tcrs = cat(3,movement(:,ones(1,numRoiUnique)),tcrs);
    strDesc = cat(1,{'movement'},strDesc);
end
numRow = size(tcrs,3);

%%% whether merge plot
flagMerge = get(findobj('Tag','check_analyze_merge'),'Value')==1;
if flagMerge
    numCol = 1;
    numFigure = 1;
else
    numCol = maxCol;
    numFigure = ceil(numRoiUnique/numCol);
end

%%% main plot loop
for cf=1:numFigure
    % create figure
    hf = figure(...
        'paperorientation','landscape',...
        'paperUnits','normalized',...
        'Paperposition',[0,0,1,1],...
        'Windowstyle','docked');
    
    for cc=1:numCol
        idxRoi = (cf-1)*numCol+cc;
        if any(idxRoi > numRoiUnique); break; end
        if flagMerge; idxRoi = 1:numRoiUnique; end
        idxRoiDisp = leaforder(idxRoi);
        for cr=1:numRow
            ax(cc,cr) = subplot(numRow,numCol,numCol*(cr-1)+cc); %#ok<AGROW>
            if cr==1 % for the first row
                %                 if flagShowMove && flagMerge
                %                     plot(ax(cc,cr),tvec,movement);
                %                 else
                plot(ax(cc,cr),tvec,tcrs(:,idxRoiDisp,cr));
                %                 end
                if flagMerge
                    title(strDesc{end},'interpreter','none');
                    legend(namesUnique(idxRoiDisp),'interpreter','none');
                else
                    title(namesUnique{idxRoiDisp},'interpreter','none');
                end
            else
                plot(ax(cc,cr),tvec,tcrs(:,idxRoiDisp,cr));
            end
            if cr==numRow % for the last row
                xlabel('Time (sec)');
            else
                set(ax(cc,cr),'XTickLabel',[]);
            end
            if cc==1 % for the first column
                strYLabel = strsplit(strDesc{cr},',');
                ylabel(strYLabel{1},'interpreter','none');
            end
            xlim([min(tvec),max(tvec)]);
        end
        
        % height settings
        ap1 = get(ax(cc,1),'position');
        ap2 = get(ax(cc,numRow),'position');
        ah = (ap1(2)+ap1(4)-ap2(2))/numRow;
        for cr=1:numRow
            tmppos = ap2;
            tmppos(2) = ap2(2)+ah*(numRow-cr);
            tmppos(4) = ah;
            set(ax(cc,cr),'position',tmppos);
        end
        linkaxes(ax(cc,:),'x');
        
    end
    print(hf,'-dpsc','-append',strPathFigSave);
end

clear userdata hmtm hf ax* hObj*;
save(strPathMatSave);

end


function initAnalyze() % call once after the data table is loaded
userdata = getUserData();
hmtm = userdata.hmtm;
strColumnNames = cell(hmtm.getColumnIdentifiers().toArray());
strPopupNames = {'(none)'};
for p=1:hmtm.getColumnCount()
    if isnumeric(hmtm.getValuesAt(0,p-1))
        strPopupNames=cat(1,strPopupNames,strColumnNames(p));
    end
end

set(findobj('Tag','listbox_analyze_plot'),'String',strPopupNames(2:end));
set(findobj('Tag','popup_analyze_numerator'),'String',strPopupNames);
set(findobj('Tag','popup_analyze_denominator'),'String',strPopupNames);

end


%%%%%% visualize tracking result %%%%%%
function [tcrs,tcrsDiff,frameUnique,namesUnique,idxInMatrix] = getTcrsMu()

userdata = getUserData();
hmtm = userdata.hmtm;

% convert to matrix format
idxFull = 0:hmtm.getRowCount()-1;
frame = hmtm.getValuesAt(idxFull,'frame');
names = cell(hmtm.getValuesAt(idxFull,'Name'));
[frameUnique,~,idxFrameUnique] = unique(frame);
[namesUnique,~,idxNamesUnique] = unique(names);
numFrameUnique = numel(frameUnique);
numRoiUnique = numel(namesUnique);
idxInMatrix = sub2ind([numFrameUnique,numRoiUnique],idxFrameUnique,idxNamesUnique);

tcrs = nan(numFrameUnique,numRoiUnique,3);
flagOrigState = hmtm.flagInterpZ;
hmtm.setInterpZ(true);
tmpmu = hmtm.getValuesAt(idxFull,'mu');
hmtm.setInterpZ(flagOrigState);
for p=1:3
    tmpmat = nan(numFrameUnique,numRoiUnique);
    tmpmat(idxInMatrix) = tmpmu(:,p);
    tcrs(:,:,p) = tmpmat;
end
tcrsDiff = cat(1,nan(1,numRoiUnique,3),diff(tcrs,1,1));
tcrsDiff(cat(1,true,diff(frameUnique)~=1),:,:) = nan;

end


function updateIncoherency(tcrs,tcrsDiff,idxInMatrix)
userdata = getUserData();
hmtm = userdata.hmtm;
idxFull = 0:hmtm.getRowCount()-1;

% calc difference and incoherency
[numFrameUnique,numRoiUnique,~] = size(tcrs);
incoherency = nan(numFrameUnique,numRoiUnique);
for ct=2:numFrameUnique
    dist = pdist(permute(tcrs(ct,:,:),[2,3,1]));
    weight = exp(-(dist*0.01).^2);
    weight(isnan(weight)) = 0; % for lacked roi
    dist_move = pdist(permute(tcrsDiff(ct,:,:),[2,3,1]));
    dist_move(isnan(dist_move)) = 0; % for lacked roi
    incoherency(ct,:) = sum(squareform(weight.*dist_move.^2),2);
end

retmat = zeros(numel(idxInMatrix),4);
for p=1:3
    tmpmat = tcrsDiff(:,:,p);
    retmat(:,p) = tmpmat(idxInMatrix);
end
retmat(:,4) = incoherency(idxInMatrix);
idxcol = cat(2,hmtm.getColumnIndex('xd'),hmtm.getColumnIndex('yd'),...
    hmtm.getColumnIndex('zd'),hmtm.getColumnIndex('Incoherency'));
hmtm.setValuesAt(retmat,idxFull,idxcol,[false,true,false]);

end


function visTrackInit()
enableGUI(false);

[tcrs,tcrsDiff,frameUnique,namesUnique,idxInMatrix] = getTcrsMu();
updateIncoherency(tcrs,tcrsDiff,idxInMatrix); % update table

numRoiUnique = numel(namesUnique);
minframe = min(frameUnique);
maxframe = max(frameUnique);
sliderstep = 1/(maxframe-minframe);
strLabels = strcat(' ',namesUnique);
limUIrange = permute(cat(1,min(min(tcrs)),max(max(tcrs))),[1,3,2]);
frameRange = [min(frameUnique),max(frameUnique)];

%%% initialize control panel
hp = findobj('Tag','tab_vistrack');
if isempty(findobj('Tag','slider_vistrack_t'))
    
    visTrackCreateControlsDouble(hp,0.79,0.04,'x',limUIrange(:,1),@visTrackSetRange);
    visTrackCreateControlsDouble(hp,0.64,0.04,'y',limUIrange(:,2),@visTrackSetRange);
    visTrackCreateControlsDouble(hp,0.49,0.04,'z',limUIrange(:,3),@visTrackSetRange);
    visTrackCreateControlsSingle(hp,0.38,0.04,'t',frameRange([1,1,2]),@visTrackDraw);
    visTrackCreateControlsSingle(hp,0.27,0.04,'Marker',[eps,10,50],@visTrackSetSizeMarker);
    visTrackCreateControlsSingle(hp,0.16,0.04,'Label',[eps,6,50],@visTrackSetSizeLabel);
    
    fun_checkbox(hp,'Show Color','check_vistrack_color',[0.05,0.09,0.45,0.04],...
        'Value',1,'Callback',@visTrackColor);
    fun_checkbox(hp,'Show Label','check_vistrack_label',[0.50,0.09,0.45,0.04],...
        'Value',1,'Callback',@visTrackDraw);
    fun_label(hp,'Selection','label_vistrack_select',[0.05,0.02,0.35,0.04]);
    fun_edit( hp,'',         'edit_vistrack_select', [0.40,0.02,0.55,0.04],...
        'Callback',@visTrackSelect);
    
    set(findobj('Tag','slider_vistrack_t'),'SliderStep',[1,10]*sliderstep);
    
end

hfig = findobj('Tag','figure_vistrack');
if isempty(hfig) % initialize view panel
    hfig = figure('Tag','figure_vistrack');
end

userdata = visTrackGetUserData();
if isfield(userdata,'haxis')
    haxis = userdata.haxis;
else
    haxis = axes('Parent',hfig);
end

%%% create plot template
cla(haxis);
hold on;
stropt = {'Visible','off','Clipping','on','Parent',haxis};
for p=1:numRoiUnique
    userdata.hroi(p) = line(0,0,0,'Marker','o','DisplayName',strLabels{p},stropt{:});
    userdata.hmove(p) = line([0,1],[0,1],[0,1],stropt{:});
    userdata.hlabel(p) = text(0,0,0,strLabels{p},stropt{:});
end


%%% set mouse motion listener
hobj = findobj('-regexp','Tag','slider_vistrack_*');
hlistenerSlider = [];
for p=1:numel(hobj)
    hlistenerSlider = [hlistenerSlider,...
        addlistener(hobj(p),'ContinuousValueChange',hobj(p).Callback)]; %#ok<AGROW>
end

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(handle(hfig),'JavaFrame');
hlistenerAxis = addlistener(jFrame.fHG2Client.getAxisComponent,'MouseWheelMoved',@visTrackWheel);
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');


%%% save data and handles as userdata of vistrack figure
userdata.tcrs = tcrs;
userdata.tcrsDiff = tcrsDiff;
userdata.frameUnique = frameUnique;
userdata.namesUnique = namesUnique;
userdata.hlistenerAxis = hlistenerAxis;
userdata.hlistenerSlider = hlistenerSlider;
userdata.haxis = haxis;
visTrackSetUserData(userdata);

visTrackColor();
visTrackDraw();

% title(strFileName,'interpreter','none');
xlabel('x');
ylabel('y');
zlabel('z');

figure(hfig);
cameratoolbar(hfig,'Show');
cameratoolbar(hfig,'SetMode','orbit');
cameratoolbar(hfig,'SetCoordSys','none');
camproj(haxis,'perspective')
box(haxis,'on');
grid(haxis,'on');
hold(haxis,'on');
axis(haxis,'equal');
hold(haxis,'on');
set(haxis,'Xlim',limUIrange(:,1),'Ylim',limUIrange(:,2),'Zlim',limUIrange(:,3));

enableGUI(true);
end


function visTrackCreateControlsSingle(parent,start,height,label,range,fun)
fun_label(parent,label,['label_vistrack_',label],[0.05,start+height,0.45,height]);
fun_edit(parent,range(2),['edit_vistrack_',label],[0.50,start+height,0.45,height],...
    'CallBack',fun);
fun_slider(parent,['slider_vistrack_',label],[0.05,start,0.90,height],...
    range(1),range(3),range(2),'CallBack',fun);
end


function visTrackCreateControlsDouble(parent,start,height,label,range,fun)
fun_label(parent,label,['label_vistrack_',label],[0.05,start+height*2,0.2,height]);
fun_edit(parent,num2str(range(1)),['edit_vistrack_',label,'_min'],...
    [0.25,start+height*2,0.35,height],'CallBack',fun);
fun_edit(parent,num2str(range(2)),['edit_vistrack_',label,'_max'],...
    [0.60,start+height*2,0.35,height],'CallBack',fun);
fun_slider(parent,['slider_vistrack_',label,'_min'],[0.05,start+height,0.90,height],...
    range(1),range(2),range(1),'CallBack',fun);
fun_slider(parent,['slider_vistrack_',label,'_max'],[0.05,start,0.90,height],...
    range(1),range(2),range(2),'CallBack',fun);
end


function visTrackDraw(obj,~,~)

if exist('obj','var') && strcmp(obj.Style,'edit')
    t = str2double(obj.String);
    set(findobj('Tag','slider_vistrack_t'),'Value',t); % update
else % from slider or axis callback
    t = round(get(findobj('Tag','slider_vistrack_t'),'Value'));
    set(findobj('Tag','edit_vistrack_t'),'String',num2str(t)); % update
end

userdata = visTrackGetUserData();
tcrs = userdata.tcrs;
tcrsDiff = userdata.tcrsDiff;
idxt = find(userdata.frameUnique==t);
numrois = size(tcrs,2);

% tv = [t,min(t+1,size(tcrs,1))];
for p=1:numrois
    if ~isempty(idxt) && ~isnan(tcrs(idxt,p,1))
        set(userdata.hroi(p),...
            'XData',tcrs(idxt,p,1),...
            'YData',tcrs(idxt,p,2),...
            'ZData',tcrs(idxt,p,3),...
            'Visible','on');
        set(userdata.hmove(p),...
            'XData',[tcrs(idxt,p,1),tcrs(idxt,p,1)-tcrsDiff(idxt,p,1)],...
            'YData',[tcrs(idxt,p,2),tcrs(idxt,p,2)-tcrsDiff(idxt,p,2)],...
            'ZData',[tcrs(idxt,p,3),tcrs(idxt,p,3)-tcrsDiff(idxt,p,3)],...
            'Visible','on');
        set(userdata.hlabel(p),...
            'Position',permute(tcrs(idxt,p,:),[3,1,2]),...
            'Visible','on');
    else
        userdata.hroi(p).Visible = 'off';
        userdata.hmove(p).Visible = 'off';
        userdata.hlabel(p).Visible = 'off';
    end
end

if get(findobj('Tag','check_vistrack_label'),'Value') == 0 % do not display label
    set(userdata.hlabel,'Visible','off');
end

end


function visTrackColor(~,~,~)
userdata = visTrackGetUserData();
if get(findobj('Tag','check_vistrack_color'),'Value')==1
    numrois = size(userdata.tcrs,2);
    clmap = hsv2rgb(cat(2,rand(numrois,1),ones(numrois,2))); % hsv
    userdata.clmap = clmap;
    visTrackSetUserData(userdata);
end
visTrackSelect();
end


function visTrackSelect(~,~,~)
userdata = visTrackGetUserData();
strRegexp = ['^',get(findobj('Tag','edit_vistrack_select'),'String'),'$'];
namesUnique = userdata.namesUnique;
clmap = userdata.clmap;
numrois = size(userdata.tcrs,2);

idxRegexp = regexp(namesUnique,strRegexp);
idxSelect = find(cellfun(@numel,idxRegexp));

if isempty(idxSelect) && get(findobj('Tag','check_vistrack_color'),'Value')==1
    idxSelect = 1:numrois;
end

%%% set default colors for all objects
set(userdata.hroi,'Color','k');
set(userdata.hmove,'Color','k');
set(userdata.hlabel,'Color','k');

%%% coloring
for p=idxSelect(:)'
    userdata.hroi(p).Color = clmap(p,:);
    userdata.hmove(p).Color = clmap(p,:);
    userdata.hlabel(p).Color = clmap(p,:);
end

end


function visTrackWheel(~,ev)
userdata = visTrackGetUserData;
haxis = userdata.haxis;
cp = haxis.CameraPosition;
ct = haxis.CameraTarget;
cuv = haxis.CameraUpVector;
cva = haxis.CameraViewAngle;

modifier = ev.getModifiersEx();
mask_ctrl  = java.awt.event.MouseEvent.CTRL_DOWN_MASK;
mask_shift = java.awt.event.MouseEvent.SHIFT_DOWN_MASK;
mask_alt = java.awt.event.MouseEvent.ALT_DOWN_MASK;
flag_ctrl  = bitand(modifier,mask_ctrl )==mask_ctrl;
flag_shift = bitand(modifier,mask_shift)==mask_shift;
flag_alt = bitand(modifier,mask_alt)==mask_alt;
scroll = ev.getPreciseWheelRotation()*ev.getScrollAmount();

if flag_ctrl % horizontally scrolling
    transloc = cross(ct-cp,cuv)*cva*0.0001*scroll;
    haxis.CameraPosition = cp + transloc;
    haxis.CameraTarget   = ct + transloc;
elseif flag_shift % vertically scrolling
    transloc = sqrt(sum((ct-cp).^2))*cuv*cva*0.0001*scroll;
    haxis.CameraPosition = cp + transloc;
    haxis.CameraTarget   = ct + transloc;
elseif flag_alt % temporally scrolling
    hobj = findobj('Tag','slider_vistrack_t');
    tmpvalue = hobj.Value + hobj.SliderStep(1)*(hobj.Max-hobj.Min)*sign(scroll);
    hobj.Value = max(hobj.Min,min(hobj.Max,tmpvalue));
    hobj.Callback();
else % zoom
    haxis.CameraViewAngle = cva*(1+0.01*scroll);
end

end


function visTrackSetRange(obj,~)
userdata = visTrackGetUserData();
hobjSlider = flipud(findobj('-regexp','Tag','slider_vistrack_[xyz]_(min|max)'));
hobjEdit = flipud(findobj('-regexp','Tag','edit_vistrack_[xyz]_(min|max)'));

switch obj.Style
    case {'slider'}
        objlist = hobjSlider;
        values = cell2mat(get(hobjSlider,'Value'));
    case {'edit'}
        objlist = hobjEdit;
        values = str2num(char(get(hobjEdit,'String'))); %#ok<ST2NM>
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
    minvalue = get(hobjSlider(p),'Min');
    maxvalue = get(hobjSlider(p),'Max');
    values(p) = max(min(values(p),maxvalue),minvalue);
end

%%% update components
for p=1:6
    set(hobjSlider(p),'Value',values(p));
    set(hobjEdit(p),'String',num2str(values(p)));
end

%%% update view
set(userdata.haxis,'xlim',values(1:2),'ylim',values(3:4),'zlim',values(5:6));

end


function visTrackSetSizeMarker(obj,~)
userdata = visTrackGetUserData();
switch obj.Style
    case {'slider'}
        value = get(obj,'Value');
    case {'edit'}
        value = str2double(obj.String);
end
hobjSlider = findobj('Tag','slider_vistrack_Marker');
hobjEdit = findobj('Tag','edit_vistrack_Marker');
value = max(min(value,hobjSlider.Max),hobjSlider.Min);
set(hobjSlider,'Value',value);
set(hobjEdit,'String',num2str(value));
set(userdata.hroi,'MarkerSize',value);
end


function visTrackSetSizeLabel(obj,~)
userdata = visTrackGetUserData;
switch obj.Style
    case {'slider'}
        value = get(obj,'Value');
    case {'edit'}
        value = str2double(obj.String);
end
hobjSlider = findobj('Tag','slider_vistrack_Label');
hobjEdit = findobj('Tag','edit_vistrack_Label');
value = max(min(value,hobjSlider.Max),hobjSlider.Min);
set(hobjSlider,'Value',value);
set(hobjEdit,'String',num2str(value));
set(userdata.hlabel,'FontSize',value);
end


function userdata = visTrackGetUserData()
userdata = get(findobj('Tag','figure_vistrack'),'UserData');
end


function visTrackSetUserData(userdata)
set(findobj('Tag','figure_vistrack'),'UserData',userdata);
end



%%%%%% tracking using qpipe %%%%%%
function trackqp_setworkdir()
tmpdir = uigetdir();
if tmpdir~=0
    set(findobj('Tag','edit_trackqp_workdir'),'String',tmpdir);
end
end


function trackqp_export()
%%% export files required by qpipe tracking
run_id = get(findobj('Tag','edit_trackqp_runid'),'String');
frames = eval(get(findobj('Tag','edit_trackqp_frames'),'String'));
if ~iscell(frames); frames={frames}; end
strTargetDir = get(findobj('Tag','edit_trackqp_workdir'),'String');

disp('Exporting files required by QPipe tracking...');
enableGUI(false);

%%% fixed parameters
strSummaryFile = 'summary.csv';
strPifDir = 'doAff';
strPosDir = 'uniqPos';
strPosFile = 'unique_positions_0000.txt';

%%% obtain data
userdata = getUserData();
if ~exist('strTargetDir','var')||isempty(strTargetDir)
    strTargetDir = userdata.PathName;
end
himp = userdata.himp;
hmtm = userdata.hmtm;
dims = himp.getDimensions()';
idxc = himp.getChannel();

%%% create output
strSummaryPath = fullfile(strTargetDir,run_id,strSummaryFile);
disp(['Writing summary file : ',strSummaryPath]);
if ~exist(fileparts(strSummaryPath),'dir'); mkdir(fileparts(strSummaryPath)); end
fid = fopen(strSummaryPath,'w');
nsubs = numel(frames);
for p=1:nsubs
    
    %%% prepare pif headers
    PIF_SIGNATURE  = int8(sprintf('%s%c', 'PIF', 0));
    PIF_OFFSET     = uint32(48);
    PIF_PIXEL_TYPE = uint32(2); % image will be converted to uint16;
    PIF_PIXEL_SIZE = uint32(2);
    dims_sub = uint32([dims([1,2,4]), 0, numel(frames{p})]);
    header1 = PIF_SIGNATURE;
    header2 = [PIF_OFFSET, dims_sub, PIF_PIXEL_TYPE, PIF_PIXEL_SIZE, ...
        zeros(1, PIF_OFFSET/4 - 9)];
    
    %%% write pif file
    imname = strsplit(userdata.imname,'_');
    sub_id = num2str(p,'%04d');
    data_id = [imname{1},'_sub',sub_id];
    strPifName = sprintf('%s_mCherry_%dx%d_z%d_t%d-%d.pif',data_id,dims_sub);
    strPifPath = fullfile(strTargetDir,run_id,sub_id,strPifDir,strPifName);
    disp(['Writing pif file : ',strPifPath]);
    if ~exist(fileparts(strPifPath),'dir'); mkdir(fileparts(strPifPath)); end
    fd = fopen(strPifPath,'w');
    fwrite(fd, header1, 'uint8');
    fwrite(fd, header2, 'uint32');
    for q=1:numel(frames{p})
        fwrite(fd,uint16(getImageMiji(himp,[],idxc,frames{p}(q))),'uint16');
    end
    fclose(fd);
    
    %%% write uniquepos.txt
    %     pos = hmtm.getValuesInFrame(frames{p}(1),'mu');
    pos = hmtm.getValuesInFrame(frames{p}(1),'mu') - 1; % for 0-based index
    pos(:,4) = 0.1; % add required elements;; not important?
    strPosPath = fullfile(strTargetDir,run_id,sub_id,strPosDir,strPosFile);
    if ~exist(fileparts(strPosPath),'dir'); mkdir(fileparts(strPosPath)); end
    dlmwrite(strPosPath,size(pos,1));
    dlmwrite(strPosPath,pos,'delimiter',' ','-append');
    
    %%% write summary.csv
    fprintf(fid,'%s,%s,%d,%d,%d,%d,%d\n',sub_id,data_id,dims_sub);
    
end
fclose(fid);

disp('Exporting files required by QPipe tracking... done!');
enableGUI(true);

end


function trackqp_import(strin)
%%% import qpipe tracking results

run_id = get(findobj('Tag','edit_trackqp_runid'),'String');
frames = eval(get(findobj('Tag','edit_trackqp_frames'),'String'));
if ~iscell(frames); frames={frames}; end
strTargetDir = get(findobj('Tag','edit_trackqp_workdir'),'String');

disp('Importing QPipe tracking results...');
switch strin
    case 'spf'
        disp('Import SPF results...');
        strResult = 'spf_%04d_position.txt';
    case 'modetrack01'
        disp('Import MODETRACK01 results...');
        strResult = 'MODETRACK01_%04d_position0.txt';
    otherwise
        disp('Error: strin is unusual.');
        return;
end

enableGUI(false);

userdata = getUserData();
if ~exist('strTargetDir','var')||isempty(strTargetDir)
    strTargetDir = userdata.PathName;
end
hmtm = userdata.hmtm;
nsubs = numel(frames);
flagFrameUsed = false(max(cat(2,frames{:})),1);
numcell = numel(hmtm.convertFrameToRows(frames{1}(1)));
dv = userdata.defaultdata(ones(numcell*numel(flagFrameUsed),1),:);
cr = 0;
for p=1:nsubs
    numt = numel(frames{p});
    numcell = numel(hmtm.convertFrameToRows(frames{p}(1)));
    dd = userdata.defaultdata(ones(numcell,1),:);
    dd(:,hmtm.getColumnIndex('Name')+1) = cell(hmtm.getValuesInFrame(frames{p}(1),'Name'));
    dd(:,hmtm.getColumnIndex('params_em')+1) ...
        = num2cell(hmtm.getValuesInFrame(frames{p}(1),'params_em'));
    
    sub_id = num2str(p,'%04d');
    strSubDir = fullfile(strTargetDir,run_id,sub_id);
    disp(['Reading from ',strSubDir]);
    strCentAveDir = fullfile(strSubDir,'centAve');
    strTrackDir = fullfile(strSubDir,'trackPos');
    
    centave = load(fullfile(strCentAveDir,['shiftsAve_t0^%',num2str(numt),'.txt']));
    for q=1:numt
        if flagFrameUsed(frames{p}(q)); continue; end % skip overwrapped frame
        flagFrameUsed(frames{p}(q)) = true;
        %         tracked = load(fullfile(strTrackDir,sprintf(strResult,q-1)));
        tracked = load(fullfile(strTrackDir,sprintf(strResult,q-1))) + 1; % convert from 0-based to 1-based
        dd(:,hmtm.getColumnIndex('mu')+1) ...
            = num2cell(bsxfun(@minus,tracked(:,1:3),centave(q,:)));
        dd(:,hmtm.getColumnIndex('frame')+1) = {frames{p}(q)};
        if (numcell+cr)<=size(dv,1)
            dv((1:numcell)+cr,:) = dd;
        else
            dv = cat(1,dv(1:cr,:),dd);
        end
        cr = cr + numcell;
    end
end
dv = dv(1:cr,:);
dv(:,hmtm.getColumnIndex('uniqueID')+1) = num2cell(1:size(dv,1));

%%% set imported data to RoiManager
disp('Update ROI information based on imported data.');
hmtm.removeRowsAt(0:hmtm.getRowCount()-1,[false,true,false]);
hmtm.insertRowsAt(dv,0,[false,false,false]);

hmtm.convertCharToString(); % call directly
userdata.frame.repaint(); % update roi manager display
hmrm = userdata.hmrm;
hmrm.invalidateAll(); % force update roi display % call directly
autosave(); % call directly
clearSelection();

disp('Importing QPipe tracking results... done!');
enableGUI(true);

end
