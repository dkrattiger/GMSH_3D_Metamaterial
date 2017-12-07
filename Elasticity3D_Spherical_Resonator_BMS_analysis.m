%% 3D_Elasticity_FE_Code
% ======================================================================= %
% Dimitri Krattiger
% 10-12-2012

close all; clear; clc; 
format compact
format shortg
warning off all
tic

addpath(genpath('Global_Files'))

%% Check for & Create Save directories
% ======================================================================= %

save_results = false;
load_results = true;

% if save_data folder does not exist, make it
if save_results
    
    % make save_data directory
    if ~exist('save_data','dir')
        mkdir('save_data');
    end
    
    % make save_data/models subdirectory
    if ~exist('save_data/models','dir')
        mkdir('save_data/models');
    end
    
    % make save_data/solutions subdirectory
    if ~exist('save_data/solutions','dir')
        mkdir('save_data/solutions');
    end
    
    % make save_data/BMS_solutions subdirectory
    if ~exist('save_data/BMS_solutions','dir')
        mkdir('save_data/BMS_solutions');
    end
    
    % make figures directory
    if ~exist('figures','dir')
        mkdir('figures');
    end
    
end

%% User Selected Parameters
% ======================================================================= %

% Solution Options
n_curves = 150;         % number of dispersion curves to plot
full_soln = true;       % compute full dispersion
BMS_LIR_soln = true;    % compute LIR-BMS dispersion
use_v0 = false;         % specify starting vector for eigenvalue computation

% % save results
% save_results = true;
% date_string = datestr(now,'dd-mmm-yyyy HH-MM');

% Mode Animations
make_vid = false;
show_stress = true;
mode_plot = [1:20];     % select modes to plot. [] =  don't animate any modes
% mode_plot = [];    % select modes to plot. [] =  don't animate any modes
k_sel = 16;          % k-point to plot modes at

date_string = datestr(now,'dd-mmm-yyyy HH-MM');


%% Model Definition
% ======================================================================= %
% Define Unit Cell Geometry
r_scube = 2.5e-3;               % small cube radius
r_bcube = 7.75e-3;              % big cube radius
r_spher = 5e-3;                 % internal sphere radius
t_layer = 2.5e-3;               % sphere layer thickness
% t_layer = 0;                  % sphere layer outer radius

r_layer = r_spher + t_layer;    % sphere layer outer radius

Lx = 2*r_bcube;
Ly = 2*r_bcube;
Lz = 2*r_bcube;

r1 = [Lx;0;0];
r2 = [0;Ly;0]; 
r3 = [0;0;Lz];
R = [r1,r2,r3];

%% Define Mesh Size
% ======================================================================= %

% n = 2 (linear)
% =====
% md = 1; n = 2; % mesh density (1257 dof, 1074 periodic dof)
% md = 2; n = 2;% mesh density (9135 dof, 8484 periodic dof)
% md = 3; n = 2; % mesh density (29973 dof, 28566 periodic dof)
% md = 4; n = 2; % mesh density (70107 dof, 67656 periodic dof)

% n = 3 (quadratic)
% =====
md = 1; n = 3; % mesh density (9135 dof, 8484 periodic dof)
% md = 2; n = 3;% mesh density (70107 dof, 67656 periodic dof)

% n = 4 (cubic)
% =====
% md = 1; n = 4; % mesh density (29973 dof, 28566 periodic dof)

% n = 5 (quartic)
% =====
%md = 1; n = 5; % mesh density (72123 dof, 69546 periodic dof)

% number of divisions per region edge
n_ele_sc = 4*md;     % number of elements per edge of internal cube
n_ele_sph = md;      % number of radial elements from internal cube to sphere edge
n_ele_layer = md;    % number of radial elements in spherical layer
n_ele_bc = md;       % number of radial elements from layer to outer cube

%% Lead Properties
% ======================================================================= %
rho1    = 11.6e3;  % density
% E1      = 200e9;     % Young's Modulus
% v1      = 0.33;     % Poisson's Ratio

% Lame's Parameters
lam1    = 4.23e10;
mu1     = 1.49e10;

% Isotropic Elasticity Tensor
D1 = [lam1+2*mu1, lam1, lam1, 0, 0, 0; 
      lam1, lam1+2*mu1, lam1, 0, 0, 0; 
      lam1, lam1, lam1+2*mu1, 0, 0, 0; 
      0, 0, 0, mu1, 0, 0; 
      0, 0, 0, 0, mu1, 0; 
      0, 0, 0, 0, 0, mu1];

%% epoxy Properties
% ======================================================================= %
rho2    = 1.18e3;  	% density

% Lame's Parameteres  
kappa2 = 7.61e9;
mu2     = 1.59e9;
lam2 = kappa2-2*mu2;

% Isotropic Elasticity Tensor
D2 = [lam2+2*mu2, lam2, lam2, 0, 0, 0;
      lam2, lam2+2*mu2, lam2, 0, 0, 0; 
      lam2, lam2, lam2+2*mu2, 0, 0, 0; 
      0, 0, 0, mu2, 0, 0; 
      0, 0, 0, 0, mu2, 0; 
      0, 0, 0, 0, 0, mu2]; 
  
%% Rubber Properties
% ======================================================================= %
rho3    = 1.3e3;  	% density kg/m^3

% Lame's Parameteres
lam3    = 6e5;      
mu3     = 4e4;

% Isotropic Elasticity Tensor
D3 = [lam3+2*mu3, lam3, lam3, 0, 0, 0;
      lam3, lam3+2*mu3, lam3, 0, 0, 0; 
      lam3, lam3, lam3+2*mu3, 0, 0, 0; 
      0, 0, 0, mu3, 0, 0; 
      0, 0, 0, 0, mu3, 0; 
      0, 0, 0, 0, 0, mu3]; 
  
% % lower rubber mass
% rho3 = rho3/1e3;
  
  
%% store material properties in cell arrays;
% ======================================================================= %
Ds = {D1,D2,D3};
rhos = {rho1,rho2,rho3};

%% Create Wave vector path
% ======================================================================= %
% n_bz = 1;             % number of brillouin Zones to cover

% "Corners" in wave vector path
% sym_pts = {'Gamma','X','M','R','Gamma'};
sym_pts = {'M_10','Gamma','X_5'};

% integer refinement level of wave-vector path
kap_ref = 5;    

% create set of wave-vectors at which to compute frequencies
[kappa,kappa_plot] = wave_vector(sym_pts,kap_ref,R);
n_kap = size(kappa,2);


% secondary kappa vector
kap_ref2 = 9;    

% create set of wave-vectors at which to compute frequencies
[kappa2,kappa_plot2] = wave_vector(sym_pts,kap_ref2,R);
n_kap2 = size(kappa2,2);

%% define 3d mesh
% ======================================================================= %

[xlocs,ylocs,zlocs,element_node_table,pattern_vec,patchfaces_lite,C_lite,fedges] = ...
     mesh_spherical_resonator(n,n_ele_sc,n_ele_sph,n_ele_layer,...
         n_ele_bc,r_scube,r_bcube,r_spher,r_layer);
     
n_dpn = 3;                  % number of DOFs per node
n_nodes = length(xlocs);    % number of nodes
n_dof = n_nodes*n_dpn;      % number of total DOFs
n_elements = length(element_node_table);    % number of elements

%% Find node sets for boundaries
% ======================================================================= %

% collect unit cell coordinates
coordinates = [xlocs,ylocs,zlocs];

% compute node sets for overall unit cell
[node_sets] = find_node_sets(coordinates,R);
dof_sets = node2dof(node_sets,n_dpn);


%% Display unit cell geometry
% ======================================================================= %

figure(2);clf
plot_FEM_model(coordinates,patchfaces_lite,C_lite,fedges);
axis equal
view(3)
drawnow
% pause

%% Form Transformation to enforce periodicity
% ======================================================================= %

T_per = Periodic_Boundary_Conditions(dof_sets);
n_dof_per = size(T_per.s0,2);
n_dof = size(T_per.s0,1);


%% Form Master Mass and Stiffness
% ======================================================================= %
% store mesh info in a structure
Mesh_info.Ds            = Ds;
Mesh_info.rhos          = rhos;
Mesh_info.pattern_vec   = pattern_vec;
Mesh_info.n             = n;
Mesh_info.n_dof         = n_dof;
Mesh_info.elenodes      = element_node_table;
Mesh_info.coordinates   = coordinates;

% define model save path
pathstring = 'save_data/models/';

% define element order string
orderstrings = {'1st','2nd','3rd','4th','5th',...
                '6th','7th','8th','9th','10th'};   
orderstring = orderstrings{n-1};

% full model description string
modeldescription = [sprintf('SphericalResonator_%iDOF_',n_dof_per),...
                    orderstring,'OrderEles'];

model_savestring = [pathstring,modeldescription,'_FreeModel'];% form master mass and stiffness matrices


if exist([model_savestring,'.mat'],'file') && load_results
    load([model_savestring,'.mat'])
    save(model_savestring,'K_free','M_free','coordinates')
else
    [K_free,M_free] = master_mass_stiffness_3DE(Mesh_info);    
    if save_results
        save(model_savestring,'K_free','M_free','coordinates')
    end
end

%% compute full system dispersion
% ======================================================================= %

if full_soln
    
    % if full dispersion results exist for specified model size, then load
    % them. Otherwise compute them.
    solutionpathstring = 'save_data/solutions/';

    solutiondescription = [sprintf('%i',n_kap),sprintf('kpts_%iBands',n_curves)];

    solution_savestring = [solutionpathstring,modeldescription,'_',solutiondescription];   
      
                
    if exist([solution_savestring,'.mat'],'file') && load_results
        load([solution_savestring,'.mat'])
    else
        
        % display up front timing info
        fprintf('\nFull Dispersion Computation\n') 
        fprintf('======================================\n') 
        fprintf('No. DOF (after periodicity enforced): %i\n\n',n_dof_per) 
        t_full_start = tic;
        [omega_full,PHI_full,t_kloop_full] = dispersion_solver_w_k(kappa,K_free,M_free,dof_sets,R,n_curves);
        t_full = toc(t_full_start);
        
        if save_results
            size_PHI = n_dof*n_curves*n_kap*8*2;
            if size_PHI > 2e9
                save(solution_savestring,'omega_full','PHI_full','t_kloop_full','-v7.3')
            else
                save(solution_savestring,'omega_full','PHI_full','t_kloop_full')            
            end
        end
    end
    f_full = real(omega_full)/(2*pi);
    t_full = sum(t_kloop_full);


    figure(3);clf('reset')
    x_tick_labels = {'M/10','\Gamma','X/5'};
    h_disp1 = dispersion_plot(kappa_plot,{f_full},x_tick_labels);
    set(h_disp1{1},'linestyle','-','color','k','linewidth',2,'marker','.');
end



%% Perform BMS Plus Reduction
% =================================================================== %

do_BMS = true;
if do_BMS
    
    %w_cut_max = max(max(f_full*2*pi));
    clear options_BMS
    
    %w_cut_mult = 1.25;

    % interior reduction parameters
    options_BMS.InteriorMethod       = 'CB+';
    options_BMS.n_FI                 = 250;
    %options_BMS_pl.w_i                  = 1.0*w_cut_max*w_cut_mult;

    % boundary reduction parameters
    options_BMS.BoundaryMethod       = 'exact';
    options_BMS.n_CC                 = 40;
    
    
%     options_BMS.BoundaryMethod       = 'hybrid';
%     options_BMS.n_CC                 = 38;
    %options_BMS_pl.w_b                  = 24.0*w_cut_max*w_cut_mult;

    % additional options
    options_BMS_pl.verbose              = true;

    % perform BMS reduction 
    % [K_BMSpl,M_BMSpl,dof_sets_BMS,t_up_front_plus,T_BMS] = BMS_plus(K_free,M_free,coordinates,w_cut*0.5,R,n_FI,n_LI);
    [K_BMS,M_BMS,dof_sets_BMS,t_up_front,T_BMS] = BMS(K_free,M_free,coordinates,R,options_BMS);

    % solve for dispersion
    [w_BMS,PHI_BMS,t_kloop_BMS] = dispersion_solver_w_k(kappa,K_BMS,M_BMS,dof_sets_BMS,R,n_curves);
    f_BMS = w_BMS/(2*pi);

    % norm_fac = sqrt(rho1/E1)*Lx;
    figure(3);hold on
    dispersion_plot(kappa_plot,{f_BMS})
    drawnow

    tol = 1e-3;
    e_freq2 = 100*abs(f_BMS-f_full)./f_full;
    e_freq_max2 = max(max(e_freq2(f_full>tol*max(max(f_full)))))
    % pause

    t_BMS = sum(t_kloop_BMS) + t_up_front

    t_full = sum(t_kloop_full)
   
end


%% Compute BMS dispersion using k-w method
% ======================================================================= %
k_of_w = false;
if k_of_w
    n_om = 200;
    omega = linspace(0,w_cut_max,n_om);
    [kappa_BMS,~,t_wloop_BMS] = dispersion_solver_k_w(omega,K_BMS,[],M_BMS,dof_sets_BMS,R);

    figure(5);clf;hold on;view(3)
    plot3(real(kappa_BMS),imag(kappa_BMS),omega,'g.');hold on

    % plot BMS w(k) solution
    k_plot = ceil(n_kap/2):n_kap;
    xlim([-pi/Lx,pi/Lx]*0.20);ylim([-pi/Lx,pi/Lx]*2);


    kappa_BMS_real = abs(real(kappa_BMS));
    kappa_BMS_real(abs(imag(kappa_BMS))>0.001*pi/Lx) = nan;
    kappa_BMS_imag = abs(imag(kappa_BMS));
    kappa_BMS_imag(abs(real(kappa_BMS))>0.001*pi/Lx) = nan;

    figure(4);clf
    subplot(1,2,2)
    plot(kappa_BMS_real,omega/(2*pi),'g.','linewidth',2);hold on
    % plot(kappa(1,k_plot),w_BMS(:,k_plot),'ro','markersize',4)
    plot(kappa(1,k_plot),f_full(1:9,k_plot),'k-')
    xlim([0,pi/Lx]*0.20)
    ylim([0,max(omega)/(2*pi)])

    subplot(1,2,1)
    plot(-kappa_BMS_imag,omega/(2*pi),'g.','linewidth',2);hold on
    xlim([-pi/Lx,0]*0.2)
    ylim([0,max(omega)/(2*pi)])
    
    % save figure?
    if save_results
        figurepathstring = 'figures/';
        solutiondescription = [sprintf('kofw_%i',n_om),sprintf('Freqs_%iBands',n_curves)];
        figure_savestring = [figurepathstring,modeldescription,'_',solutiondescription];   

        saveas(gcf,figure_savestring)
    end
end

%% Sort Dispersion Branches
% ======================================================================= %

no_flat_lines = false;

PHI_full_sorted = PHI_full;
f_full_sorted = f_full;
PHI_BMS_sorted = PHI_BMS;
f_BMS_sorted = f_BMS;

if no_flat_lines
    tol = 5e-2;
    n_back = 5;
    
    if full_soln
        for i = 2:n_kap
            fprintf('full dispersion sorting: k-pt %i of %i\n',i,n_kap) 
            norm_mat = 1;
            
            % compute aggragate cross orthogonality (for current k-point
            % with previous 5)
            PHI2 = PHI_full_sorted(:,:,i);
            PHI2 = PHI2*diag(diag(PHI2'*norm_mat*PHI2).^(-0.5));
            
            Cmat = zeros(n_curves);
            
            mode_list = max(i-n_back,1):(i-1);
            for j = mode_list
                PHI1 = PHI_full_sorted(:,:,j);
                PHI1 = PHI1*diag(diag(PHI1'*norm_mat*PHI1).^(-0.5));
                Cmat = Cmat-abs(PHI2'*norm_mat*PHI1);
            end
            
            % use linear programming algorithm to find optimal matching
            [i_sort] = match_points(Cmat);
            
            % re-sort frequencies and modes
            f_full_sorted(:,i) = f_full_sorted(i_sort,i);
            PHI_full_sorted(:,:,i) = PHI_full_sorted(:,i_sort,i);
        end
        
        % put interesting modes first and flat modes second
        i_keep_full = (max(f_full_sorted,[],2)-min(f_full_sorted,[],2)) > tol*max(f_full_sorted,[],2);
        f_full_sorted = [f_full_sorted(i_keep_full,:);f_full_sorted(~i_keep_full,:)];
        PHI_full_sorted = cat(2,PHI_full_sorted(:,i_keep_full,:),PHI_full_sorted(:,~i_keep_full,:));
    end
    
    
    norm_mat = 1;
    for i = 2:n_kap                
        fprintf('BMS dispersion sorting: k-pt %i of %i\n',i,n_kap) 

        % compute aggragate cross orthogonality (for current k-point
        % with previous few)
        %PHI2 = PHI_BMS(:,:,i);
        kvec = kappa(:,i);
        PHI2 = BMS_Plus_Mode_Expansion(PHI_BMS_sorted(:,:,i),...
                        dof_sets_BMS,kvec,R,T_BMS,K_BMS,M_BMS);
        PHI2 = PHI2*diag(diag(PHI2'*norm_mat*PHI2).^(-0.5));
    
        Cmat = zeros(n_curves);
%         n_back = 2;
        mode_list = max(i-n_back,1):(i-1);
        for j = mode_list
            kvec = kappa(:,j);
            %PHI1 = PHI_BMS(:,:,j);            
            PHI1 = BMS_Plus_Mode_Expansion(PHI_BMS_sorted(:,:,j),...
                        dof_sets_BMS,kvec,R,T_BMS,K_BMS,M_BMS);
            PHI1 = PHI1*diag(diag(PHI1'*norm_mat*PHI1).^(-0.5));
            Cmat = Cmat-abs(PHI2'*norm_mat*PHI1);
        end

        % use linear programming algorithm to find optimal matching
        [i_sort] = match_points(Cmat);

        % re-sort frequencies and modes
        f_BMS_sorted(:,i) = f_BMS_sorted(i_sort,i);
        PHI_BMS_sorted(:,:,i) = PHI_BMS_sorted(:,i_sort,i);
    end

    % put interesting modes first and flat modes second
    i_keep_BMS = (max(f_BMS_sorted,[],2)-min(f_BMS_sorted,[],2)) > tol*max(f_BMS_sorted,[],2);
    f_BMS_sorted = [f_BMS_sorted(i_keep_BMS,:);f_BMS_sorted(~i_keep_BMS,:)];
    PHI_BMS_sorted = cat(2,PHI_BMS_sorted(:,i_keep_BMS,:),PHI_BMS_sorted(:,~i_keep_BMS,:));
end

%% Plot Dispersion comparison for BMS
% ======================================================================= %

% modes to keep in foreground (not faded)
i_keep_full = (max(f_full_sorted,[],2)-min(f_full_sorted,[],2)) > tol*max(f_full_sorted,[],2);
i_keep_BMS = (max(f_BMS_sorted,[],2)-min(f_BMS_sorted,[],2)) > tol*max(f_BMS_sorted,[],2);

figure(4);clf('reset')
x_tick_labels = {'M/10','\Gamma','X/5'};

% plot flat branches
omega_group2 = {f_full_sorted(~i_keep_full,:)/1e3,f_BMS_sorted(~i_keep_BMS,:)/1e3};    
h_disp2 = dispersion_plot(kappa_plot,omega_group2,x_tick_labels);
    set(h_disp2{1},'linestyle','-','color',[1,1,1]*0.8,'linewidth',2,'marker','none');
    set(h_disp2{2},'linestyle',':','color',[1,0.8,0.8],'linewidth',2,'marker','.');
    set(h_disp2{2},'markersize',4,'markerfacecolor',1-(1-[1,0.8,0.8])*1);
hold on 

% plot interesting branches
omega_group1 = {f_full_sorted(i_keep_full,:)/1e3,f_BMS_sorted(i_keep_BMS,:)/1e3};
h_disp1 = dispersion_plot(kappa_plot,omega_group1,x_tick_labels);
    set(h_disp1{1},'linestyle','-','color',[1.0,1.0,1.0]*0,'linewidth',2,'marker','none');
    set(h_disp1{2},'linestyle',':','color',[1.0,0.2,0.2]*1.0,'linewidth',2,'marker','.');
    set(h_disp1{2},'markersize',4,'markerfacecolor',[1,0.2,0.2]*1);

% plot legend
legendstrings = {'full','BMS'};
legend([h_disp1{1}(1),h_disp1{2}(1)],legendstrings)

% save figure?
if save_results
    figurepathstring = 'figures/';
    solutiondescription = [sprintf('%i',n_kap),sprintf('kPts_%iBands',n_curves)];
    figuredescription = 'DispersionComparison';
    figure_savestring = [figurepathstring,modeldescription,'_',solutiondescription,'_',figuredescription];

    saveas(gcf,figure_savestring)
end

%% Compute Error in Mode shapes and Frequencies
% ======================================================================= %

% tolerance for minimum amount of variation in dispersion curve
n_curves_keep = 6;

% preallocate arrays
e_f_BMS = zeros(n_curves,n_kap);
e_PHI_BMS = zeros(n_curves,n_kap);

for j2 = 1:n_kap
    tstart = tic;
    % wave vector at current k-point
    kvec = kappa(:,j2);
    
    % phase modification at boundaries due to wavevector
    lam  = exp(-1i*kvec'*R);
    lam(end+1:3) = 0;
    
    T_per_k = T_per.s0 + T_per.s1*lam(1) + T_per.s2*lam(2) + T_per.s3*lam(3)...
            + T_per.s12*lam(1)*lam(2) + T_per.s23*lam(2)*lam(3) + T_per.s13*lam(1)*lam(3)...
            + T_per.s123*lam(1)*lam(2)*lam(3);
        
    f_1 = f_full(:,j2);
    PHI_1 = T_per_k*PHI_full(:,:,j2);

    f_2 = f_BMS(:,j2);
    PHI_2 = BMS_Plus_Mode_Expansion(PHI_BMS(:,:,j2),...
                        dof_sets_BMS,kvec,R,T_BMS,K_BMS,M_BMS);

    % compute error in frequencies and mode shapes
    [e_PHI,e_f,~,~] = eig_error(PHI_1,f_1,PHI_2,f_2);
    e_PHI_BMS(:,j2) = e_PHI;
    e_f_BMS(:,j2) = e_f;
    
    
    [j2,toc(tstart)]
    
end

i_zero = f_full < max(f_full(:))*1e-3;
disp('mode error')
max(e_PHI_BMS(:))
ephimaxs = max(e_PHI_BMS,[],2);
[ephimaxs,isort] = sort(ephimaxs,'descend');
[ephimaxs,isort]


disp('frequency error')
max(e_f_BMS(:))
max(e_f_BMS(~i_zero))

% pause
figure(88);clf
plot3(kappa_plot,f_full(1:149,:),e_PHI_BMS(1:149,:))
xlabel('wave vector');ylabel('frequency (Hz)');zlabel('mode_error');
% saveas('test_error_figure');
% pause
%% Compute Strain displacement matrix (B) for each element
% ======================================================================= %
Bs = zeros(6,3*n^3,n^3,n_elements);

for i = 1:n_elements
    i
    coords_ele = coordinates(element_node_table(i,:),:);
    Bs(:,:,:,i) = element_BMatrix_3DE(n,coords_ele);
end

% determine number of elements each node is a part of
element_table_vec = (element_node_table(:));
node_uses = histc(element_table_vec,[1:n_nodes]);
element_dof_table = zeros(n_elements,3*n^3);
for i = 1:n_elements
    for j = 1:n^3
        element_dof_table(i,3*j-[2,1,0]) = 3*element_node_table(i,j)-[2,1,0];
    end
end

%% plot mode shapes
% ======================================================================= %
%load('data/Spherical_Resonator_modes28566DOF_10-2_15-21.mat') 
%load('data/Spherical_Resonator_workspace28566DOF_10-2_15-21.mat')

[~,~,~,~,~,patchfaces_lite_b,~,~] = ...
     mesh_spherical_resonator(n,n_ele_sc,n_ele_sph,n_ele_layer,...
         n_ele_bc,r_scube,r_bcube,r_spher,r_layer,false);
     
show_stress = true;
make_vid = false; 
% k_sel = (n_kap+1)/2-1; 
k_sel = 65;
mode_plot = 1:12;
% mode_plot = 5;

% p = 22;
xrange = max(xlocs)-min(xlocs);
yrange = max(ylocs)-min(ylocs);
zrange = max(zlocs)-min(zlocs);

xlimits = [min(xlocs)-0.15*xrange,max(xlocs)+0.15*xrange];
ylimits = [min(xlocs)-0.15*xrange,max(xlocs)+0.15*xrange];
zlimits = [min(xlocs)-0.15*xrange,max(xlocs)+0.15*xrange];

if length(mode_plot) >= 1
    
    mode_amp = Lx*0.1;
    PHI_plots = BMS_Plus_Mode_Expansion(PHI_BMS_sorted(:,mode_plot,k_sel),...
                    dof_sets_BMS,kappa(:,k_sel),R,T_BMS,K_BMS,M_BMS);

    for q = 1:length(mode_plot)


        if make_vid
            % video save string
            dtime = date;
%             ctime = clock;
            vidstring = sprintf('videos/mode%i_kpt_%i_of_%i_%i-%i_%i-%i',...
                mode_plot(q),k_sel,n_kap,ctime(2),ctime(3),ctime(4),ctime(5));
%             vidobj = VideoWriter(vidstring,'Motion JPEG AVI');
%             vidobj = VideoWriter(vidstring,'Uncompressed AVI');
            vidobj = VideoWriter(vidstring,'MPEG-4');
            vidobj.FrameRate = 15;
            vidobj.Quality = 100;
            open(vidobj);
        end

        PHI_plot = PHI_plots(:,q);
        f_plot = f_BMS_sorted(mode_plot(q),k_sel);
%         end
        [~,i_max] = max(PHI_plot);
        PHI_plot = mode_amp*PHI_plot/PHI_plot(i_max);

        figure(600);clf
        view(3)
%         title(['mode ',num2str(mode_plot(q)),', w = ',num2str(f_plot,3),' Hz'])

         xyzlocs = [xlocs,ylocs,zlocs];    
        [h_patch,h_line] = plot_FEM_model(xyzlocs,patchfaces_lite,[],fedges);

        axis equal
        axis manual
        set(gcf,'color','w')

        delete(h_patch)
        delete(h_line)
        
        % add title
        title(['Mode ',num2str(mode_plot(q)),', ',...
              num2str(f_BMS_sorted(mode_plot(q),k_sel)/1000),' kHz'])
          
        if make_vid
            nt = 60;
            tvec = linspace(0,(nt/15)*pi,nt);
            tvec(end) = [];
            nt = nt-1;
        else
            nt = 1;
            tvec = linspace(0,pi,nt);
        end
        
        
        % set colormap (get rid of parula)
        colormap(viridis);
        cmap = colormap;
        for p = 1:nt
            
            % define new node positions based on mode displacements
            xnewlocs = xlocs + real(PHI_plot(1:3:end)*exp(1i*tvec(p)));
            ynewlocs = ylocs + real(PHI_plot(2:3:end)*exp(1i*tvec(p)));
            znewlocs = zlocs + real(PHI_plot(3:3:end)*exp(1i*tvec(p)));

            xyznewlocs = [xnewlocs,ynewlocs,znewlocs];      
            
            [h_patch,h_line] = plot_FEM_model(xyznewlocs,patchfaces_lite_b,[],fedges);
            set(h_patch,'edgecolor','none')
            [h_patch_2,h_line] = plot_FEM_model(xyznewlocs,patchfaces_lite,[],fedges);
            set(h_patch_2,'facecolor','none')
%             break
            
            color_option = 2;
            
            if color_option == 0
                set(h_patch,'FaceColor','flat')
                set(h_patch,'FaceVertexCData',C_lite)
                set(h_patch,'facealpha',1)
                
            else
                
                if color_option == 1
                    % use displacement magnitude to determine colors
                    disp_vec = real(PHI_plot(1:3:end)*exp(1i*tvec(p))).^2 + ...
                               real(PHI_plot(2:3:end)*exp(1i*tvec(p))).^2 + ...
                               real(PHI_plot(3:3:end)*exp(1i*tvec(p))).^2;
                    disp_vec = sqrt(disp_vec);
                    color_vec = disp_vec;
                    
                elseif color_option == 2 
                % use von-mises stress to determine colors
                
                    % compute new stresses based on current displacements
                    sigma = zeros(n_nodes,6);
                    for r = 1:n_elements
                        for k = 1:n^3
                            B_element = Bs(:,:,k,r);
                            displacement_element = real(PHI_plot(element_dof_table(r,:))*exp(1i*tvec(p)));
                            sigmas_element = Ds{pattern_vec(r)}*B_element*displacement_element;
                            sigma(element_node_table(r,k),:) = ...
                            sigma(element_node_table(r,k),:) + ...
                            sigmas_element'/node_uses(element_node_table(r,k));
                        end
                    end

                    % compute von-mises stress from current stress vector
                    sigma_von_mises = sqrt((1/2)*((sigma(:,1)-sigma(:,2)).^2 + ...
                             (sigma(:,2)-sigma(:,3)).^2 + ...
                             (sigma(:,3)-sigma(:,1)).^2 + ...
                             6*(sigma(:,4).^2 + sigma(:,5).^2 + sigma(:,6).^2)));
                    color_vec = sigma_von_mises;                
                end
                
                
                color_vec = log(color_vec);
                
                % find min and max values of color range for plotting
                if p == 1
                    color_max = max((color_vec(unique(patchfaces_lite(:)))));
%                     color_min = min((color_vec(unique(patchfaces_lite(:)))));
                    color_avg = sum((color_vec(unique(patchfaces_lite(:)))))/...
                        length(color_vec(unique(patchfaces_lite(:))));
                    color_min = color_avg-(color_max-color_avg);
                end 
                
                % map color range to colormap scale
                l_cmap = length(colormap);
                color_vals = l_cmap*((color_vec)-...
                    color_min)/(color_max-color_min);

                color_vals(color_vals>size(cmap,1)) = size(cmap,1);
                color_vals(color_vals<1) = 1;
    %             color_vals_lite = color_vals(i_patch_nodes_lite,:);
                
                set(h_patch,'FaceColor','interp')
                set(h_patch,'FaceVertexCData',color_vals)
                set(h_patch,'CDataMapping','direct')
                
                % add colorbar to plot
                cb = colorbar;
                ylimitscb = get(cb, 'YLim');
                n_ticks = 5;
                set(cb, 'YTick', linspace(ylimitscb(1), ylimitscb(2), n_ticks));
                cb_tickvals = exp(color_min+linspace(0, l_cmap, n_ticks)*(color_max-color_min)/l_cmap);
                cb_ticklabels = cell(1,n_ticks);
                for i = 1:n_ticks
                    if i == 1
                        cb_ticklabels{i} = ...
                            ['<',num2str(cb_tickvals(i),'%4.2e')];
                    elseif i == n_ticks
                        cb_ticklabels{i} = ...
                            ['>',num2str(cb_tickvals(i),'%4.2e')];
                    else
                        cb_ticklabels{i} = ...
                            num2str(cb_tickvals(i),'%4.2e');
                    end
                end
                set(cb, 'YTickLabel', cb_ticklabels)
            end
            

            % remove axis and make background color to white
%             axis off
%             set(gcf,'color','w')

            % set renderer
            set(gcf,'Renderer','opengl');

            % write current frame to video file
            if make_vid
                set(gcf,'Renderer','zbuffer');
                frame = getframe(gcf);
                writeVideo(vidobj,frame);
            end

            % display figure then delete current patches
            drawnow
            if p ~=nt                    
                delete(h2)             
                delete(hline)
            end
        end
        
        if make_vid
            close(vidobj);
        end
        
        
        % set axis scaling so whole mode is visible
        xlim(xlimits)
        ylim(ylimits)
        zlim(zlimits)
        
        % save figure?
        if save_results
            figurepathstring = 'figures/';
            solutiondescription = [sprintf('%i',n_kap),sprintf('kPts_%iBands',n_curves)];
            figuredescription = sprintf('Mode%i_kpt%iof%i',mode_plot(q),k_sel,n_kap);
            figure_savestring = [figurepathstring,modeldescription,'_',figuredescription];

            saveas(gcf,figure_savestring)
        end
    end
end


%% Evaluate BMS AMLS+ LIR Reduction
% ======================================================================= %
eval_BMS = false;
if eval_BMS
    
    f_full_fsort = sort(f_full);
    
    % run parameters
    run_select = 'Full';
    switch run_select
        case 'Full' % full run values 
            
            % Reduction types
            InteriorMethods = {'CB'     ,'CB'       ,'CB+'      ,'CB+'};
            BoundaryMethods = {'exact'  ,'none'     ,'exact'    ,'none'};
            
            % multiplier to n_FI for number of CC interface modes to calculate
            %b_mults = [0.2,-1];        
            b_mults =  0.15*ones(size(InteriorMethods));
            
            % number of fixed-interface modes to try
            %n_FIs = ceil(logspace(log10(5),log10(500),10));
            %n_FIs = [125,150,200,250,300,400,500];
            n_FIs = [125,150,200,250,500,1000,1500];
                        
        case 'Small' % small test values
    
            % Reduction types
            InteriorMethods = {'CB'     ,'CB'       ,'CB'       ,'CB+'      ,'CB+'      ,'CB+'};
            BoundaryMethods = {'exact'  ,'hybrid'   ,'none'     ,'exact'    ,'hybrid'   ,'none'};
            
                        
            % multiplier to n_FI for number of CC interface modes to calculate
            %b_mults = [0.2,-1];        
            b_mults =  0.15*ones(size(InteriorMethods));
            
            % number of fixed-interface modes to try
            n_FIs = ceil(logspace(log10(150),log10(300),5));
            n_FIs = [125,150,200];
            n_FIs = [125,150,200,300];
            
            % multiplier to n_FI for number of CC interface modes to calculate
            %b_mults = [0.2,-1];        
%             b_mults = [0.15,-1,0.15,0.15,-1];   
    end

    % sizes needed for preallocation
    n_b_mults = length(b_mults);
    n_n_FIs = length(n_FIs);
    n_i_meths = length(InteriorMethods);
    
    % plotting values
    linestyles = {'-','--',':','-.'};
    colors = get(groot,'factoryAxesColorOrder');

    % frequencies to check in Dispersion calculation
    i_check = abs(f_full_fsort)>1e-1;

    % preallocate arrays
    n_FI_save               = nan(n_i_meths,n_n_FIs);
    n_LI_save               = nan(n_i_meths,n_n_FIs);
    e_f_maxs                = nan(n_i_meths,n_n_FIs);
    
    t_BMSs                  = zeros(n_i_meths,n_n_FIs);
    t_BMS_ks                = zeros(n_i_meths,n_n_FIs);
    t_kloop_BMSs            = zeros(n_i_meths,n_n_FIs,n_kap);
    t_up_fronts             = zeros(n_i_meths,n_n_FIs);    
    
    legendstrings_BMS       = cell(n_i_meths,1);
    
    for k = 1:n_i_meths
        for j = 1:n_n_FIs

            % multiplier for defining number of boundary modes
            b_mult = b_mults(k);            

            % Setup BMS options
            clear options_BMS
            options_BMS.InteriorMethod      = InteriorMethods{k};
            if b_mult<=0
                options_BMS.BoundaryMethod  = 'none';
            else
                options_BMS.BoundaryMethod  = BoundaryMethods{k};

                if strcmpi(BoundaryMethods{k},'hybrid')
                    options_BMS.n_CC = ceil(n_FIs(j)*b_mult*2);
                else
                    options_BMS.n_CC = ceil(n_FIs(j)*b_mult);
                end
            end
            options_BMS.n_FI                = n_FIs(j);
            options_BMS.verbose             = true;
            options_BMS.plots               = true;
            options_BMS.orthoTypeLIRExact   = 'qr';
            options_BMS.orthoTypeLIRExact   = 'svd';


            % Define a Save path
            % All this just define a string and path for saving?
            % sadly yes
            BMSsolutionpathstring = 'save_data/BMS_solutions/';

            solutiondescription = sprintf('%ikpts_%iBands',n_kap,n_curves);

            interior_type = [options_BMS.InteriorMethod,'_'];
            interior_size = sprintf('%iFImodes_',options_BMS.n_FI);
%             if b_mult<=0
%                 boundary_type = 'NoBoundReduction';
%                 boundary_size = [];
%                 boundary_orthog = [];
%             else
            boundary_type = [upper(options_BMS.BoundaryMethod(1)),lower(upper(options_BMS.BoundaryMethod(2:end))),'BoundReduct_'];
            boundary_size = sprintf('%iInitialCCModes',options_BMS.n_CC);
%             if isfield(options_BMS,'orthoTypeLIRExact');
%                 boundary_orthog = ['BoundOrthog_',upper(options_BMS.orthoTypeLIRExact)];
%             else
%                 boundary_orthog = 'BoundOrthogQR';
%             end
%             end

%             BMSdescription = [interior_type,interior_size,boundary_type,boundary_size,...
%                               boundary_orthog]; 
                          
            BMSdescription = [interior_type,interior_size,boundary_type,boundary_size];

            BMSsolution_savestring = [BMSsolutionpathstring,modeldescription,...
                            '_',solutiondescription,'_',BMSdescription];


            if exist([BMSsolution_savestring,'.mat'],'file')  && load_results
                load([BMSsolution_savestring,'.mat'])
            else

                %[K_BMSpl,M_BMSpl,dof_sets_BMSpl,t_up_front_plus_save(i,j),T_BMS_plus] = BMS_plus(K_free,M_free,coordinates,R,opts_BMS);
                [K_BMS,M_BMS,dof_sets_BMS,t_up_front_BMS,T_BMS_temp] = BMS(K_free,M_free,coordinates,R,options_BMS);
                %[K_BMS,M_BMS,dof_sets_BMS,t_up_fronts(k,i,j)] = BMS(K_free,M_free,coordinates,R,options_BMS);

                % compute BMS dispersion
                [w_BMS,PHI_BMS,t_kloop_BMS] = ...
                 dispersion_solver_w_k(kappa,K_BMS,M_BMS,dof_sets_BMS,R,n_curves);

                % convert to Hz
                f_BMS = w_BMS/(2*pi);

                % save BMS results
                if save_results
                    save(BMSsolution_savestring,'f_BMS','PHI_BMS','t_kloop_BMS','t_up_front_BMS')
                end
            end

            % store results
            f_BMS_save{k,j} = f_BMS;
            PHI_BMS_save{k,j} = PHI_BMS;
            t_kloop_BMSs(k,j,:) = t_kloop_BMS;
            t_up_fronts(k,j) = t_up_front_BMS;

            % timing data
            t_BMSs(k,j) = sum(t_kloop_BMSs(k,j,:)) + t_up_fronts(k,j);
            t_BMS_ks(k,j) = sum(t_kloop_BMSs(k,j,:))/n_kap;

            % model dimensions
            n_FI_save(k,j) = n_FIs(j);
            n_LI_save(k,j) = size(PHI_BMS,1)-n_FI_save(k,j);

            % evaluate error and store maximum error
            e_f_BMS = zeros(size(f_full_fsort));
            e_f_BMS(i_check) = 100*(f_BMS(i_check)-f_full_fsort(i_check))./f_full_fsort(i_check);


            e_f_maxs(k,j) = max(max(abs(e_f_BMS)));

            % create legend strings
            legendstrings_BMS{k} = ['Interior: ',InteriorMethods{k},', Bound: ',BoundaryMethods{k},' (b mult. = ',num2str(b_mult),')'];

            % plot intermediate results
            figure(12);clf
            for qq = 1:k
                h1 = semilogy(squeeze(t_BMSs(qq,:))/t_full,squeeze(e_f_maxs(qq,:)));hold on

                set(h1,'marker','.',...
                       'linewidth',2,...
                       'markersize',12,...
                       'color',colors(qq,:));
                xlabel('Time Fraction, R');ylabel('Maximum Frequency Error (%)')
            end
            legendstrings_temp = legendstrings_BMS(1:k);
            legend(legendstrings_temp{:},'location','northeast')
            drawnow

            % save figure
            if save_results
                figure_pathstring = 'figures/';
                figure_description = ['PerformancePlot_',run_select,'Run'];
                performance_figure_savestring = [figure_pathstring,modeldescription,...
                                '_',solutiondescription,'_',figure_description];

                saveas(gcf,performance_figure_savestring)
            end

            % output some model-size info
            n_FI_save
            n_LI_save
        end
    end
end
