clear
clc

d_ac = zeros(10,251);
d_rt = zeros(10,251);
d_tv = zeros(10,251);
msd_rt=zeros(10,50);
msd_tv = zeros(10,50);
for qq = 1:1
load('Autocorrelation_velmag3.mat')
kkk =qq
Cvv = Cvv_mean10;
dl = 1;
lc = find(Cvv<0,1);
lc1 = (lc)-Cvv(lc)/(Cvv(lc)-Cvv(lc-1));
lc =lc1;

load('Pr_p10.mat')
remain = eval(strcat('Pr',num2str(kkk)));

%thing = 'C:/Users/tjshe/OneDrive/Documents/MATLAB/InverseLANL/p3_part_1/p3_part_1/x01/control_planes';
if qq<10
    %thing = strcat('C:/Users/tjshe/OneDrive/Documents/MATLAB/InverseLANL/p3_part_1/p3_part_1/x0', num2str(kkk), '/control_planes');
    thing = strcat('C:/Users/tjshe/OneDrive/Documents/MATLAB/InverseLANL/p10/p10_x0', num2str(kkk), '/control_planes');
    %thing = strcat('C:/Users/tjshe/OneDrive/Documents/MATLAB/InverseLANL/p5_part_4/part_4/x0', num2str(kkk), '/control_planes');
else
    %thing = strcat('C:/Users/tjshe/OneDrive/Documents/MATLAB/InverseLANL/p3_part_1/p3_part_1/x', num2str(kkk), '/control_planes');
    thing = strcat('C:/Users/tjshe/OneDrive/Documents/MATLAB/InverseLANL/p10/p10_x', num2str(kkk), '/control_planes');
    %thing = strcat('C:/Users/tjshe/OneDrive/Documents/MATLAB/InverseLANL/p5_part_4/part_4/x', num2str(kkk), '/control_planes');
end
%% Different Models
%bctrw - classic bernouli
%bctrw_rt - bernouli with random tortuosity
%bctrw_tv - bernouli with tortuosity and velocity coupled
%ctrw_rt - tortuosity and velocity uncoupled
%ctrw_tv - tortuosity and velocity couple

%% Load in data
sub2 =strcat(thing,'/Control_00.dat' );
A=importdata(sub2);
[ll,~] = size(A);
cell_length =1;%length of a cell
cells = floor(50/cell_length);
b =zeros(ll,cells);
v =zeros(ll,cells);
vx =zeros(ll,cells);
dist =zeros(ll,cells);
tt=zeros(ll,cells);
TM_num =10;
interval = floor(50/TM_num);
v_initial = sqrt(A(:,5).^2+A(:,6).^2+A(:,7).^2);
vx_initial = A(:,5);
N = length(v_initial);
aa = find(vx_initial ==0);
vx_initial(aa) = vx_initial(randi(N,length(aa),1));
 %number of particles
%% Sort data
for i = 1:cells
    j = i*cell_length;
    if j<10
        sub = strcat('0',num2str(j),'.dat');
    else
        sub = strcat(num2str(j),'.dat');
    end
    sub2 = strcat(thing, '/Control_', sub);
    %sub2  = strcat('C:/Users/tjshe/OneDrive/Documents/MATLAB/InverseLANL/p3_part_1/p3_part_1/x01/control_planes/Control_',sub);
    %sub2 = strcat('C:/Users/tjshe/OneDrive/Documents/MATLAB/InverseLANL/p10/p10_x01/control_planes/Control_', sub);
    A=importdata(sub2);
    b(:,i)= A(:,10); 
    v(:,i) = sqrt(A(:,5).^2+A(:,6).^2+A(:,7).^2);
    vx(:,i) = A(:,5);
    dist(:,i) = A(:,8);
    tt(:,i) = A(:,1);
end

vx = reshape(vx, [],1);
vx_l = length(vx);
cell_dist = dist;
transition_time =tt;
for i = 2:cells
   cell_dist(:,i) = cell_dist(:,i)-dist(:,i-1);
   transition_time(:,i) = transition_time(:,i)-tt(:,i-1);
end
figure(101)
imagesc(cell_dist)
colorbar

effective_vel = cell_dist./transition_time;
v_initial = effective_vel(:,1);
tourosity = mean(mean(cell_dist))/dl;
cell_dist(cell_dist <1) =1;
%v = reshape(v,[],1);
vv = reshape(v,[],1);
ttbins = logspace(-1,6,1000);
v_true = reshape(effective_vel,[],1);
tour = reshape(cell_dist,[],1);
TL = length(tour);
%% Velocity Lagrangian
load('FullLagrangianPDFs.mat')
vbins = vbins_full3;
vlag = vmat_lag_full3(kkk,:);
rr = randsample(100, 10^7, true, vlag);
v_samp = vbins_full3(rr); 
[~, v_true_ind] = histc(v_true,vbins);
tour_samp = zeros(10^7,1);
for i = 1:100
    aa = v_true_ind == i;
    bb = rr == i;
    if sum(aa)>0
    tour_temp = tour(aa);
    tour_samp(bb) = tour_temp(randi(length(tour_temp), sum(bb),1));
    else
    tour_samp(bb) = tour(randi(length(tour), sum(bb),1));
    end
    
end

%% Initialize total distance travel for each particle
d_rt = zeros(N,50);
d_tv = zeros(N,50);

%% Bernouli Walk
t_bctrw = zeros(N,2);
t_bctrw_tv = zeros(N,2);
t_bctrw_rt = zeros(N,2);
vt = v_initial; %travel velocity
tour_cor = cell_dist(:,1);
%remain = exp(-dl/lc); %probability we remain
for i = 1:50/dl
   P = rand(N,1);
   aa = P>remain;
   %ii = randi(10^7, sum(aa),1);
   %vt(aa) = v_samp(ii);
   ii = randi(length(vv), sum(aa),1);
   vt(aa) = vv(ii);
   %classic bernouli
   t_bctrw(:,1) = t_bctrw(:,1)+tourosity*dl./vt;
   t_bctrw(:,2) = t_bctrw(:,2)+tourosity*dl;
   %bernouli random tortuosity
   rrr =tour(randi(TL,N,1)); 
   t_bctrw_rt(:,1) = t_bctrw_rt(:,1) + rrr./vt;
   t_bctrw_rt(:,2) = t_bctrw_rt(:,2) + rrr.*dl;
   d_rt(:,i) = t_bctrw_rt(:,2);
   %bernouli correlated velocity and velocity
   %tour_cor(aa) = tour_samp(ii);
   bb = logical(abs(aa-1));
   rrr =tour(randi(TL,sum(bb),1));
   %tour_cor(bb) =rrr;
   %t_bctrw_tv(:,1) = t_bctrw_tv(:,1) + tour_cor./vt;
   %t_bctrw_tv(:,2) = t_bctrw_tv(:,2) + tour_cor.*dl;
   [~, vt_ind]= histc(vt, vbins);
   for mm = 1:100
       aa = vt_ind ==mm;
       if sum(aa)>0
       bb = rr == mm;
       if sum(bb)>0
          tour_temp = tour_samp(bb);
          tour_cor(aa) =tour_temp(randi(sum(bb), sum(aa),1));
       else
           tour_temp = tourosity;
           tour_cor(aa) =tourosity;
       end
       %tour_temp = tour_samp(bb);
       %tour_cor(aa) =tour_temp(randi(sum(bb), sum(aa),1));
       end
   end
   t_bctrw_tv(:,1) = t_bctrw_tv(:,1) + tour_cor./vt;
   t_bctrw_tv(:,2) = t_bctrw_tv(:,2) + tour_cor.*dl;
   d_tv(:,i) = t_bctrw_tv(:,2);
   
   % store different lengths
   if i ==15/dl
      t_bctrw15 = t_bctrw(:,1); 
      t_bctrw_tv15 = t_bctrw_tv(:,1);
      t_bctrw_rt15 = t_bctrw_rt(:,1);

   end
   if i ==30/dl
      t_bctrw30 = t_bctrw(:,1); 
      t_bctrw_tv30 = t_bctrw_tv(:,1);
      t_bctrw_rt30 = t_bctrw_rt(:,1);
   end
end
%% Calculate Mean displacement
msd_rt2 = zeros(50,1);
msd_tv2 = zeros(50,1);
for i =1:50
    msd_rt2(i) = sum((d_rt(:,i)-mean(d_rt(:,i))).^2)/N;
    msd_tv2(i) = sum((d_tv(:,i)-mean(d_tv(:,i))).^2)/N;
end

msd_rt(qq,:) = msd_rt2';
msd_tv(qq,:) = msd_tv2';



%% Figures
t_actual = tt(:,50);
t_actual15 = tt(:,15);
t_actual30 = tt(:,30);
x = logspace(0,8,1000);
%classic bernouli
c_bctrw = histc(t_bctrw(:,1),x);
c_bctrw15 = histc(t_bctrw15(:,1),x);
c_bctrw30 = histc(t_bctrw30(:,1),x);
eval(strcat('c_bctrw_',num2str(kkk), '= c_bctrw'));
eval(strcat('c_bctrw15_',num2str(kkk), '= c_bctrw15'));
eval(strcat('c_bctrw30_',num2str(kkk), '= c_bctrw30'));
%bernouli random tortuosity
c_bctrw_rt = histc(t_bctrw_rt(:,1),x);
c_bctrw15_rt = histc(t_bctrw_rt15(:,1),x);
c_bctrw30_rt = histc(t_bctrw_rt30(:,1),x);
eval(strcat('c_bctrw_rt_',num2str(kkk), '= c_bctrw_rt'));
eval(strcat('c_bctrw15_rt_',num2str(kkk), '= c_bctrw15_rt'));
eval(strcat('c_bctrw30_rt_',num2str(kkk), '= c_bctrw30_rt'));
%bernouli correlated velocity tortuosity
c_bctrw_tv = histc(t_bctrw_tv(:,1),x);
c_bctrw15_tv = histc(t_bctrw_tv15(:,1),x);
c_bctrw30_tv = histc(t_bctrw_tv30(:,1),x);
eval(strcat('c_bctrw_tv_',num2str(kkk), '= c_bctrw_tv'));
eval(strcat('c_bctrw15_tv_',num2str(kkk), '= c_bctrw15_tv'));
eval(strcat('c_bctrw30_tv_',num2str(kkk), '= c_bctrw30_tv'));
%Actual BTCs
c_actual = histc(t_actual(:,1),x);
c_actual15 = histc(t_actual15(:,1),x);
c_actual30 = histc(t_actual30(:,1),x);
eval(strcat('c_actual_',num2str(kkk), '= c_actual'));
eval(strcat('c_actual15_',num2str(kkk), '= c_actual15'));
eval(strcat('c_actual30_',num2str(kkk), '= c_actual30'));
 

%% FIGURES
figure(1)
loglog(x,c_actual,'.',x,c_bctrw, x, c_bctrw_rt, x, c_bctrw_tv);
legend('DFN', 'bctrw', 'RT', 'TV')

figure(2)
cdf_bctrw = c_bctrw;
cdf_bctrw_rt = c_bctrw_rt;
cdf_bctrw_tv = c_bctrw_tv;
cdf_actual = c_actual;
for i = 2:length(x)
    cdf_bctrw(i) = cdf_bctrw(i)+cdf_bctrw(i-1);
    cdf_bctrw_tv(i) = cdf_bctrw_tv(i)+cdf_bctrw_tv(i-1);
    cdf_bctrw_rt(i) = cdf_bctrw_rt(i)+cdf_bctrw_rt(i-1);
    cdf_actual(i) = cdf_actual(i)+cdf_actual(i-1);
end
cdf_bctrw = cdf_bctrw/N;
cdf_bctrw_tv = cdf_bctrw_tv/N;
cdf_bctrw_rt = cdf_bctrw_rt/N;
cdf_actual = cdf_actual/N;
semilogx(x, cdf_bctrw, x, cdf_bctrw_rt, x, cdf_bctrw_tv, ...
     x,cdf_actual, '.', 'Linewidth', 1.5);
legend('BCTRW', 'RT', 'TV', 'DFN')

figure(3*qq)
loglog(x, 1-cdf_bctrw,'*', x, 1-cdf_bctrw_rt, '*', x, 1-cdf_bctrw_tv,'*', ...
    x,1-cdf_actual, '.');
legend('BCTRW', 'RT', 'TV', 'DFN')
drawnow
%%
figure(4)
dist_bins = [50:1:300];
d_actual = histc(dist(:,50),dist_bins);
d_bctrw = histc(t_bctrw(:,2),dist_bins);
d_bctrw_rt = histc(t_bctrw_rt(:,2),dist_bins);
d_bctrw_tv = histc(t_bctrw_tv(:,2),dist_bins);

dc_ac(kkk,:)= d_actual;
dc_rt(kkk,:)= d_bctrw_rt;
dc_tv(kkk,:)= d_bctrw_tv;

plot(dist_bins, d_bctrw_rt, dist_bins, d_bctrw_tv, dist_bins, d_actual, '.')
end
