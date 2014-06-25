%LOAD DATA
close all, clear all;
Data = load('02D.txt.new');

% Define Variables

CC       = Data(:,1); %CycleCount/n

Y_force  = Data(:,7);
Y_force = abs(Y_force);

A_dist_Y = Data(:,9);

% Rescale distance axis
A_dist_Y = A_dist_Y-min(A_dist_Y);

Y_force = abs(Y_force);

% Find number of extension/relaxion cycles

[Peaks,Index] = findpeaks( Y_force,'minpeakdistance',100);
 
distlocs = A_dist_Y(Index);

M = [Index,Peaks];

[values,order] = sort(Peaks,'descend');

sX = M(order,:);

N = 20; % Alter this number to ensure all peaks are marked
 
maxValues = sX(:,2);
MaxValueIndices = sX(:,1);
 
maxValues = maxValues(1:N);
MaxValueIndices = MaxValueIndices(1:N);
  
figure;

plot(Y_force,'-b'); hold on;
plot(MaxValueIndices,Y_force(MaxValueIndices),'k^','markerfacecolor',[1 0 0]);
set(gca,'Linewidth',2,'FontSize',16,'FontWeight','bold','XMinorTick','on');
xlabel('Time [ms]','FontSize',16,'FontWeight','bold');
ylabel('Force [pN]','FontSize',16,'FontWeight','bold');

 sorted = sort(MaxValueIndices(:));
 
 lowforce = 2; % remove events below this force level
 highforce = 12;
% RLOESS FILTER
coef = 50; % 75 standard
coef2 = 100; % dx peak finding slips

% Preallocate output arrays 
dx_rlx = [2,1:N-2];
dx_ext = [2,1:N-2];

P1_rlx = [2,1:N-2];
P1_ext = [2,1:N-2];
P2_rlx = [2,1:N-2];
P2_ext = [2,1:N-2];


%Y_force = smooth(Y_force, 70, 'rlowess');

hold on;
for i=1:N-1
%for i = 3
    % Parse into Experiments
    f = Y_force(sorted(i):sorted(i+1));
    x = A_dist_Y(sorted(i):sorted(i+1));
    % SPLIT DATA into extension relaxation events 
    [minf_x,minf_y] = min(f);
    f_rlx = f(1:minf_y);
    f_ext = f(minf_y:end);
    x_rlx = x(1:minf_y);
    x_ext = x(minf_y:end);
    % Shorten ARRAYS to remove low force section
    pl1 = size(find(f_rlx >= lowforce),1);
    pl2 = size(find(f_ext <= lowforce),1);
    % Shorten ARRAYS to remove high force section
    pl11 = size(find(f_rlx >= highforce),1);
    pl22 = size(find(f_ext <= highforce),1);
    f_rlxs = f_rlx(pl11:pl1);
    x_rlxs = x_rlx(pl11:pl1);
    f_exts = f_ext(pl2:pl22);
    x_exts = x_ext(pl2:pl22);
    % FILTER DATA
    f_rlxss = smooth(f_rlxs,coef,'rlowess');% 
    f_extss = smooth(f_exts,coef,'rlowess');
    f_rlxs2 = smooth(f_rlxs,coef2,'rlowess');% 
    f_exts2 = smooth(f_exts,coef2,'rlowess');
    f_rlxds = diff(f_rlxs2);% Change from f_rlxss!!!!!
    f_extds = diff(f_exts2);    
%     [f_rlxdminx,f_rlxdminy] = min((f_rlxds));%min(abs(f_rlxds));
%     [f_rlxdmaxx,f_rlxdmaxy] = max((f_rlxds));
%     [f_extdminx,f_extdminy] = min((f_extds));
%     [f_extdmaxx,f_extdmaxy] = max((f_extds));
    % Find multiple cutting points
    [cut_rlx,cut_Index_rlx] = findpeaks( f_rlxds,'minpeakdistance',100); % 100 distance
    [cut_ext,cut_Index_ext] = findpeaks( -1 * f_extds,'minpeakdistance',100); % THIS IS THE -1 PART
    cut_distlocs_rlx = A_dist_Y(cut_Index_rlx);
    cut_distlocs_ext = A_dist_Y(cut_Index_ext);
    M_rlx = [cut_Index_rlx,cut_rlx];
    M_ext = [cut_Index_ext,cut_ext];
    [values_rlx,order_rlx] = sort(cut_rlx,'descend');
    [values_ext,order_ext] = sort(cut_ext,'descend'); 
    sX_rlx = M_rlx(order_rlx,:);
    sX_ext = M_ext(order_ext,:);
    N_cut = 2; % Alter this number to ensure all peaks are marked
    maxValues_rlx = sX_rlx(:,2);
    MaxValueIndices_rlx = sX_rlx(:,1);
    maxValues_rlx = maxValues_rlx(1:N_cut);
    MaxValueIndices_rlx = MaxValueIndices_rlx(1:N_cut);
    maxValues_ext = sX_ext(:,2);
    MaxValueIndices_ext = sX_ext(:,1);
    maxValues_ext = maxValues_ext(1:N_cut);
    MaxValueIndices_ext = MaxValueIndices_ext(1:N_cut);
    sorted_rlx = sort(MaxValueIndices_rlx(:),'ascend');
    sorted_ext = sort(MaxValueIndices_ext(:),'descend');

          for j = 1: N_cut
          %for j = 2
          f_rlxscut = f_rlxss(sorted_rlx(j):end); % WHAT ABOUT f_rlxss!!!!!!
          [f_rlxscutx, f_rlxscuty] = max(f_rlxscut);
          f_rlxp1 = (sorted_rlx(j) + f_rlxscuty -1);
          f_extscut = f_extss(1:sorted_ext(j));
          [f_extscutx, f_extscuty] = max(f_extscut);
          f_extp1 = f_extscuty;    
          % Find Equal Position
          f_rlxscont = f_rlxss(1:sorted_rlx(j));
          [f_rlxseqx, f_rlxseqy] = find(f_rlxscont >= f_rlxscutx); % >=
          f_rlxp2 = size(f_rlxseqy,1);
          f_extscont = f_extss(sorted_ext(j):end);
          [f_extseqx, f_extseqy] = find(f_extscont <= f_extscutx); % <=
          f_extp2 = (size(f_extseqy,1) + sorted_ext(j) -1);
          % Store P1 and P2 values
          P1_rlx(j,i) = x_rlxs(f_rlxp1);
          P1_ext(j,i) = x_exts(f_extp1);
          P2_rlx(j,i) = x_rlxs(f_rlxp2);
          P2_ext(j,i) = x_exts(f_extp2);
          % Calculate DELTA X in nm          
          dx_rlx(j,i) = x_rlxs(f_rlxp2) - x_rlxs(f_rlxp1);
          dx_ext(j,i) = x_exts(f_extp2) - x_exts(f_extp1);
          
          figure;
          plot(x_rlxs,f_rlxss,'-r',x_exts,f_extss,'-k','linewidth',2); hold on;
          %plot(x_rlxs,f_rlxss/max(f_rlxss),'-r',x_exts,f_extss/max(f_extss),'-k','linewidth',2); hold on; 
          %plot(x_rlxs(1:end-1),f_rlxds,'--r',x_exts(1:end-1),f_extds,'--k','linewidth',2); hold on; 
          
          plot(x_rlxs(f_rlxp1),(f_rlxss(f_rlxp1)),'r*',x_exts(f_extp1),(f_extss(f_extp1)),'k+','markerfacecolor',[0 0 1],'markersize',10); hold on;  
          
          plot(x_rlxs(f_rlxp2),(f_rlxss(f_rlxp2)),'r*',x_exts(f_extp2),(f_extss(f_extp2)),'k+','markerfacecolor',[0 0 1],'markersize',10); hold on;
          plot(x_rlxs(sorted_rlx),f_rlxss(sorted_rlx),'r^','markerfacecolor',[1 0 0]); hold on;
          %plot(x_rlxs(sorted_rlx),f_rlxds(sorted_rlx),'k^','markerfacecolor',[1 0 0]); hold on;
          plot(x_exts(sorted_ext),f_extss(sorted_ext),'k^','markerfacecolor',[0 0 1]);
          set(gca,'Linewidth',2,'FontSize',16,'FontWeight','bold','XMinorTick','on');% 300 420
          
          xlabel('Extension [ nm ]','FontSize',16,'FontWeight','bold');
          ylabel('Force [ pN ]','FontSize',16,'FontWeight','bold');
          
          figure;
          plot(x_rlxs(2:end),f_rlxds,'-r',x_exts(2:end),f_extds,'-k','linewidth',2); hold on;
          %plot(x_rlxs,f_rlxss/max(f_rlxss),'-r',x_exts,f_extss/max(f_extss),'-k','linewidth',2); hold on; 
          %plot(x_rlxs(1:end-1),f_rlxds,'--r',x_exts(1:end-1),f_extds,'--k','linewidth',2); hold on; 
          
          plot(x_rlxs(f_rlxp1),(f_rlxds(f_rlxp1)),'g*',x_exts(f_extp1),(f_extds(f_extp1)),'g+','markerfacecolor',[0 0 1],'markersize',10); hold on;  
          
          plot(x_rlxs(f_rlxp2),(f_rlxds(f_rlxp2)),'c*',x_exts(f_extp2),(f_extds(f_extp2)),'c+','markerfacecolor',[0 0 1],'markersize',10); hold on;
          plot(x_rlxs(sorted_rlx),f_rlxds(sorted_rlx),'k^','markerfacecolor',[1 0 0]); hold on;
          %plot(x_rlxs(sorted_rlx),f_rlxds(sorted_rlx),'k^','markerfacecolor',[1 0 0]); hold on;
          plot(x_exts(sorted_ext),f_extds(sorted_ext),'b^','markerfacecolor',[0 0 1]);
          set(gca,'Linewidth',2,'FontSize',16,'FontWeight','bold','XMinorTick','on');% 300 420
          
          xlabel('Extension [ nm ]','FontSize',16,'FontWeight','bold');
          ylabel('Force [ pN ]','FontSize',16,'FontWeight','bold');
          
          


%           % Calculate DELTA X in nm
%           dx_rlx(j,i) = x_rlxs(f_rlxp2) - x_rlxs(f_rlxp1);
%           dx_ext(j,i) = x_exts(f_extp2) - x_exts(f_extp1);
          end
          
  end





 