%% Code  for pubication RLTS Project :
clear all;
clc;
close all;

%% Figure 2 in the paper:
 
%the data is organized for each experimantal day at a cell, when each cell is organized by experimantal groups: undamaged tissue, sutured and soldered. 
data(:,1)={[196,208,208];[148,140,120];[64,106,118,26,40]}; % zigzag soldering 
data(:,2)={[196,156,194];[80,70,110];[50,40,40,40]}; % fast discrete soldering 
data(:,3)={[284,270,236];[110,150,152];[120,144,70]};% slow discrete soldering
data(:,4)={[254,360,272,246,272,272,220];[82,110,100];[132,132,76,136,118]}; % double slow discrete soldering

s_d=size(data); % get the data size for the loop 

figure('position',[100 100 350 500]);
AxPos=[0.2 0.3 0.75 0.6];
axes('position',AxPos);
hold on;

for i=1:s_d(2) % runs for the different experimantal days
    median_IT(i)=median(cell2mat(data(1,i))); % the median value for the undamaged tissue group
    
    % The data analysis for the suturing group:
    suture=100*cell2mat(data(2,i))./median_IT(i);
    meanS(i)=mean(suture);
    stdS=std(suture);
    b=bar(i-0.125,meanS(i),0.25,'FaceColor','flat');  % Plot the mean value of the suturing group
    hold on
    b.CData =[0.3010 0.7450 0.9330]; % suturing group color
    line([i-0.125 i-0.125],[meanS(i)-stdS meanS(i)+stdS],'color',[0 0 0],'LineWidth',3) % Plot the standard deviation value of the suturing group
   
    for j=1:length(suture) % Plot each indevidual value of the suturing group
    plot(i+rand(1)/10-0.125,suture(j),'Marker','o','MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none')
    end

    %The data analysis for the soldering group
    for n=3:s_d(1)
        if ~isempty(cell2mat(data(n,i)))
            solder(n-2,:)=100*cell2mat(data(n,i))./median_IT(i);
            meanRLTS(n-2,i)=mean(solder(n-2,:));
            stdSo=std(solder(n-2,:));
            b2=bar(i+0.125,meanRLTS(n-2,i),0.25,'FaceColor','flat');% Plot the mean value of the soldering group
            b2.CData =[40,106,70]./255; % soldering (albumin 800 and ICG 0.6) color 
            line([i+0.125 i+0.125],[meanRLTS(i)-stdSo meanRLTS(i)+stdSo],'color',[0 0 0],'LineWidth',3) %Plot the standard deviation value of the soldering group
            for j=1:length(solder) % Plot each indevidual value of the suturing group
            plot(i+rand(1)/10+0.125,solder(n-2,j),'Marker','>','MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none')
            end
            clear solder
        end
    end

end
ylabel('Normalized Burst Pressure  [%]')
set(gca,'xtick',1:4,'xticklabel',{'Zig-Zag','Fast discrete','Slow discrete','Double slow discrete'},'xlim',[0.5 4.5],'ylim',[0 100],'FontSize',16,'FontName','Times new Roman','Box','off')
xlabel('Experimental protocol')
%legend('Manual sutures','Robotic solder')
%% Statistical analysis for figure 2 dataset:

allM_data=[];
for i=1:s_d(2)
    median_IT(i)=median(cell2mat(data(1,i)));
    suture=100*cell2mat(data(2,i))./median_IT(i);
    solder=100*cell2mat(data(3,i))./median_IT(i);
    M_data(:,1)=[suture';solder'];
    M_data(:,2)=[-1*ones(length(suture),1);ones(length(solder),1)];
    M_data(:,3)=i.*ones(length(suture)+length(solder),1);
    allM_data=[allM_data; M_data];
    clear M_data
end
for con=2:s_d(2)
    C(allM_data(:,3)==con,con-1)=1;
    C(allM_data(:,3)==1,con-1)=-1;
end
tab=table(allM_data(:,1),allM_data(:,2),C(:,1),C(:,2),C(:,3),'VariableNames',{'y','method','exp2','exp3','exp4'}); % sumerized table with the data organized by  method factor ('sutuered' or 'soldered') and experiment factor(1-3) code using dummy codeing  

% matlab code analysis
[p,tbl,stats,terms]=anovan(allM_data(:,1),{allM_data(:,2),allM_data(:,3)},'model',2,'varnames',{'method','exp'});
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

% dummy codeing varuble analysis using linear model
mdl=fitlm(tab,'y~exp2+exp3+exp4+method+exp2:method+exp3:method+exp4:method');
mdlNOexp=fitlm(tab,'y~method+exp2:method+exp3:method+exp4:method');
mdlNOm=fitlm(tab,'y~exp2+exp3+exp4+exp2:method+exp3:method+exp4:method');
mdlNoI=fitlm(tab,'y~exp2+exp3+exp4+method');

yhat=mdl.Fitted;
yhatNOexp=mdlNOexp.Fitted;% model fit without the exp factor
yhatNOm=mdlNOm.Fitted; % model fit without the method factor
yhatNOI=mdlNoI.Fitted;% model fit without the interaction between exp factor and method factor 

b=mdl.Coefficients.Estimate;
nu_res=mdl.DFE;
MS_res=mdl.MSE;

[p_all,F_all] = coefTest(mdl);
dfm=1;dfExp=3; dfI=3;

SS_NOexp=sum((yhat-yhatNOexp).^2);
SS_NOm=sum((yhat-yhatNOm).^2);
SS_NOI=sum((yhat-yhatNOI).^2);

MS_NOexp=SS_NOexp/dfExp;
F_NOexp=MS_NOexp/MS_res;
p_EXp=1-fcdf(F_NOexp,dfExp,nu_res);

MS_NOm=SS_NOm/dfm;
F_NOm=MS_NOm/MS_res;
p_NOm=1-fcdf(F_NOm,dfm,nu_res);

MS_NOI=SS_NOI/dfI;
F_NOI=MS_NOI/MS_res;
p_NOI=1-fcdf(F_NOI,dfI,nu_res);

% % Post- hoc analysis for experiment factor:- Intagration factor is statiscaly sagnificent 
diff12=b(6);
diff13=b(7);
diff14=b(8);
diff23=b(7)-b(6);
diff24=b(8)-b(6);
diff34=b(8)-b(7);
% 
Sb2=mdl.CoefficientCovariance;
Sdiff12=sqrt(Sb2(6,6));
Sdiff13=sqrt(Sb2(7,7));
Sdiff14=sqrt(Sb2(8,8));
Sdiff23=sqrt(Sb2(7,7)+Sb2(6,6)-2*Sb2(6,7));
Sdiff24=sqrt(Sb2(8,8)+Sb2(6,6)-2*Sb2(6,8));
Sdiff34=sqrt(Sb2(8,8)+Sb2(7,7)-2*Sb2(7,8));

p12=2*(1-tcdf(abs(diff12/Sdiff12),nu_res));
p13=2*(1-tcdf(abs(diff13/Sdiff13),nu_res));
p14=2*(1-tcdf(abs(diff14/Sdiff14),nu_res));

p23=2*(1-tcdf(abs(diff23/Sdiff23),nu_res));
p24=2*(1-tcdf(abs(diff24/Sdiff24),nu_res));
p34=2*(1-tcdf(abs(diff34/Sdiff34),nu_res));

% Perform Bonferroni correction
p_value=[p12,p13,p14,p23,p24,p34];
num_comparisons = length(p_value);
alpha = 0.05;
bonferroni_alpha = alpha / num_comparisons;
[p_values_sorted, sorted_indices] = sort(p_value);
significant_indices = sorted_indices(p_values_sorted < bonferroni_alpha);
disp('Significant indices:');
p_value(significant_indices)
disp(['Bonferroni corrected alpha: ', num2str(bonferroni_alpha)]);
%%
clear all;
clc;
%% Figure 3 in the paper:

data(:,1)={[254,360,272,246,272,272,220];[82,110,100];[132,132,76,136,118];[68,102,84,140,96];[122,108,92,68]};% solder cell3 :albomin 800 cell4 albomin 735 cell5 albomin 670%21_10

figure('position',[100 100 400 550]);
AxPos=[0.2 0.45 0.7 0.5];
axes('position',AxPos);
median_data(1)=median(cell2mat(data(1,1)));
s_d=size(data);
for i=2:s_d(1)
    mean_data(i-1)=100.*mean(cell2mat(data(i,1)))./median_data(1);
    std_data(i-1)=std(100.*cell2mat(data(i,1))./median_data(1));
end
b=bar(1:s_d(1)-1,mean_data);
b.FaceColor = 'flat';

b.CData(1,:) =[0.3010 0.7450 0.9330]; % suturing group color
b.CData(2,:) =[40,106,70]./255; % albumin800 group color
b.CData(3,:) =[48,132,128]./255; % albumin 735 group color
b.CData(4,:) =[56,148,180]./255; % albumin 670 group color


hold on

for j=2:s_d(1)
p(1)=line([j-1.1,j-1.1],[mean_data(j-1)-std_data(j-1)  mean_data(j-1)+std_data(j-1)]);
p(1).Color=[ 0 0 0];
p(1).LineWidth=3;
end
marker={'o';'>';'>';'>'};
for n=2:s_d(1)
    n_d=length(cell2mat(data(n,1)));
    plot((n-1)*ones(1,n_d)+randi([-5 5],1,n_d)./100,100.*(cell2mat(data(n,1))./median_data(1)),'Marker',marker{n-1},'MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none','LineWidth',1)
hold on 
end
text(-1,-5,'Normalized Burst Pressure  [%]','FontSize',18,'FontName','Times new Roman','Rotation',90)
xcond={'Sutures','800 mg/ml Albumin','735 mg/ml Albumin','670 mg/ml Albumin'};
set(gca,'xtick',1:s_d(1)-1,'xticklabel',[],'xlim',[0 s_d(1)],'ylim',[0 100],'Fontsize',16,'FontName','Times new Roman','Box','off')
text(1,-30, 'Sutures','FontSize',16,'Rotation',90,'FontName','Times new Roman')
text(2,-70, '800 mg/ml Albumin','FontSize',16,'Rotation',90,'FontName','Times new Roman')
text(3,-70, '735 mg/ml Albumin','FontSize',16,'Rotation',90,'FontName','Times new Roman')
text(4,-70, '670 mg/ml Albumin','FontSize',16,'Rotation',90,'FontName','Times new Roman')
text(0.8,-80,'Experimntal condition','FontSize',18,'FontName','Times new Roman')
%% Statistical analysis for figure 3 dataset:
median_IT=median(cell2mat(data(1,1)));
allM_data=[];
con=[];
for i=2:s_d(1)
 allM_data=[allM_data;100*cell2mat(data(i,1))'./median_IT];
 con=[con;i*ones(length(cell2mat(data(i,1))),1)];
end

[p,tbl,stats,terms]=anovan(allM_data,con,'model',1,'varnames',{'condition'});
[results,~,~,gnames] = multcompare(stats);

%%
clear all;
clc;
%% Figure 4 in the paper:

data(:,1)={[152,166,166,172,150];[94,140,170,175,134];[164,76,40,80,88];[120,106,176];[106,136,152]};% solder cell3 :ICG0.6 cell4 ICG0.45 cell5 ICG 0.310_11mean_IT=median(cell2mat(data(1,1)));


figure('position',[100 100 400 600]);
AxPos=[0.2 0.4 0.7 0.5];
axes('position',AxPos);

median_data(1)=median(cell2mat(data(1,1)));
s_d=size(data);

for i=2:s_d(1)
    mean_data(i-1)=100.*mean(cell2mat(data(i,1)))./median_data(1);
    std_data(i-1)=std(100.*cell2mat(data(i,1))./median_data(1));
end
b=bar(1:s_d(1)-1,mean_data);
b.FaceColor = 'flat';

b.CData(1,:) =[0.3010 0.7450 0.9330]; % suturing color
b.CData(2,:) =[40,106,70]./255; % albumin800 ICG0.6
b.CData(3,:) =[57,151,100]./255; % ICG0.45
b.CData(4,:) =[113,201,153]./255; % ICG0.3

hold on

for j=2:s_d(1)
p(1)=line([j-1.1,j-1.1],[mean_data(j-1)-std_data(j-1)  mean_data(j-1)+std_data(j-1)]);
p(1).Color=[ 0 0 0];
p(1).LineWidth=3;
end
marker={'o';'>';'>';'>'};
for n=2:s_d(1)
    n_d=length(cell2mat(data(n,1)));
    plot((n-1)*ones(1,n_d)+randi([-5 5],1,n_d)./100,100.*(cell2mat(data(n,1))./median_data(1)),'Marker',marker{n-1},'MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none','LineWidth',1)
hold on 
end

ylabel('Normalized Burst Pressure  [%]','FontSize',18,'FontName','Times ')
xcond={'Sutures','ICG0.6','ICG0.45','ICG0.3'}
set(gca,'xtick',1:s_d(1)-1,'xticklabel',[],'xlim',[0 s_d(1)],'ylim',[0 120],'FontSize',16,'FontName','Times new Roman','Box','off')
text(1,-35, 'Sutures','FontSize',16,'Rotation',90,'FontName','Times new Roman')
text(2,-65, '0.6 mg/ml ICG','FontSize',16,'Rotation',90,'FontName','Times new Roman')
text(3,-70, '0.45 mg/ml ICG','FontSize',16,'Rotation',90,'FontName','Times new Roman')
text(4,-65, '0.3 mg/ml ICG','FontSize',16,'Rotation',90,'FontName','Times new Roman')
text(0.8,-80,'Experimntal condition','FontSize',18,'FontName','Times new Roman')
%% Statistical analysis for figure 4 dataset:

median_IT=median(cell2mat(data(1,1)));
allM_data=[];
con=[];
for i=2:s_d(1)
 allM_data=[allM_data;100*cell2mat(data(i,1))'./median_IT];
 con=[con;i*ones(length(cell2mat(data(i,1))),1)];
end

[p,tbl,stats,terms]=anovan(allM_data,con,'model',1,'varnames',{'condition'});
[results,~,~,gnames] = multcompare(stats);

%%
clear all;
clc;
%% Figure 5 in the paper:

data(:,1)={[284,304,356,304];[68,108,122,118,106];[188,168,180,210]};%first day of HM
data(:,2)={[292,238,210,286,244];[136,108,82,114,112];[136,104,142,152,216]}; %second day of HM
data(:,3)={[170,162,172,162];[];[60,120,110,136]}; % cadaver

s_d=size(data); % get the data size for the loop 

figure('position',[100 100 300 550]);
AxPos=[0.25 0.35 0.7 0.6];
axes('position',AxPos);


for i=1:s_d(2) % runs for the different experimantal days
    median_IT(i)=median(cell2mat(data(1,i))); % the median value for the undamaged tissue group
    
    % The data analysis for the suturing group:
    suture=100*cell2mat(data(2,i))./median_IT(i);
    meanS(i)=mean(suture);
    stdS=std(suture);
    b=bar(i-0.125,meanS(i),0.25,'FaceColor','flat');  % Plot the mean value of the suturing group
    hold on
    b.CData =[0.3010 0.7450 0.9330]; % suturing group color
    line([i-0.125 i-0.125],[meanS(i)-stdS meanS(i)+stdS],'color',[0 0 0],'LineWidth',3) % Plot the standard deviation value of the suturing group
   
    for j=1:length(suture) % Plot each indevidual value of the suturing group
    plot(i+rand(1)/12-0.125,suture(j),'Marker','o','MarkerEdgeColor','k','MarkerSize',8,'LineWidth',1) 
    end

    %The data analysis for the soldering group
    for n=3:s_d(1)
        if ~isempty(cell2mat(data(n,i)))
            solder(n-2,:)=100*cell2mat(data(n,i))./median_IT(i);
            meanRLTS(n-2,i)=mean(solder(n-2,:));
            stdSo=std(solder(n-2,:));
            b2=bar(i+0.125,meanRLTS(n-2,i),0.25,'FaceColor','flat');% Plot the mean value of the soldering group
            b2.CData =[57,151,100]./255;  % soldering (albumin 800 and ICG 0.45) color 
            line([i+0.125 i+0.125],[meanRLTS(i)-stdSo meanRLTS(i)+stdSo],'color',[0 0 0],'LineWidth',3) %Plot the standard deviation value of the soldering group
            for j=1:length(solder) % Plot each indevidual value of the suturing group
            plot(i+rand(1)/10+0.125,solder(n-2,j),'Marker','>','MarkerEdgeColor','k','MarkerSize',8,'LineWidth',1) 
            end
            clear solder
        end
    end

end
line([2.5,2.5],[0 120],'LineStyle','--','color','k','LineWidth',3)
ylabel('Normalized Burst Pressure  [%]','FontName','Times new Roman','FontSize',18)
set(gca,'xtick',1:3,'xticklabel',[],'xlim',[0.5 3.5],'ylim',[0 120],'FontSize',16,'FontName','Times new Roman','Box','off')
text(0.9,-45, 'HM at BGU','FontSize',16,'Rotation',90,'FontName','Times new Roman')  
text(1.15,-35, 'lab','FontSize',16,'Rotation',90,'FontName','Times new Roman') 
text(1.9,-45, 'HM at Lahav','FontSize',16,'Rotation',90,'FontName','Times new Roman')
text(2.15,-35, 'C.R.O','FontSize',16,'Rotation',90,'FontName','Times new Roman')
text(3,-40, 'Cadaver','FontSize',16,'Rotation',90,'FontName','Times new Roman')
text(0.5,-55, 'Experiment protocol','FontSize',16,'FontName','Times new Roman')

%% Statistical analysis for figure 5 dataset:

s_d=size(data);
allM_data=[];
for i=1:s_d(2)-1
    median_IT(i)=median(cell2mat(data(1,i)));
    suture=100*cell2mat(data(2,i))./median_IT(i);
    solder=100*cell2mat(data(3,i))./median_IT(i);
    M_data(:,1)=[suture';solder'];
    M_data(:,2)=[-1*ones(length(suture),1);ones(length(solder),1)];
   M_data(:,3)=i.*ones(length(suture)+length(solder),1);
   allM_data=[allM_data; M_data];
   clear M_data
end

for con=2
    C(allM_data(:,3)==con,con-1)=1;
    C(allM_data(:,3)==1,con-1)=-1;
end
tab=table(allM_data(:,1),allM_data(:,2),C(:,1),'VariableNames',{'y','method','exp2'});
[p,tbl,stats,terms]=anovan(allM_data(:,1),{allM_data(:,2),allM_data(:,3)},'model',2,'varnames',{'method','exp'});
%[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);