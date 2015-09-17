clc
clear
% theta(1:12)=1;alpha(1:12)=0.8;%%test
%%test--att
% theta=[7.86484464247420,10.9256019189366,12.1639000338355,10.9873719336213,-0.208406840480289,6.36056234583259,9.22535631659349,4.12304325661133,6.19609773486933,11.0675249646027,8.65058054206778,7.75996773151805];
% alpha=[7.12048698166934,7.27345437750985,9.14607740143361,7.18443870429327,0.367035470826866,1.38776393513730,8.01847057947578,0.391638658994368,1.50557516560250,8.43804934276927,3.09135615475917,2.69970015360378];
%%test--rest
% theta=[12.7751974845368,12.7873935565209,12.8962042319162,19.1385031468905,12.7100661834402,12.7459554655817,34.1315627259312,14.3064915024945,13.6681888999645,14.2612456514516,13.4717222197866,14.0619653548311];
% alpha=[7.40583912968185,6.92579144986660,7.36381238428012,14.8612708129992,6.92286899295752,8.07671026769833,27.5326382587131,7.29869436822047,6.33531563960417,8.23798840487257,6.32831416267540,7.17045838303781];
%% this part aims to get the threshold----alpha , theta, and their ratio  (this part should be at starting part)
load('zj_rest2');
% load('rest_2.mat');
delta_rest=rest_2(:,1:3);
theta_rest=rest_2(:,4:7);
alpha_rest=rest_2(:,8:12);
beta_rest=rest_2(:,13:23);
aa=0;bb=0;cc=0;dd=0; ave_theta_rest=0;ave_alpha_rest=0;ratio_rest=0;
for i=1:12
   for j=1:length(theta_rest)/12
       aa(i,j)=mean(delta_rest(j*12-12+i,1:3));
     bb(i,j)=mean(theta_rest(j*12-12+i,1:4)) ;
      cc(i,j)=mean(alpha_rest(j*12-12+i,1:5)) ;
      dd(i,j)=mean(beta_rest(j*12-12+i,1:11));
   end
   ave_delta_rest(i)=mean(aa(i,:));
   ave_beta_rest(i)=mean(dd(i,:));
   ave_theta_rest(i)=mean(bb(i,:));
   ave_alpha_rest(i)=mean(cc(i,:));
   
%    ratio_rest(i)=ave_theta_rest(i)/ave_alpha_rest(i);%%%ratio
p_delta_rest(i)=10^(ave_delta_rest(i)/10);
p_beta_rest(i)=10^(ave_beta_rest(i)/10);
p_theta_rest(i)=10^(ave_theta_rest(i)/10);
    p_alpha_rest(i)=10^(ave_alpha_rest(i)/10);
    p_ratio_thetaoveralpha_rest(i)=p_theta_rest(i)/p_alpha_rest(i);
    p_ratio_betaoveralpha_rest(i)=p_beta_rest(i)/p_alpha_rest(i);
    p_ratio_thetaalphaoverbeta_rest(i)=(p_theta_rest(i)+p_alpha_rest(i))/p_beta_rest(i);
end
%% normalize
for k=1:12
summ(k)=p_delta_rest(k)+p_alpha_rest(k)+p_theta_rest(k)+p_beta_rest(k);
% summ(k)=ave_delta_rest(k)+ave_alpha_rest(k)+ave_theta_rest(k)+ave_beta_rest(k);

norm_delta_rest(k)=p_delta_rest(k)/summ(k);
norm_alpha_rest(k)=p_alpha_rest(k)/summ(k);
norm_beta_rest(k)=p_beta_rest(k)/summ(k);
norm_theta_rest(k)=p_theta_rest(k)/summ(k);

% norm_delta_rest(k)=ave_delta_rest(k)/summ(k);
% norm_alpha_rest(k)=ave_alpha_rest(k)/summ(k);
% norm_beta_rest(k)=ave_beta_rest(k)/summ(k);
% norm_theta_rest(k)=ave_theta_rest(k)/summ(k);
end
%%
load('zj_att2.mat');


delta_att=att_2(:,1:3);
theta_att=att_2(:,4:7);
alpha_att=att_2(:,8:12);
beta_att=att_2(:,13:23);
aa=0;bb=0;cc=0;dd=0; 
ave_theta_att=0;ave_alpha_att=0;ratio_att=0;

for i=1:12
   for j=1:length(theta_att)/12
       aa(i,j)=mean(delta_att(j*12-12+i,1:3));
     bb(i,j)=mean(theta_att(j*12-12+i,1:4)) ;
      cc(i,j)=mean(alpha_att(j*12-12+i,1:5)) ;
      dd(i,j)=mean(beta_att(j*12-12+i,1:11));
   end
   ave_delta_att(i)=mean(aa(i,:));
   ave_beta_att(i)=mean(dd(i,:));
   ave_theta_att(i)=mean(bb(i,:));
   ave_alpha_att(i)=mean(cc(i,:));
   
%    ratio_att(i)=ave_theta_att(i)/ave_alpha_att(i);%%%ratio
p_delta_att(i)=10^(ave_delta_att(i)/10);
p_beta_att(i)=10^(ave_beta_att(i)/10);
p_theta_att(i)=10^(ave_theta_att(i)/10);
    p_alpha_att(i)=10^(ave_alpha_att(i)/10);
      p_ratio_thetaoveralpha_att(i)=p_theta_att(i)/p_alpha_att(i);
    p_ratio_betaoveralpha_att(i)=p_beta_att(i)/p_alpha_att(i);
    p_ratio_thetaalphaoverbeta_att(i)=(p_theta_att(i)+p_alpha_att(i))/p_beta_att(i);
end

for l=1:12
summ(l)=p_delta_att(l)+p_alpha_att(l)+p_theta_att(l)+p_beta_att(l);

norm_delta_att(l)=p_delta_att(l)/summ(l);
norm_alpha_att(l)=p_alpha_att(l)/summ(l);
norm_beta_att(l)=p_beta_att(k)/summ(l);
norm_theta_att(l)=p_theta_att(l)/summ(l);
% norm_ratio_att(l)=norm_theta_att(l)/norm_alpha_att(l);
end

% %% this part  aims to determine user's attention level
% % 
% % theta=mean(freq_data(:,4:7)');
% % alpha=mean(freq_data(:,8:12)');
% % theta=mean(rest(:,4:7)');
% % alpha=mean(rest(:,8:12)');
% theta_rate=0;alpha_rate=0;ratio_rate=0;ratio=0;
% for k=1:12
% theta_rate(k)=(theta(k)-ave_theta_rest(k))/(ave_theta_att(k)-ave_theta_rest(k));
% if theta_rate(k)>1 
%     theta_rate(k)=1;
% end
%     if theta_rate(k)<0
%         theta_rate(k)=0;
%     end
% 
% 
% alpha_rate(k)=(alpha(k)-ave_alpha_rest(k))/(ave_alpha_att(k)-ave_alpha_rest(k));
% if alpha_rate(k)>1 
%     alpha_rate(k)=1;
% end
%     if alpha_rate(k)<0
%         alpha_rate(k)=0;
%     end
% 
% 
% ratio(k)=theta(k)/alpha(k);
% ratio_rate(k)=(ratio(k)-ratio_rest(k))/(ratio_att(k)-ratio_rest(k));
% if ratio_rate(k)>1 
%     ratio_rate(k)=1;
% end
%     if ratio_rate(k)<0
%         ratio_rate=0;
%     end
% end
% 
% att_lvl=sum(ratio_rate/3+alpha_rate/3+theta_rate/3)/12;%%the weight may change
% 
% if att_lvl<0.2
%     disp('1')
% end 
% 
% if att_lvl>0.2&att_lvl<0.4
%     disp('2')
% end 
% 
% if att_lvl>0.4&att_lvl<0.6
%     disp('3')
% end 
% 
% if att_lvl>0.6&att_lvl<0.8
%     disp('4')
% end 
% 
% if att_lvl>0.8
%     disp('5')
% end 


