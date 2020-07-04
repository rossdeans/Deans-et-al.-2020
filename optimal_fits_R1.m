%script to fit optimal model to A, Kleaf and gs data

%Now includes elastic modulus epsilon, to take into account dilution
%effects on osmotic pressure
%now also includes modelling for how A, gs and Kleaf should vary with dw
%and psi_s given a Vcmax (which is calculated from a given A at dw=0.015
%and psi_s=0).

%also does relative costs now relative to total costs, not A


raw_data_cal=readtable('cal_set.xls')
raw_data_test=readtable('test_set.xls')
%%
dw=1.5e-2;
psi_s=0;
Ca=400;
g_star=48.4;
epsilon=10; %elastic modulus in MPa

%convert to structure for manipulations
data_cal=table2struct(raw_data_cal);
data_test=table2struct(raw_data_test);

%to calculate variables
for i=1:size(data_cal,1)
    %data_cal(i).chi=data_cal(i).gs*data_cal(i).Kleaf*1e-3/(data_cal(i).Kleaf*1e-3*(data_cal(i).TLP+psi_s)-data_cal(i).gs*dw);
    data_cal(i).chi=data_cal(i).gs*data_cal(i).Kleaf*1e-3*(1+data_cal(i).TLP/epsilon)/(data_cal(i).Kleaf*1e-3*(data_cal(i).TLP+psi_s)-data_cal(i).gs*dw);
    
    
    %data_cal(i).mu_chi=1.6*(data_cal(i).A)^2/((Ca-g_star)*(data_cal(i).chi)^2*(data_cal(i).TLP+psi_s));
    %data_cal(i).mu_K=1.6*dw*(data_cal(i).A)^2/((Ca-g_star)*1e-6*(data_cal(i).Kleaf)^2*(data_cal(i).TLP+psi_s));
    %data_cal(i).mu_piv=1.6*data_cal(i).A*(1e-3*data_cal(i).Kleaf+data_cal(i).chi*dw)*((Ca-g_star)*data_cal(i).chi*1e-3*data_cal(i).Kleaf*(data_cal(i).TLP+psi_s)-1.6*data_cal(i).A*(1e-3*data_cal(i).Kleaf+data_cal(i).chi*dw))/((Ca-g_star)*(data_cal(i).chi)^2*(1e-3*data_cal(i).Kleaf)^2*(data_cal(i).TLP+psi_s)^3);
    
    data_cal(i).mu_chi=1.6*(1+data_cal(i).TLP/epsilon)*(data_cal(i).A)^2/((Ca-g_star)*(data_cal(i).chi)^2*(data_cal(i).TLP+psi_s));
    data_cal(i).mu_K=1.6*dw*(data_cal(i).A)^2/((Ca-g_star)*1e-6*(data_cal(i).Kleaf)^2*(data_cal(i).TLP+psi_s));
    data_cal(i).mu_piv=1.6*data_cal(i).A*(1e-3*data_cal(i).Kleaf*(1-psi_s/epsilon)+data_cal(i).chi*dw)*((Ca-g_star)*data_cal(i).chi*1e-3*data_cal(i).Kleaf*(data_cal(i).TLP+psi_s)-1.6*data_cal(i).A*(1e-3*data_cal(i).Kleaf*(1+data_cal(i).TLP/epsilon)+data_cal(i).chi*dw))/((Ca-g_star)*(data_cal(i).chi)^2*(1e-3*data_cal(i).Kleaf)^2*(data_cal(i).TLP+psi_s)^3);
    
    
    %relative costs, relative to total costs, as a percentage
    %data_cal(i).p_chi=100*1e-3*data_cal(i).Kleaf/(1e-3*data_cal(i).Kleaf+data_cal(i).chi*dw)*(data_cal(i).TLP+psi_s)/(2*data_cal(i).TLP+psi_s);
    %data_cal(i).p_K=100*data_cal(i).chi*dw/(1e-3*data_cal(i).Kleaf+data_cal(i).chi*dw)*(data_cal(i).TLP+psi_s)/(2*data_cal(i).TLP+psi_s);
    %data_cal(i).p_pi=100*data_cal(i).TLP/(2*data_cal(i).TLP+psi_s);
    
    data_cal(i).p_chi=100*1e-3*data_cal(i).Kleaf*(data_cal(i).TLP+psi_s)/(data_cal(i).TLP*(1e-3*data_cal(i).Kleaf*(1+(1-psi_s/epsilon)/(1+data_cal(i).TLP/epsilon))+2*data_cal(i).chi*dw/(1+data_cal(i).TLP/epsilon))+psi_s*(1e-3*data_cal(i).Kleaf+data_cal(i).chi*dw/(1+data_cal(i).TLP/epsilon)));
    data_cal(i).p_K=100*data_cal(i).chi*dw*(data_cal(i).TLP+psi_s)/(1+data_cal(i).TLP/epsilon)/(data_cal(i).TLP*(1e-3*data_cal(i).Kleaf*(1+(1-psi_s/epsilon)/(1+data_cal(i).TLP/epsilon))+2*data_cal(i).chi*dw/(1+data_cal(i).TLP/epsilon))+psi_s*(1e-3*data_cal(i).Kleaf+data_cal(i).chi*dw/(1+data_cal(i).TLP/epsilon)));
    data_cal(i).p_pi=100/((1+psi_s/data_cal(i).TLP)*((1e-3*data_cal(i).Kleaf*(1+data_cal(i).TLP/epsilon)+data_cal(i).chi*dw)/(1e-3*data_cal(i).Kleaf*(1-psi_s/epsilon)+data_cal(i).chi*dw))+1);

end

%calculate mean cost parameters
mu_chi_mean=mean([data_cal.mu_chi]);
mu_K_mean=mean([data_cal.mu_K]);
mu_piv_mean=mean([data_cal.mu_piv]);

%calculate mean TLP
TLP_mean=mean([data_cal.TLP]);

%find TLP calculated from mus
f=funct_epsilon_fminsearch(Ca,g_star,dw,mu_K_mean,mu_chi_mean,mu_piv_mean,psi_s,epsilon);
TLP_calc=fminsearch(f,-psi_s+2);

%calculate predicted psi
%psi_calc=psi_s-(TLP_calc+psi_s)/((mu_chi_mean*(1+TLP_calc/epsilon)/(mu_K_mean*dw))^0.5+1);
%use delta psi instead
D_psi_calc=(TLP_calc+psi_s)/((mu_chi_mean*(1+TLP_calc/epsilon)/(mu_K_mean*dw))^0.5+1);


%test whether TLPs are significantly different from calculated TLP
[h_TLP,p_TLP]=ttest([data_cal.TLP],TLP_calc);


%%
%to compare models (calibration data and test data) against data
%calculate predicted chi, Kleaf and gs for calibration data
for i=1:size(data_cal,1)
    %data_cal(i).chi_calc=(1.6/(mu_chi_mean*(Ca-g_star)*(TLP_calc+psi_s)))^0.5*data_cal(i).A;
    %data_cal(i).Kleaf_calc=1e3*(1.6*dw/(mu_K_mean*(Ca-g_star)*(TLP_calc+psi_s)))^0.5*data_cal(i).A;
    %data_cal(i).gs_calc=data_cal(i).chi_calc*1e-3*data_cal(i).Kleaf_calc*(TLP_calc+psi_s)/(1e-3*data_cal(i).Kleaf_calc+data_cal(i).chi_calc*dw);
    %data_cal(i).psi_calc=psi_s-(data_cal(i).gs_calc*dw/(1e-3*data_cal(i).Kleaf_calc));
    
    data_cal(i).chi_calc=(1.6*(1+TLP_calc/epsilon)/(mu_chi_mean*(Ca-g_star)*(TLP_calc+psi_s)))^0.5*data_cal(i).A;
    data_cal(i).Kleaf_calc=1e3*(1.6*dw/(mu_K_mean*(Ca-g_star)*(TLP_calc+psi_s)))^0.5*data_cal(i).A;
    data_cal(i).gs_calc=data_cal(i).chi_calc*1e-3*data_cal(i).Kleaf_calc*(TLP_calc+psi_s)/(1e-3*data_cal(i).Kleaf_calc*(1+TLP_calc/epsilon)+data_cal(i).chi_calc*dw);
    %data_cal(i).psi_calc=psi_s-(data_cal(i).gs_calc*dw/(1e-3*data_cal(i).Kleaf_calc));
    %use delta psi instead
    data_cal(i).D_psi_calc=(data_cal(i).gs_calc*dw/(1e-3*data_cal(i).Kleaf_calc));
    
    
    %calculate psi (observed LWP from gs and Kleaf)
    %data_cal(i).psi=psi_s-(data_cal(i).gs*dw/(1e-3*data_cal(i).Kleaf));
    %use delta psi instead
    data_cal(i).D_psi=(data_cal(i).gs*dw/(1e-3*data_cal(i).Kleaf));
end

%calculate predicted chi, Kleaf and gs for test data
for i=1:size(data_test,1)
    %data_test(i).chi_calc=(1.6/(mu_chi_mean*(Ca-g_star)*(TLP_calc+psi_s)))^0.5*data_test(i).A;
    %data_test(i).Kleaf_calc=1e3*(1.6*dw/(mu_K_mean*(Ca-g_star)*(TLP_calc+psi_s)))^0.5*data_test(i).A;
    %data_test(i).gs_calc=data_test(i).chi_calc*1e-3*data_test(i).Kleaf_calc*(TLP_calc+psi_s)/(1e-3*data_test(i).Kleaf_calc+data_test(i).chi_calc*dw);
    %data_test(i).psi_calc=psi_s-(data_test(i).gs_calc*dw/(1e-3*data_test(i).Kleaf_calc));
    
    data_test(i).chi_calc=(1.6*(1+TLP_calc/epsilon)/(mu_chi_mean*(Ca-g_star)*(TLP_calc+psi_s)))^0.5*data_test(i).A;
    data_test(i).Kleaf_calc=1e3*(1.6*dw/(mu_K_mean*(Ca-g_star)*(TLP_calc+psi_s)))^0.5*data_test(i).A;
    data_test(i).gs_calc=data_test(i).chi_calc*1e-3*data_test(i).Kleaf_calc*(TLP_calc+psi_s)/(1e-3*data_test(i).Kleaf_calc*(1+TLP_calc/epsilon)+data_test(i).chi_calc*dw);
    %data_test(i).psi_calc=psi_s-(data_test(i).gs_calc*dw/(1e-3*data_test(i).Kleaf_calc));
    %use D_psi instead
    data_test(i).D_psi_calc=(data_test(i).gs_calc*dw/(1e-3*data_test(i).Kleaf_calc));
    
    
    %calculate psi (observed LWP from gs and Kleaf)
    %data_test(i).psi=psi_s-(data_test(i).gs*dw/(1e-3*data_test(i).Kleaf));
    %use D_psi instead
    data_test(i).D_psi=(data_test(i).gs*dw/(1e-3*data_test(i).Kleaf));
    
    
    %calculate mu for Cowan-Farquhar model
    data_test(i).mu_E=1.6*(data_test(i).A)^2/(dw*(Ca-g_star)*(data_test(i).gs)^2);
    
    %calculate k for each species (following Reviewer 2 comments
    data_test(i).k=data_test(i).A*data_test(i).gs/(data_test(i).gs*(Ca-g_star)-1.6*data_test(i).A);
    
    
end

%find mean mu_E for Cowan-Farquhar model
mu_E_mean=mean([data_test.mu_E]);

%find mean k
mean_k=mean([data_test.k]);

%calculate residuals etc for R^2 for full data
SSE_cal_Kleaf=sum(([data_cal.Kleaf_calc]-[data_cal.Kleaf]).^2);
SST_cal_Kleaf=sum(([data_cal.Kleaf]-mean([data_cal.Kleaf])).^2);
Rsquared_cal_Kleaf=1-SSE_cal_Kleaf/SST_cal_Kleaf;

SSE_cal_gs=sum(([data_cal.gs_calc]-[data_cal.gs]).^2);
SST_cal_gs=sum(([data_cal.gs]-mean([data_cal.gs])).^2);
Rsquared_cal_gs=1-SSE_cal_gs/SST_cal_gs;

%calculate residuals etc for R^2 for partial data
SSE_test_Kleaf=sum(([data_test.Kleaf_calc]-[data_test.Kleaf]).^2);
SST_test_Kleaf=sum(([data_test.Kleaf]-mean([data_test.Kleaf])).^2);
Rsquared_test_Kleaf=1-SSE_test_Kleaf/SST_test_Kleaf;

SSE_test_gs=sum(([data_test.gs_calc]-[data_test.gs]).^2);
SST_test_gs=sum(([data_test.gs]-mean([data_test.gs])).^2);
Rsquared_test_gs=1-SSE_test_gs/SST_test_gs;



%%
%to model the case for soil drying and different vpd

%pick an A at dw=0.015 and psi_s=0 and calculate the equivalent V, then
%simulate optimal solution for varying dw and psi_s

%initial A under the usual conditions
A0=15;

%set ranges of simulations
dw_range=(0.1:0.1:3)*1e-2; %old range was a vpd of 0.1 to 5
psi_s_range=-(0:0.1:2); %old range was from 0 to 3

I=length(dw_range);
J=length(psi_s_range);

%declare vectors to hold solutions
TLP_dw=zeros(I,1);
chi_dw=zeros(I,1);
Kleaf_dw=zeros(I,1);
A_dw=zeros(I,1);
gs_dw=zeros(I,1);
%psi_dw=zeros(I,1);
D_psi_dw=zeros(I,1); %use delta psi now
A_CF=zeros(I,1);
gs_CF=zeros(I,1);

TLP_psi=zeros(J,1);
chi_psi=zeros(J,1);
Kleaf_psi=zeros(J,1);
A_psi=zeros(J,1);
gs_psi=zeros(J,1);
%psi_psi=zeros(J,1); %use delta psi now
D_psi_psi=zeros(J,1);

%calculate initial optimal values
%find initial TLP
f=funct_epsilon_fminsearch(Ca,g_star,dw,mu_K_mean,mu_chi_mean,mu_piv_mean,psi_s,epsilon);
TLP0=fminsearch(f,-psi_s+2);

%calculate initial chi
%chi0=(1.6/(mu_chi_mean*(Ca-g_star)*(TLP0+psi_s)))^0.5*A0;
chi0=(1.6*(1+TLP0/epsilon)/(mu_chi_mean*(Ca-g_star)*(TLP0+psi_s)))^0.5*A0;

%calculate initial Kleaf
%Kleaf0=1e3*(1.6*dw/(mu_K_mean*(Ca-g_star)*(TLP0+psi_s)))^0.5*A0;
Kleaf0=1e3*(1.6*dw/(mu_K_mean*(Ca-g_star)*(TLP0+psi_s)))^0.5*A0;

%calculate effective V (now k in model)
%V=(chi0*1e-3*Kleaf0*(TLP0+psi_s)*A0)/((Ca-g_star)*chi0*1e-3*Kleaf0*(TLP0+psi_s)-1.6*(1e-3*Kleaf0+chi0*dw)*A0);
V=(chi0*1e-3*Kleaf0*(TLP0+psi_s)*A0)/((Ca-g_star)*chi0*1e-3*Kleaf0*(TLP0+psi_s)-1.6*(1e-3*Kleaf0*(1+TLP0/epsilon)+chi0*dw)*A0);

%calculate initial gs for Cowan-Farquhar model
gs0=(1.6/(mu_E_mean*dw*(Ca-g_star)))^0.5*A0;

%calculate k for Cowan-Farquhar model
k=gs0*A0/(gs0*(Ca-g_star)-1.6*A0);

%do calculations for variation in dw
for i=1:I
    
    %calculate TLP
    f=funct_epsilon_fminsearch(Ca,g_star,dw_range(i),mu_K_mean,mu_chi_mean,mu_piv_mean,psi_s,epsilon);
    TLP_dw(i)=fminsearch(f,-psi_s+2);
    
    %calculate other variables
    %chi_dw(i)=1.6*V/(TLP_dw(i)+psi_s)*(((Ca-g_star)*(TLP_dw(i)+psi_s)/(1.6*mu_chi_mean))^0.5-(1+(mu_K_mean*dw_range(i)/mu_chi_mean)^0.5));
    %Kleaf_dw(i)=1e3*chi_dw(i)*(mu_chi_mean*dw_range(i)/mu_K_mean)^0.5;
    %A_dw(i)=V*(Ca-g_star)*chi_dw(i)*1e-3*Kleaf_dw(i)*(TLP_dw(i)+psi_s)/(1.6*V*(1e-3*Kleaf_dw(i)+chi_dw(i)*dw_range(i))+chi_dw(i)*1e-3*Kleaf_dw(i)*(TLP_dw(i)+psi_s));
    %gs_dw(i)=chi_dw(i)*1e-3*Kleaf_dw(i)*(TLP_dw(i)+psi_s)/(1e-3*Kleaf_dw(i)+chi_dw(i)*dw_range(i));
    
    chi_dw(i)=V*(1.6*(1+TLP_dw(i)/epsilon)/(TLP_dw(i)+psi_s))^0.5*(((Ca-g_star)/mu_chi_mean)^0.5-(1.6/(TLP_dw(i)+psi_s))^0.5*((1+TLP_dw(i)/epsilon)^0.5+(mu_K_mean*dw_range(i)/mu_chi_mean)^0.5));
    Kleaf_dw(i)=1e3*V*(1.6*mu_chi_mean*dw_range(i)/(mu_K_mean*(TLP_dw(i)+psi_s)))^0.5*(((Ca-g_star)/mu_chi_mean)^0.5-(1.6/(TLP_dw(i)+psi_s))^0.5*((1+TLP_dw(i)/epsilon)^0.5+(mu_K_mean*dw_range(i)/mu_chi_mean)^0.5));
    A_dw(i)=V*(Ca-g_star)*chi_dw(i)*1e-3*Kleaf_dw(i)*(TLP_dw(i)+psi_s)/(1.6*V*(1e-3*Kleaf_dw(i)*(1+TLP_dw(i)/epsilon)+chi_dw(i)*dw_range(i))+chi_dw(i)*1e-3*Kleaf_dw(i)*(TLP_dw(i)+psi_s));
    gs_dw(i)=chi_dw(i)*1e-3*Kleaf_dw(i)*(TLP_dw(i)+psi_s)/(1e-3*Kleaf_dw(i)*(1+TLP_dw(i)/epsilon)+chi_dw(i)*dw_range(i));
    %psi_dw(i)=psi_s-gs_dw(i)*dw_range(i)/(1e-3*Kleaf_dw(i));
    D_psi_dw(i)=gs_dw(i)*dw_range(i)/(1e-3*Kleaf_dw(i)); %now using delta psi
    
    gs_CF(i)=1.6*k*(((Ca-g_star)/(1.6*mu_E_mean*dw_range(i)))^0.5-1);
    A_CF(i)=k*gs_CF(i)*(Ca-g_star)/(1.6*k+gs_CF(i));
    
    
end

%do calculations for variation in psi_s
for j=1:J
    
    %calculate TLP
    f=funct_epsilon_fminsearch(Ca,g_star,dw,mu_K_mean,mu_chi_mean,mu_piv_mean,psi_s_range(j),epsilon);
    TLP_psi(j)=fminsearch(f,-psi_s_range(j)+2);
    
    %calculate other variables
    %chi_psi(j)=1.6*V/(TLP_psi(j)+psi_s_range(j))*(((Ca-g_star)*(TLP_psi(j)+psi_s_range(j))/(1.6*mu_chi_mean))^0.5-(1+(mu_K_mean*dw/mu_chi_mean)^0.5));
    %Kleaf_psi(j)=1e3*chi_psi(j)*(mu_chi_mean*dw/mu_K_mean)^0.5;
    %A_psi(j)=V*(Ca-g_star)*chi_psi(j)*1e-3*Kleaf_psi(j)*(TLP_psi(j)+psi_s_range(j))/(1.6*V*(1e-3*Kleaf_psi(j)+chi_psi(j)*dw)+chi_psi(j)*1e-3*Kleaf_psi(j)*(TLP_psi(j)+psi_s_range(j)));
    %gs_psi(j)=chi_psi(j)*1e-3*Kleaf_psi(j)*(TLP_psi(j)+psi_s_range(j))/(1e-3*Kleaf_psi(j)+chi_psi(j)*dw);
    
    chi_psi(j)=V*(1.6*(1+TLP_psi(j)/epsilon)/(TLP_psi(j)+psi_s_range(j)))^0.5*(((Ca-g_star)/mu_chi_mean)^0.5-(1.6/(TLP_psi(j)+psi_s_range(j)))^0.5*((1+TLP_psi(j)/epsilon)^0.5+(mu_K_mean*dw/mu_chi_mean)^0.5));
    Kleaf_psi(j)=1e3*V*(1.6*mu_chi_mean*dw/(mu_K_mean*(TLP_psi(j)+psi_s_range(j))))^0.5*(((Ca-g_star)/mu_chi_mean)^0.5-(1.6/(TLP_psi(j)+psi_s_range(j)))^0.5*((1+TLP_psi(j)/epsilon)^0.5+(mu_K_mean*dw/mu_chi_mean)^0.5));
    A_psi(j)=V*(Ca-g_star)*chi_psi(j)*1e-3*Kleaf_psi(j)*(TLP_psi(j)+psi_s_range(j))/(1.6*V*(1e-3*Kleaf_psi(j)*(1+TLP_psi(j)/epsilon)+chi_psi(j)*dw)+chi_psi(j)*1e-3*Kleaf_psi(j)*(TLP_psi(j)+psi_s_range(j)));
    gs_psi(j)=chi_psi(j)*1e-3*Kleaf_psi(j)*(TLP_psi(j)+psi_s_range(j))/(1e-3*Kleaf_psi(j)*(1+TLP_psi(j)/epsilon)+chi_psi(j)*dw);
    %psi_psi(j)=psi_s_range(j)-gs_psi(j)*dw/(1e-3*Kleaf_psi(j));
    D_psi_psi(j)=gs_psi(j)*dw/(1e-3*Kleaf_psi(j)); %now using delta psi
    
end


%%
%to plot the data and fits
fontsize=10;
markersize=5;
linewidth=1;

% to calculate modelled outputs for plotting
A=0:1:20;
chi=(1.6*(1+TLP_calc/epsilon)/(mu_chi_mean*(Ca-g_star)*(TLP_calc+psi_s)))^0.5*A;
Kleaf=1e3*(1.6*dw/(mu_K_mean*(Ca-g_star)*(TLP_calc+psi_s)))^0.5*A;
gs=chi*1e-3.*Kleaf.*(TLP_calc+psi_s)./(1e-3.*Kleaf.*(1+TLP_calc/epsilon)+chi.*dw);
D_psi=gs*dw./(1e-3*Kleaf);




%%
%Figure 1 model and data outputs vs A
f1=figure('Units','centimeters','Position',[1 1 17.3 14]);
set(f1,'Units','normalized');




%Fig 1a (gs vs A)
ax_1a=subplot(2,2,1)
plot(A,gs,'--k','LineWidth',linewidth)
hold on
plot([data_test.A],[data_test.gs],'ro','MarkerSize',markersize)
plot([data_cal.A],[data_cal.gs],'ko','MarkerSize',markersize)
xlabel('\itA\rm_{sat} (\mumol m^{-2} s^{-1})','Interpreter', 'tex','FontSize',fontsize)
ylabel('\itg\rm_{sat} (mol m^{-2} s^{-1})','Interpreter', 'tex','FontSize',fontsize)
box('on')
xl=[0 20];
yl=[0 0.71];
axis([xl yl])
text(xl(1)+0.05*(xl(2)-xl(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfa\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)

%Fig 1b (Kleaf vs A)
ax_1b=subplot(2,2,2)
plot(A,Kleaf,'--k','LineWidth',linewidth)
hold on
plot([data_test.A],[data_test.Kleaf],'ro','MarkerSize',markersize)
plot([data_cal.A],[data_cal.Kleaf],'ko','MarkerSize',markersize)
xlabel('\itA\rm_{sat} (\mumol m^{-2} s^{-1})','Interpreter', 'tex','FontSize',fontsize)
ylabel('\itK\rm_{leaf} (mmol m^{-2} s^{-1} MPa^{-1})','Interpreter', 'tex','FontSize',fontsize)
box('on')
xl=[0 20];
yl=[0 25];
axis([xl yl])
text(xl(1)+0.05*(xl(2)-xl(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfb\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)

%Fig 1c (psi vs A)
ax_1c=subplot(2,2,3)
D_psi_A=fitlm([data_test.A],[data_test.D_psi])
[h_psi,p_psi]=ttest([data_test.D_psi],D_psi_calc)
plot(A,D_psi,'--k','LineWidth',linewidth)
hold on
plot([data_test.A],[data_test.D_psi],'ro','MarkerSize',markersize)
plot([data_cal.A],[data_cal.D_psi],'ko','MarkerSize',markersize)
xlabel('\itA\rm_{sat} (\mumol m^{-2} s^{-1})','Interpreter', 'tex','FontSize',fontsize)
ylabel('\Delta\Psi (MPa)','Interpreter', 'tex','FontSize',fontsize)
box('on')
xl=[0 20];
yl=[0 1];
axis([xl yl])
text(xl(1)+0.05*(xl(2)-xl(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfc\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)

%Fig 1d (TLP vs A)
ax_1d=subplot(2,2,4)
TLP_A=fitlm([data_cal.A],[data_cal.TLP])
h=plot(TLP_A)
set(h,'Marker','none')
set(h,'LineWidth',linewidth)
set(h,'Color','black')
hold on
plot(A,TLP_calc*ones(size(A)),'--k','LineWidth',linewidth)
hold on
plot([data_cal.A],[data_cal.TLP],'ko','MarkerSize',markersize)
xlabel('\itA\rm_{sat} (\mumol m^{-2} s^{-1})','Interpreter', 'tex','FontSize',fontsize)
ylabel('\it\pi\rm_{tlp} (MPa)','Interpreter', 'tex','FontSize',fontsize)
box('on')
xl=[0 20];
yl=[0 2.8];
axis([xl yl])
text(xl(1)+0.05*(xl(2)-xl(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfd\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)
title([])
legend('off')



%%
%Figure 2 (A, gs, Kleaf and psi)
f2=figure('Units','centimeters','Position',[1 1 10.4 24]);
set(f2,'Units','normalized');

xl_dw=[0 0.03]; %max used to be 0.05
xl_psi_s=[-2 0]; %min used to be -3

%A
%Fig 2a (predicted A vs dw)
ax_2a=subplot(5,2,1)
plot(dw_range,A_dw,'--k','LineWidth',linewidth)
hold on
plot(dw_range,A_CF,':k','LineWidth',linewidth)
%xlabel('\Delta\itw\rm','Interpreter', 'tex','FontSize',fontsize)
ylabel('\itA\rm_{sat} (\mumol m^{-2} s^{-1})','Interpreter', 'tex','FontSize',fontsize)
box('on')
yl=[0 25];
axis([xl_dw yl])
xticks([0 0.015 0.03])
text(xl_dw(1)+0.05*(xl_dw(2)-xl_dw(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfa\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)
title([])
legend('off')

%Fig 2b (predicted A vs psi_s)
ax_2b=subplot(5,2,2)
plot(psi_s_range,A_psi,'--k','LineWidth',linewidth)
hold on
%xlabel('\Psi_s (MPa)','Interpreter', 'tex','FontSize',fontsize)
%ylabel('\itA\rm (\mumol m^{-2} s^{-1}','Interpreter', 'tex','FontSize',fontsize)
box('on')
yl=[0 25];
axis([xl_psi_s yl])
text(xl_psi_s(1)+0.05*(xl_psi_s(2)-xl_psi_s(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfb\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)
title([])
legend('off')


%gs
%Fig 2c (predicted gs vs dw)
ax_2c=subplot(5,2,3)
plot(dw_range,gs_dw,'--k','LineWidth',linewidth)
hold on
plot(dw_range,gs_CF,':k','LineWidth',linewidth)
%xlabel('\Delta\itw\rm','Interpreter', 'tex','FontSize',fontsize)
ylabel('\itg\rm_{sat} (mol m^{-2} s^{-1})','Interpreter', 'tex','FontSize',fontsize)
box('on')
yl=[0 1.3];
axis([xl_dw yl])
xticks([0 0.015 0.03])
text(xl_dw(1)+0.05*(xl_dw(2)-xl_dw(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfc\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)
title([])
legend('off')

%Fig 2d (predicted gs vs psi_s)
ax_2d=subplot(5,2,4)
plot(psi_s_range,gs_psi,'--k','LineWidth',linewidth)
hold on
%xlabel('\Psi_s (MPa)','Interpreter', 'tex','FontSize',fontsize)
%ylabel('\itA\rm (\mumol m^{-2} s^{-1}','Interpreter', 'tex','FontSize',fontsize)
box('on')
yl=[0 1.3];
axis([xl_psi_s yl])
text(xl_psi_s(1)+0.05*(xl_psi_s(2)-xl_psi_s(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfd\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)
title([])
legend('off')


%Kleaf
%Fig 2e (predicted Kleaf vs dw)
ax_2e=subplot(5,2,5)
plot(dw_range,Kleaf_dw,'--k','LineWidth',linewidth)
hold on
%xlabel('\Delta\itw\rm','Interpreter', 'tex','FontSize',fontsize)
ylabel('\itK\rm_{leaf} (mmol m^{-2} s^{-1} MPa^{-1})','Interpreter', 'tex','FontSize',fontsize)
box('on')
yl=[0 25];
axis([xl_dw yl])
xticks([0 0.015 0.03])
text(xl_dw(1)+0.05*(xl_dw(2)-xl_dw(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfe\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)
title([])
legend('off')

%Fig 2f (predicted Kleaf vs psi_s)
ax_2f=subplot(5,2,6)
plot(psi_s_range,Kleaf_psi,'--k','LineWidth',linewidth)
hold on
%xlabel('\Psi_s (MPa)','Interpreter', 'tex','FontSize',fontsize)
%ylabel('\itK\rm_{leaf} (mmol m^{-2} s^{-1} MPa^{-1})','Interpreter', 'tex','FontSize',fontsize)
box('on')
yl=[0 25];
axis([xl_psi_s yl])
text(xl_psi_s(1)+0.05*(xl_psi_s(2)-xl_psi_s(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bff\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)
title([])
legend('off')


%TLP
%Fig 2g (predicted TLP vs dw)
ax_2g=subplot(5,2,7)
plot(dw_range,TLP_dw,'--k','LineWidth',linewidth)
hold on
%xlabel('\Delta\itw\rm','Interpreter', 'tex','FontSize',fontsize)
ylabel('\it\pi\rm_{tlp} (MPa)','Interpreter', 'tex','FontSize',fontsize)
box('on')
yl=[0 6];
axis([xl_dw yl])
xticks([0 0.015 0.03])
text(xl_dw(1)+0.05*(xl_dw(2)-xl_dw(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfg\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)
title([])
legend('off')

%Fig 2h (predicted TLP vs psi_s)
ax_2h=subplot(5,2,8)
plot(psi_s_range,TLP_psi,'--k','LineWidth',linewidth)
hold on
%plot(psi_s_range,TLP_calc-psi_s_range,':k','LineWidth',linewidth) remove
%this for now
%xlabel('\Psi_s (MPa)','Interpreter', 'tex','FontSize',fontsize)
%ylabel('\it\pi\rm (MPa)','Interpreter', 'tex','FontSize',fontsize)
box('on')
yl=[0 6];
axis([xl_psi_s yl])
text(xl_psi_s(1)+0.05*(xl_psi_s(2)-xl_psi_s(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfh\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)
title([])
legend('off')


%psi
%Fig 2i (predicted psi vs dw)
ax_2i=subplot(5,2,9)
plot(dw_range,D_psi_dw,'--k','LineWidth',linewidth)
hold on
xlabel('\Delta\itw\rm','Interpreter', 'tex','FontSize',fontsize)
ylabel('\Delta\Psi (MPa)','Interpreter', 'tex','FontSize',fontsize)
box('on')
yl=[0 0.5]; %need to adjust this
axis([xl_dw yl])
xticks([0 0.015 0.03])
text(xl_dw(1)+0.05*(xl_dw(2)-xl_dw(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfi\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)
title([])
legend('off')

%Fig 2j (predicted psi vs psi_s)
ax_2j=subplot(5,2,10)
plot(psi_s_range,D_psi_psi,'--k','LineWidth',linewidth)
hold on
%plot(psi_s_range,psi_s_range,':k','LineWidth',linewidth)
xlabel('\Psi_s (MPa)','Interpreter', 'tex','FontSize',fontsize)
%ylabel('\Psi (MPa)','Interpreter', 'tex','FontSize',fontsize)
box('on')
yl=[0 0.5];
axis([xl_psi_s yl])
text(xl_psi_s(1)+0.05*(xl_psi_s(2)-xl_psi_s(1)),yl(2)-0.1*(yl(2)-yl(1)),'\bfj\rm','FontSize',fontsize,'Interpreter', 'tex')
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)
title([])
legend('off')


%%
%Figure 3 (relative costs of traits)
l=size(data_cal,1);
for i=1:l
    rel_costs(i).p=data_cal(i).p_chi;
    rel_costs(i).name='chi';
end

for i=1:l
    rel_costs(i+l).p=data_cal(i).p_K;
    rel_costs(i+l).name='K';
end
for i=1:l
    rel_costs(i+2*l).p=data_cal(i).p_pi;
    rel_costs(i+2*l).name='pi';
end

p_chi_mean=mean([data_cal.p_chi]);
p_chi_se=std([data_cal.p_chi])/sqrt(l);
p_K_mean=mean([data_cal.p_K]);
p_K_se=std([data_cal.p_K])/sqrt(l);
p_pi_mean=mean([data_cal.p_pi]);
p_pi_se=std([data_cal.p_pi])/sqrt(l);

%[p_costs,t,stats_costs]=anova1([rel_costs.p],{rel_costs.name})
%[c_costs,m_costs]=multcompare(stats_costs)

[p_costs,t,stats_costs]=kruskalwallis([rel_costs.p],{rel_costs.name})
[c_costs,m_costs]=multcompare(stats_costs)

f3=figure('Units','centimeters','Position',[1 1 8.9 7]);
set(f3,'Units','normalized');
boxplot([rel_costs.p],{rel_costs.name},'Labels',{'stomata','hydraulics', 'osmotic pressure'})
hold on
errorbar([1 2 3],[p_chi_mean p_K_mean p_pi_mean],[p_chi_se p_K_se p_pi_se],'ko','MarkerSize',markersize)
ylabel('Proportion of total cost (%)','Interpreter', 'tex','FontSize',fontsize)
text([1 2 3],[max([data_cal.p_chi])+5 max([data_cal.p_K])+5 max([data_cal.p_pi])+5],{'a','b','a'},'FontSize',10,'HorizontalAlignment','center')
set(gca,'FontSize',10)
set(gca,'LineWidth',linewidth)
xl=get(gca,'XLim');
yl=[0 60];
axis([xl yl])

%%
%to plot A vs g relation with constant k_fit

%to fit gs vs A data with constant k
f=funct_k(Ca,g_star,[data_test.gs],[data_test.A]);
k_fit=fminsearch(f,0.15)

gs_k_fit=0:0.01:1;
A_k_fit=k_fit*gs_k_fit*(Ca-g_star)./(1.6*k_fit+gs_k_fit);

A_k_fit_data=k_fit*[data_test.gs]*(Ca-g_star)./(1.6*k_fit+[data_test.gs]);

%to fit gs vs A data with constant Vcmax
%f=funct_V(Ca,g_star,509.5,[data_test.gs],[data_test.A]);
%V_fit=fminsearch(f,200)

%A_V_fit=k_fit*gs_k_fit*(Ca-g_star)./(1.6*k_fit+gs_k_fit);

%A_V_fit_data=k_fit*[data_test.gs]*(Ca-g_star)./(1.6*k_fit+[data_test.gs]);

%for constant k but gs varying (Reviewer 2) using mean of each species
%gs_k=(1.6*mean_k*A)./(mean_k*(Ca-g_star)-A);

%for constant k but gs varying (Reviewer 2)
%gs_k_fit=(1.6*k_fit*A)./(k_fit*(Ca-g_star)-A);

%k_test=0.055;
%gs_k_test=(1.6*k_test*A)./(k_test*(Ca-g_star)-A);

SSE_A_k_fit=sum((A_k_fit_data-[data_test.A]).^2);
SST_A_k_fit=sum(([data_test.A]-mean([data_test.A])).^2);
Rsquared_A_k_fit=1-SSE_A_k_fit/SST_A_k_fit


figure(1)
plot(gs_k_fit,A_k_fit,'--k','LineWidth',linewidth)
hold on
plot([data_test.gs],[data_test.A],'ro','MarkerSize',markersize)
xlabel('\itg\rm_{sat} (mol m^{-2} s^{-1})','Interpreter', 'tex','FontSize',fontsize)
ylabel('\itA\rm_{sat} (\mumol m^{-2} s^{-1})','Interpreter', 'tex','FontSize',fontsize)
box('on')
xl=[0 0.71];
yl=[0 20];
axis([xl yl])
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',linewidth)

%%



%%
%Figure S1 (mus)
l=size(data_cal,1);
for i=1:l
    mus(i).mu=data_cal(i).mu_chi;
    mus(i).name='chi';
end

for i=1:l
    mus(i+l).mu=data_cal(i).mu_K;
    mus(i+l).name='K';
end
for i=1:l
    mus(i+2*l).mu=data_cal(i).mu_piv;
    mus(i+2*l).name='piv';
end


mu_chi_se=std([data_cal.mu_chi])/sqrt(l);
mu_K_se=std([data_cal.mu_K])/sqrt(l);
mu_piv_se=std([data_cal.mu_piv])/sqrt(l);

%[p_costs,t,stats_costs]=anova1([rel_costs.p],{rel_costs.name})
%[c_costs,m_costs]=multcompare(stats_costs)

[p_mus,t,stats_mus]=kruskalwallis([mus.mu],{mus.name})
[c_mus,m_mus]=multcompare(stats_mus)

fS1=figure('Units','centimeters','Position',[1 1 8.9 7]);
set(fS1,'Units','normalized');
boxplot([mus.mu],{mus.name},'Labels',{'\it\mu_{\chi}\rm','\it\mu_{K}\rm','\it\mu_{\pi}\rm/\itv\rm'})
hold on
errorbar([1 2 3],[mu_chi_mean mu_K_mean mu_piv_mean],[mu_chi_se mu_K_se mu_piv_se],'ko','MarkerSize',markersize)
ylabel('Proportion of total cost (%)','Interpreter', 'tex','FontSize',fontsize)
text([1 2 3],[max([data_cal.p_chi])+5 max([data_cal.p_K])+5 max([data_cal.p_pi])+5],{'a','b','a'},'FontSize',10,'HorizontalAlignment','center')
set(gca,'FontSize',10)
set(gca,'LineWidth',linewidth)
xl=get(gca,'XLim');
yl=[0 60];
axis([xl yl])

%%
figure(1)
plot(gs_dw,A_dw)
axis([0 0.3 0 20])

figure(2)
plot(gs_psi,A_psi)
axis([0 0.3 0 20])

figure(3)
plot(dw_range,gs_dw)
axis([0 0.05 0 0.3])
%%
%to compare relative cost of Kleaf between single and multiple venation
l=size(data_cal,1);
p_K_s=[];
p_K_m=[];
for i=1:l
    if strcmp(data_cal(i).venation,'s')
        p_K_s=[p_K_s data_cal(i).p_K];
    elseif strcmp(data_cal(i).venation,'m')
        p_K_m=[p_K_m data_cal(i).p_K];
    end
end
        
p_K_s_mean=mean(p_K_s);
p_K_s_se=std(p_K_s)/sqrt(length(p_K_s));
p_K_m_mean=mean(p_K_m);
p_K_m_se=std(p_K_m)/sqrt(length(p_K_m));

[h_venation,p_venation]=ttest2(p_K_s,p_K_m,'Vartype','unequal','Tail','right')

boxplot([data_cal.p_K],{data_cal.venation},'Labels',{'single','multiple'})
hold on
errorbar([1 2],[p_K_s_mean p_K_m_mean],[p_K_s_se p_K_m_se],'ko')
xlabel('Venation type')
ylabel('Percentage of total cost')

%%
%{
%to calculate and plot mean A vs Kleaf, A vs gs and A vs TLP curves
Kleaf=0:1:20;
gs=0:0.01:0.5;

%for Kleaf
A_Kleaf=(mu_K_mean*(Ca-g_star)/(1.6*dw))^0.5*(TLP_calc+psi_s)^0.5*1e-3*Kleaf;
A_gs=(mu_chi_mean*(Ca-g_star)/1.6)^0.5*(1+(mu_K_mean*dw/mu_chi_mean)^0.5)*gs/(TLP_calc+psi_s)^0.5;

%%
figure(1)
plot(Kleaf,A_Kleaf,'--r')
hold on
for i=1:size(data_full,1)
    switch data_full(i).lineage
        case 'l'
            col='co'
        case 'f'
            col='bo'
        case 'g'
            col='go'
        case 'a'
            col='ro'
    end
    plot([data_full(i).Kleaf],[data_full(i).A],col)
end
for i=1:size(data_rest,1)
    switch data_rest(i).lineage
        case 'l'
            col='cx'
        case 'f'
            col='bx'
        case 'g'
            col='gx'
        case 'a'
            col='rx'
    end
    plot([data_rest(i).Kleaf],[data_rest(i).A],col)
end
xlabel('\itK\rm_{leaf} (mmol m^{-2} s^{-1} MPa^{-1})')
ylabel('\itA\rm (\mumol m^{-2} s^{-1})')

%%
figure(2)
plot(gs,A_gs,'--r')
hold on
for i=1:size(data_full,1)
    switch data_full(i).lineage
        case 'l'
            col='co'
        case 'f'
            col='bo'
        case 'g'
            col='go'
        case 'a'
            col='ro'
    end
    plot([data_full(i).gs],[data_full(i).A],col)
end
for i=1:size(data_rest,1)
    switch data_rest(i).lineage
        case 'l'
            col='cx'
        case 'f'
            col='bx'
        case 'g'
            col='gx'
        case 'a'
            col='rx'
    end
    plot([data_rest(i).gs],[data_rest(i).A],col)
end
xlabel('\itg\rm_s (mol m^{-2} s^{-1})')
ylabel('\itA\rm (\mumol m^{-2} s^{-1})')

%%
figure(3)
plot([TLP_calc TLP_calc],[0 20],'--r')
hold on
for i=1:size(data_full,1)
    switch data_full(i).lineage
        case 'l'
            col='co'
        case 'f'
            col='bo'
        case 'g'
            col='go'
        case 'a'
            col='ro'
    end
    plot([data_full(i).TLP],[data_full(i).A],col)
end

xlabel('TLP (MPa)')
ylabel('A (\mumol m^{-2} s^{-1})')
axis([0 2.8 0 inf])

%%
%look for sig differences between lineages in relative costs
[p_r_chi,t,stats_r_chi]=anova1([data_full.r_chi],{data_full.lineage})
multcompare(stats_r_chi)
figure(4)
boxplot([data_full.r_chi],{data_full.lineage})
xlabel('Lineage')
ylabel('Relative cost of stomata')



%%
[p_r_Kleaf,t,stats_r_Kleaf]=anova1([data_full.r_K],{data_full.lineage})
multcompare(stats_r_Kleaf)
figure(5)
boxplot([data_full.r_K],{data_full.lineage})
xlabel('Lineage')
ylabel('Relative cost of Kleaf')

%%
[p_r_pi,t,stats_r_pi]=anova1([data_full.r_pi],{data_full.lineage})
multcompare(stats_r_pi)
figure(6)
boxplot([data_full.r_pi],{data_full.lineage})
xlabel('Lineage')
ylabel('Relative cost of TLP')

%%
%look for differences in Kleaf costs across venations
[p_r_Kleaf_v,t,stats_r_Kleaf_v]=anova1([data_full.r_K],{data_full.venation})
multcompare(stats_r_Kleaf_v)
figure(7)
boxplot([data_full.r_K],{data_full.venation})
xlabel('Venation')
ylabel('Relative cost of Kleaf')

%%
%look for differences in total costs across lineages
p_r_total=anova1([data_full.r_total],{data_full.lineage})
figure(8)
boxplot([data_full.r_total],{data_full.lineage})

%%
plot([data_full.r_total],[data_full.A],'ko')

%%
%boxplots of relative costs
figure(1)
boxplot([data_full.r_K],{data_full.lineage})

figure(2)
plot([data_full.A],[data_full.cost_total],'ko')
hold on
plot([0 20],[0 20],'--r')
%}
%{
boxplot([data.mu_chi],{data.lineage})

boxplot([data.mu_piv],{data.lineage})
%}

%plot(data.gs,data.A, 'ko')

