%%
close all
clc
clear all
A=ones(1,200);
MatrixSocial=diag(A);
Sum_GSA=zeros(1,1);
Sum_MEC_GSA=zeros(1,1);
Sum_V2V_GSA=zeros(1,1);
Sum_exh1=zeros(1,1);
Sum_Rand1=zeros(1,1);
Sum_Rand2=zeros(1,1);
Sum_GSA_DistanceVP=zeros(1,1);
Sum_HA=zeros(1,1);
loop_N_c=0;
Sum_GSA_N_c_Varying=zeros(1,1);
Sum_MEC_GSA_N_c_Varying=zeros(1,1);
Sum_V2V_GSA_N_c_Varying=zeros(1,1);
Sum_Rand1_N_c_Varying=zeros(1,1);
Sum_Rand2_N_c_Varying=zeros(1,1);
Sum_HA_N_c_Varying=zeros(1,1);
Sum_exh1_N_c_Varying=zeros(1,1);
Sum_GSA_DistanceVP_N_c_Varying=zeros(1,1);

Matrix_N_C=20:5:50;
CDF_Sum_GSA=[];
CDF_Sum_Rand2=[];
CDF_Sum_HA=[];
CDF_Sum_GSA_DistanceVP=[];
k=1e-27;
f_loc=1e5;
f_mec=1e7;
f_vec=1e6;
D_total=1e6;
C=1e4;  %1bit CPU cycle
Ttol=1;
for N_C=Matrix_N_C
    overlay_N_C=20-N_C;
     if overlay_N_C<=0
        overlay_N_C=1;
    else
    end
    loop_N_c=loop_N_c+1;
    
    for loop_n=1:50
        disp(['N_C=',num2str(N_C),',loop_n=',num2str(loop_n)])
        loop_n;
        Radius_eNodeB=250;                           
        Length_manhattan_mobility_model=500;         
        Length_vehicle=5;                            
        Width_manhattan_mobility_model=500;          
        Width_lane=4;                              
        Speed_vehicle=10;                           
        Distance_eNodeB_to_lane=10;                  
        Number_single_direction_lane=3;             
        N_M=N_C;                                      
        Number_RB=10;                                
        Carrier_frequence=2e9;                       
        Bandwidth=10e6;                              
        P_C_max=26;                                 
        P_M_max=40;                                
        RB=1e5;   
        P_noise=-104;                               
        Distance_V2V_links_threshold=130;           
        Shadowing_standard_deviation_Cellular=8;    
        Shadowing_standard_deviation_V2V=3;          
        SINR_V2V_th=5;                               
        beta_d2d=0.5.*ones(1,N_M);
        beta_mec=1-beta_d2d;
        B0=180000;
        Xi0=1e6;
        Pos_MBS=[250,250];
        overlay_N_C;
        Pos_overlay_c=[500.*rand(overlay_N_C,1),260.*ones(overlay_N_C,1)];
        Dis_overlay_c=(  ( Pos_overlay_c(:,1)-Pos_MBS(1,1) ).^2+( Pos_overlay_c(:,2)-Pos_MBS(1,2) ).^2  ).^(1/2);
        PL_overlay_c=128.1+37.6.*log10((Dis_overlay_c./1000));
        Pr_overlay_c=P_C_max-PL_overlay_c;
        SINR_overlay_c=Pr_overlay_c+104;
        Sum_overlay_c=sum(B0.*log2(1+ 10.^(SINR_overlay_c./10)));
     
        Pos_v=[500.*rand(N_M,1),250.*ones(N_M,1)];
        Pos_vv=[500.*rand(N_M*3,1),240.*ones(N_M*3,1)];
        Dis_v2vv=( ( Pos_v(:,1)-Pos_vv(:,1)' ).^2+( 10  ).^2 ).^(1/2);
        Dis_v2i=((mod(Pos_v(:,1),100)).^2+(10).^2).^(1/2);
        Pos_c=[500.*rand(N_C,1),260.*ones(N_C,1)];
        Dis_c=(  ( Pos_c(:,1)-Pos_MBS(1,1) ).^2+( Pos_c(:,2)-Pos_MBS(1,2) ).^2  ).^(1/2);
        % figure
        % figure()
        % plot(Pos_MBS(1,1),Pos_MBS(1,2),'b^')
        % hold on
        % plot(Pos_v(:,1),Pos_v(:,2),'r*')
        % hold on
        % plot(Pos_vv(:,1),Pos_vv(:,2),'r.')
        % hold on
        % plot(Pos_c(:,1),Pos_c(:,2),'k+')
        

        PL_v2i=128.1+37.6.*log10((Dis_v2i./1000));
        Pr_v2i=P_M_max-PL_v2i;
        Matrix_R=B0.*log2(1+10.^((Pr_v2i(1:20,:)+104)./10));
        R=sum(Matrix_R);
       EE_Overlay_LOOP(1,loop_n)=R./(k.*f_mec.^3+20.*1e-3.*10^(P_M_max./10));
        
       [EE_JCRAA_LOOP(1,loop_n)]=JCRAA(N_M,Dis_c,P_C_max,Dis_v2vv,P_M_max,Dis_v2i,Pos_v,Pos_c,N_C,Pos_vv,MatrixSocial,Ttol,f_mec,f_loc,f_vec,D_total,C,B0,beta_mec,beta_d2d,Sum_overlay_c,k);
       [EE_MDSS_LOOP(1,loop_n)]=MDSS(N_M,Dis_c,P_C_max,Dis_v2vv,P_M_max,Dis_v2i,Pos_v,Pos_c,N_C,Pos_vv,MatrixSocial,Ttol,f_mec,f_loc,f_vec,D_total,C,B0,beta_mec,beta_d2d,Sum_overlay_c,k);
       [EE_RSU_LOOP(1,loop_n)]=RSU(N_M,Dis_c,P_C_max,Dis_v2vv,P_M_max,Dis_v2i,Pos_v,Pos_c,N_C,Pos_vv,MatrixSocial,Ttol,f_mec,f_loc,f_vec,D_total,C,B0,beta_mec,beta_d2d,Sum_overlay_c,k);


        EE=zeros(1,50);
        for ProLoop=2:50
            W_mec=0.8.*ones(1,N_M);
            W_d2d=0.2.*ones(1,N_M);
            PL_c=128.1+37.6.*log10((Dis_c./1000));
            Pr_c=P_C_max-PL_c;
            
            PL_v2vv=128.1+37.6.*log10((Dis_v2vv./1000));
            Pr_v2vv=P_M_max.*W_d2d'-PL_v2vv;
            
            PL_v2i=128.1+37.6.*log10((Dis_v2i./1000));
            Pr_v2i=P_M_max.*W_d2d'-PL_v2i;
            
            
            Dis_v_to_c=( ( Pos_v(:,1)-Pos_c(:,1)' ).^2+( Pos_v(:,2)-Pos_c(:,2)' ).^2 ).^(1/2);   % (i,j)值表示第i个车辆用户到第j个蜂窝用户的距离
            PL_v_to_c=128.1+37.6.*log10((Dis_v_to_c./1000));
            I_v_to_c=P_M_max-PL_v_to_c;
            
       
            Dis_v2v_to_v2i=Dis_v2i.*ones(N_M,N_C);
            PL_v2v_to_v2i=128.1+37.6.*log10((Dis_v2v_to_v2i./1000));
            I_v2v_to_v2i=P_M_max.*W_d2d'-PL_v2v_to_v2i;
            
            Dis_c_to_v2i=((mod(Pos_c(:,1),100)).^2+(10).^2).^(1/2).*ones(N_C,N_M);
            PL_c_to_v2i=128.1+37.6.*log10((Dis_c_to_v2i./1000));
            I_c_to_v2i=P_C_max-PL_c_to_v2i;
            
            
            Dis_c_to_v2v=( ( Pos_c(:,1)-Pos_vv(:,1)' ).^2+( Pos_c(:,2)-Pos_vv(:,2)' ).^2 ).^(1/2);
            PL_c_to_v2v=128.1+37.6.*log10((Dis_c_to_v2v./1000));
            I_c_to_v2v=P_C_max-PL_c_to_v2v;
            
          
            SINR_c=Pr_c-10.*log10(10.^(I_v_to_c./10)+10.^(-104./10));                              % (i,j)值表示第i个车辆用户与第j个蜂窝用户共用RB时，第j个蜂窝用户的SINR值（dB）
            SINR_v2i=Pr_v2i-10.*log10(10.^(I_v2v_to_v2i./10)+10.^(I_c_to_v2i./10)+10.^(-104./10)); %（i,j）值表示第i个车辆用户与第j个蜂窝用户共用RB时，第j个v2i用户的SINR值（dB）
            SINR_v2v=Pr_v2vv-10.*log10(10.^(I_c_to_v2v./10)+10.^(-104./10));                       %
            social_v2v=MatrixSocial(1:N_M,1:N_M*3);
            SINR_v2v=SINR_v2v.*social_v2v;
            
            

           
            costMat=-SINR_v2v;
            [VP_KM,cost] = munkres(costMat);
            SINR_v2v_fvp=zeros(1,N_M);
            for i=1:N_C
                SINR_v2v_fvp(i,:)=SINR_v2v(i,VP_KM(i));
            end
            % 
            for i=1:N_M
                R_us=max(log2(1+10.^(SINR_v2i(i,:))));
                R_um=max(log2(1+10.^(SINR_v2v(i,:))));
            beta_mec_max=(Ttol.*R_us.*f_mec)./(f_mec.*D_total+R_us.*D_total.*C);
            beta_d2d_max=(Ttol.*R_um.*f_vec)./(f_vec.*D_total+R_um.*D_total.*C);
            beta_mec_d2d_min=1-(Ttol.*f_loc)./(D_total.*C);
            for x1=1:-0.01:0
                for x2=1:-0.01:0
                   if ( x1<beta_mec_max) &&  (x2<beta_d2d_max)    &&   (x1+x2>beta_mec_d2d_min)
                      beta_mec(i)=x1;
                      beta_d2d(i)=x2;
                       break
                   else
                   end
                end
            end
            end
            
            
            %-----                                   --------------------------
            SINR_V=SINR_v2v_fvp+SINR_v2i;    
            SINR_c;                         
            EE_V=B0.*log2(1+10.^(SINR_V./10));
            Utility_v=zeros(N_C,N_M);    
            for i=1:N_C
                [B,I] = sort(SINR_V(i,:),'descend');
                [~,II] = sort(I);
                Utility_v(i,:)=I;
            end
            Utility_c=zeros(N_M,N_C);    
            for i=1:N_M
                [B,I] = sort(SINR_c(i,:),'descend');
                [~,II] = sort(I);
                Utility_c(i,:)=I;
            end
            stablematch = galeshapley(N_C, Utility_v, Utility_c);
            Match=[(1:N_C)',stablematch];     %
            [W_mec,W_d2d]=PC(Match,W_mec,W_d2d,SINR_v2i,SINR_v2v_fvp,D_total,beta_mec,beta_d2d,P_M_max,PL_c_to_v2i,B0,PL_c_to_v2v);
     
            [SINR_c,SINR_v2i,SINR_v2v,SINR_v2v_Social]=SINR_social(social_v2v,P_C_max,PL_c,P_M_max,W_mec,W_d2d,...
                PL_v2vv,PL_v2i,PL_v_to_c,PL_v2v_to_v2i,PL_c_to_v2i,PL_c_to_v2v,VP_KM);
            Sum=Sum_overlay_c;
            for i=1:N_C
                Sum=Sum+B0.*log2(1+ 10.^(SINR_c(Match(i,1),Match(i,2))./10))+B0.*log2(1+ 10.^(SINR_v2i(Match(i,1),Match(i,2))./10))+B0.*log2(1+ 10.^(SINR_v2v_Social(Match(i,2),Match(i,1))./10));
            end
            P=k.*f_loc.^3+k.*f_mec.^3+k.*f_vec.^3+1e-3.*N_C.*10^(P_C_max./10)+1e-3.*sum(P_M_max.*W_d2d+P_M_max.*W_mec);
            EE(ProLoop)=Sum./P;
            while (abs(EE(ProLoop)-EE(ProLoop-1))>0.1)
                break
            end
        end
        EE_Pro(1,loop_n)=EE(ProLoop);        
    end    % for loop
    
    


    EE_Pro_N_c_Varying(1,loop_N_c)=mean(EE_Pro);
    EE_ProRandRA_N_c_Varying(1,loop_N_c)=mean(EE_JCRAA_LOOP);
    EE_MDSS_N_c_Varying(1,loop_N_c)=mean(EE_MDSS_LOOP);
    EE_RSU_N_c_Varying(1,loop_N_c)=mean(EE_RSU_LOOP);
    EE_Overlay_N_c_Varying(1,loop_N_c)=mean(EE_Overlay_LOOP);
    
    Sum_MEC_GSA_N_c_Varying(1,loop_N_c)=mean(Sum_MEC_GSA);
    Sum_V2V_GSA_N_c_Varying(1,loop_N_c)=mean(Sum_V2V_GSA);
    Sum_Rand1_N_c_Varying(1,loop_N_c)=mean(Sum_Rand1);
    Sum_Rand_N_c_Varying(1,loop_N_c)=mean(Sum_Rand2);
    Sum_HA_N_c_Varying(1,loop_N_c)=mean(Sum_HA);
    Sum_exh1_N_c_Varying(1,loop_N_c)=mean(Sum_exh1);
    Sum_GSA_DistanceVP_N_c_Varying(1,loop_N_c)=mean(Sum_GSA_DistanceVP);


end


%%

Total_user=ceil((20-Matrix_N_C))+Matrix_N_C.*2;

EE_NOMA_Pro=EE_Pro_N_c_Varying./Total_user;
EE_NOMA_RA=EE_ProRandRA_N_c_Varying./Total_user;
EE_NOMA_MDSS=EE_MDSS_N_c_Varying./Total_user;
EE_NOMA_RSU=EE_RSU_N_c_Varying./Total_user;
EE_OVERLAY=EE_Overlay_N_c_Varying./Total_user;

EE_MEC_NOMA_GSA=Sum_MEC_GSA_N_c_Varying./Total_user;
EE_V2V_NOMA_GSA=Sum_V2V_GSA_N_c_Varying./Total_user;
EE_NOMA_HA=Sum_HA_N_c_Varying./Total_user;




%%

figure_FontSize=12;

%% EE----------------------------------------不同算法
figure()
plot(Matrix_N_C,EE_NOMA_Pro,'bo-','linewidth',2)
hold on
plot(Matrix_N_C,EE_NOMA_RA,'g^-','linewidth',2)
hold on
plot(Matrix_N_C,EE_NOMA_MDSS,'r>-','linewidth',2)
hold on
plot(Matrix_N_C,EE_NOMA_RSU,'ms-','linewidth',2)
hold on
plot(Matrix_N_C,EE_OVERLAY,'k<-','linewidth',2)
grid on
set(gca,'linewidth',1.5,'fontsize',12,'fontname','Times');
legend('NOMA-JCCRAA','NOMA-JCCRAA-RA','NOMA-MDSS','NOMA-RSU','NOMA-Overlay','FontSize',figure_FontSize)
xlabel('The number of VUEs','FontSize',figure_FontSize)
ylabel('Energy Efficiency (bits/J)','FontSize',figure_FontSize)
disp('The varying number of RB Gain')
mean((Sum_GSA_N_c_Varying-Sum_Rand2_N_c_Varying)./Sum_Rand2_N_c_Varying)
mean((Sum_GSA_N_c_Varying-Sum_HA_N_c_Varying)./Sum_HA_N_c_Varying)
xlim([20 50])




%save('M20_5_50_RSU.mat','EE_NOMA_RSU')






