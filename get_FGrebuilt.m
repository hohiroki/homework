function [F_rebuilt,G_rebuilt]=get_FGrebuilt(rou,u,v,H,p,rourou,rouu,rouv,rouH,roua,iL,jL,iR,jR)
FL=[rou(iL,jL)*u(iL,jL),rou(iL,jL)*u(iL,jL)^2+p(iL,jL),rou(iL,jL)*u(iL,jL)*v(iL,jL),rou(iL,jL)*u(iL,jL)*H(iL,jL)]';
FR=[rou(iR,jR)*u(iR,jR),rou(iR,jR)*u(iR,jR)^2+p(iR,jR),rou(iR,jR)*u(iR,jR)*v(iR,jR),rou(iR,jR)*u(iR,jR)*H(iR,jR)]';
GL = [rou(iL,jL)*v(iL,jL),rou(iL,jL)*u(iL,jL)*v(iL,jL),rou(iL,jL)*v(iL,jL)^2+p(iL,jL),rou(iL,jL)*v(iL,jL)*H(iL,jL)]';
GR = [rou(iR,jR)*v(iR,jR),rou(iR,jR)*u(iR,jR)*v(iR,jR),rou(iR,jR)*v(iR,jR)^2+p(iR,jR),rou(iR,jR)*v(iR,jR)*H(iR,jR)]';
F_langta = abs([rouu-roua,rouu,rouu,rouu+roua]);
G_langta = abs([rouv-roua,rouv,rouv,rouv+roua]);
F_afa = [(p(iR,jR)-p(iL,jL)-rourou*roua*(u(iR,jR)-u(iL,jL)))/2/roua^2,rou(iR,jR)-rou(iL,jL)-(p(iR,jR)-p(iL,jL))/roua^2,...
rourou*(v(iR,jR)-v(iL,jL)),(p(iR,jR)-p(iL,jL)+rourou*roua*(u(iR,jR)-u(iL,jL)))/2/roua^2];
G_afa = [(p(iR,jR)-p(iL,jL)-rourou*roua*(v(iR,jR)-v(iL,jL)))/2/roua^2,rou(iR,jR)-rou(iL,jL)-(p(iR,jR)-p(iL,jL))/roua^2,...
-rourou*(u(iR,jR)-u(iL,jL)),(p(iR,jR)-p(iL,jL)+rourou*roua*(v(iR,jR)-v(iL,jL)))/2/roua^2];
F_r=[1,rouu-roua,rouv,rouH-rouu*roua;1,rouu,rouv,0.5*(rouu^2+rouv^2);0,0,1,rouv;1,rouu+roua,rouv,rouH+rouu*roua]';
G_r=[1,rouu,rouv-roua,rouH-rouv*roua;1,rouu,rouv,0.5*(rouu^2+rouv^2);0,-1,0,-rouu;1,rouu,rouv+roua,rouH+rouv*roua]';
F_rebuilt = 0.5*(FL+FR)-0.5* get_arfaLangtaR(F_afa,F_langta,F_r);            
G_rebuilt = 0.5*(GL+GR)-0.5* get_arfaLangtaR(G_afa,G_langta,G_r);
end