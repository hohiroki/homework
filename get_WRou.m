function [rourou,rouu,rouv,rouH,roua,roup]=get_WRou(rou_n,u_n,v_n,H_n,iL,jL,iR,jR)
gama = 1.26;
rouL = sqrt(rou_n(iL,jL));
rouR = sqrt(rou_n(iR,jR)); 
rourou=rouL*rouR;
rouu=(u_n(iL,jL)*rouL+u_n(iR,jR)*rouR)/(rouL+rouR);
rouv=(v_n(iL,jL)*rouL+v_n(iR,jR)*rouR)/(rouL+rouR);
rouH=(H_n(iL,jL)*rouL+H_n(iR,jR)*rouR)/(rouL+rouR);
roua = ((gama-1)*(rouH-0.5*(rouu^2+rouv^2)))^0.5;
roup=(rouH-0.5*(rouu^2+rouv^2))*rourou*(gama-1)/gama;
end