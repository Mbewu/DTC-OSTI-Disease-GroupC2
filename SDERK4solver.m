% Just run >> SDERK4solver

%% Params
global param;
param.gammaV = 510;
param.gammaVA = 619.2;
param.gammaVH = 1.02;
param.alphaV = 1.7;
param.aV1 = 100;
param.aV2 = 23000;
param.bHD = 4;
param.aR = 1;
param.gammaHV = 0.34;
param.bHF = 0.01;
param.bIE = 0.066;
param.aI = 1.5;
param.bMD = 1;
param.bMV = 0.0037;
param.aM = 1;
param.bF = 250000;
param.cF = 2000;
param.bFH = 17;
param.aF = 8;
param.bEM = 8.3;
param.bEI = 2.72;
param.aE = 0.4;
param.bPM = 11.5;
param.aP = 0.4;
param.bA = 0.043;
param.gammaAV = 146.2;
param.aA = 0.043;
param.r = (3 * 10 ^ (-5));

tBegin=0;
tEnd=10;
dt=.001;
t=tBegin:dt:tEnd;
N=length(t);

%% Preallocate the predictor and corrector solution vectors 
Vc=zeros(N,1);
Hc=zeros(N,1);
Ic=zeros(N,1);
Mc=zeros(N,1);
Fc=zeros(N,1);
Rc=zeros(N,1);
Ec=zeros(N,1);
Pc=zeros(N,1);
Ac=zeros(N,1);
Sc=zeros(N,1);

Vp=zeros(N,1);
Hp=zeros(N,1);
Ip=zeros(N,1);
Mp=zeros(N,1);
Fp=zeros(N,1);
Rp=zeros(N,1);
Ep=zeros(N,1);
Pp=zeros(N,1);
Ap=zeros(N,1);
Sp=zeros(N,1);

Vk1=zeros(N,1);
Hk1=zeros(N,1);
Ik1=zeros(N,1);
Mk1=zeros(N,1);
Fk1=zeros(N,1);
Rk1=zeros(N,1);
Ek1=zeros(N,1);
Pk1=zeros(N,1);
Ak1=zeros(N,1);
Sk1=zeros(N,1);

Vk2=zeros(N,1);
Hk2=zeros(N,1);
Ik2=zeros(N,1);
Mk2=zeros(N,1);
Fk2=zeros(N,1);
Rk2=zeros(N,1);
Ek2=zeros(N,1);
Pk2=zeros(N,1);
Ak2=zeros(N,1);
Sk2=zeros(N,1);

Vk3=zeros(N,1);
Hk3=zeros(N,1);
Ik3=zeros(N,1);
Mk3=zeros(N,1);
Fk3=zeros(N,1);
Rk3=zeros(N,1);
Ek3=zeros(N,1);
Pk3=zeros(N,1);
Ak3=zeros(N,1);
Sk3=zeros(N,1);

Vk4=zeros(N,1);
Hk4=zeros(N,1);
Ik4=zeros(N,1);
Mk4=zeros(N,1);
Fk4=zeros(N,1);
Rk4=zeros(N,1);
Ek4=zeros(N,1);
Pk4=zeros(N,1);
Ak4=zeros(N,1);
Sk4=zeros(N,1);

D=zeros(N,1);

tempVec = zeros(tEnd,10);

%% Initial Conditions
y0 = [ 0.01 , 1 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 0.1 ];
Vc(1)=0.01;
Hc(1)=1;
Ic(1)=0;
Mc(1)=0;
Fc(1)=0;
Rc(1)=0;
Ec(1)=1;
Pc(1)=1;
Ac(1)=1;
Sc(1)=0.1;
D(1) = 1 - Hc(1) - Rc(1) - Ic(1);

tempVec(1,:) = [ 0.01 , 1 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 0.1 ];

%% Solve the ODES

for i=2:length(Vc)
    Vp(i) = Vc(i-1)+dt*(param.gammaV * Ic(i-1) - param.gammaVA * Sc(i-1) * Ac(i-1) * Vc(i-1) - param.gammaVH * Hc(i-1) * Vc(i-1) - param.alphaV * Vc(i-1) - (param.aV1 * Vc(i-1)) / (1 + (param.aV2 * Vc(i-1))));
    Hp(i) = Hc(i-1)+dt*(param.bHD * D(i-1) * (Hc(i-1) + Rc(i-1)) + param.aR * Rc(i-1) - param.gammaHV * Vc(i-1) * Hc(i-1) - param.bHF * Fc(i-1) * Hc(i-1));
    Ip(i) = Ic(i-1)+dt*(param.gammaHV * Vc(i-1) * Hc(i-1) - param.bIE * Ec(i-1) * Ic(i-1) - param.aI * Ic(i-1));
	Mp(i) = Mc(i-1)+dt*((param.bMD * D(i-1) + param.bMV * Vc(i-1)) * (1 - Mc(i-1)) - param.aM * Mc(i-1));
	Fp(i) = Fc(i-1)+dt*(param.bF * Mc(i-1) + param.cF * Ic(i-1) - param.bFH * Hc(i-1) * Fc(i-1) - param.aF * Fc(i-1));
	Rp(i) = Rc(i-1)+dt*(param.bHF * Fc(i-1) * Hc(i-1) - param.aR * Rc(i-1));
    Ep(i) = Ec(i-1)+dt*(param.bEM * Mc(i-1) * Ec(i-1) - param.bEI * Ec(i-1) * Ic(i-1) + param.aE * (1 - Ec(i-1)));
	Pp(i) = Pc(i-1)+dt*(param.bPM * Mc(i-1) * Pc(i-1) + param.aP * (1 - Pc(i-1)));
	Ap(i) = Ac(i-1)+dt*(param.bA * Pc(i-1) - param.gammaAV * Sc(i-1) * Ac(i-1) * Vc(i-1) - param.aA * Ac(i-1));
	Sp(i) = Sc(i-1)+dt*(param.r * Pc(i-1) * (1 - Sc(i-1)));
    
        
    
    Vk1(i) = dt*(param.gammaV * Ip(i) - param.gammaVA * Sp(i) * Ap(i) * Vp(i) - param.gammaVH * Hp(i) * Vp(i) - param.alphaV * Vp(i) - (param.aV1 * Vp(i)) / (1 + (param.aV2 * Vp(i))));
    Hk1(i) = dt*(param.bHD * D(i-1) * (Hp(i) + Rp(i)) + param.aR * Rp(i) - param.gammaHV * Vp(i) * Hp(i) - param.bHF * Fp(i) * Hp(i));
    Ik1(i) = dt*(param.gammaHV * Vp(i) * Hp(i) - param.bIE * Ep(i) * Ip(i) - param.aI * Ip(i));
	Mk1(i) = dt*((param.bMD * D(i-1) + param.bMV * Vp(i)) * (1 - Mp(i)) - param.aM * Mp(i));
	Fk1(i) = dt*(param.bF * Mp(i) + param.cF * Ip(i) - param.bFH * Hp(i) * Fp(i) - param.aF * Fp(i));
	Rk1(i) = dt*(param.bHF * Fp(i) * Hp(i) - param.aR * Rp(i));
    Ek1(i) = dt*(param.bEM * Mp(i) * Ep(i) - param.bEI * Ep(i) * Ip(i) + param.aE * (1 - Ep(i)));
	Pk1(i) = dt*(param.bPM * Mp(i) * Pp(i) + param.aP * (1 - Pp(i)));
	Ak1(i) = dt*(param.bA * Pp(i) - param.gammaAV * Sp(i) * Ap(i) * Vp(i) - param.aA * Ap(i));
	Sk1(i) = dt*(param.r * Pp(i) * (1 - Sp(i)));
    
    
    
    Vk2(i) = dt*(param.gammaV * (Ip(i)+(0.5*Ik1(i))) - param.gammaVA * (Sp(i)+(0.5*Sk1(i))) * (Ap(i)+(0.5*Ak1(i))) * (Vp(i)+(0.5*Vk1(i))) - param.gammaVH * (Hp(i)+(0.5*Hk1(i))) * (Vp(i)+(0.5*Vk1(i))) - param.alphaV * (Vp(i)+(0.5*Vk1(i))) - (param.aV1 * (Vp(i)+(0.5*Vk1(i)))) / (1 + (param.aV2 * (Vp(i)+(0.5*Vk1(i))))));
    Hk2(i) = dt*(param.bHD * D(i-1) * ((Hp(i)+(0.5*Hk1(i))) + (Rp(i)+(0.5*Rk1(i)))) + param.aR * (Rp(i)+(0.5*Rk1(i))) - param.gammaHV * (Vp(i)+(0.5*Vk1(i))) * (Hp(i)+(0.5*Hk1(i))) - param.bHF * (Fp(i)+(0.5*Fk1(i))) * (Hp(i)+(0.5*Hk1(i))));
    Ik2(i) = dt*(param.gammaHV * (Vp(i)+(0.5*Vk1(i))) * (Hp(i)+(0.5*Hk1(i))) - param.bIE * (Ep(i)+(0.5*Ek1(i))) * (Ip(i)+(0.5*Ik1(i))) - param.aI * (Ip(i)+(0.5*Ik1(i))));
	Mk2(i) = dt*((param.bMD * D(i-1) + param.bMV * (Vp(i)+(0.5*Vk1(i)))) * (1 - (Mp(i)+(0.5*Mk1(i)))) - param.aM * (Mp(i)+(0.5*Mk1(i))));
	Fk2(i) = dt*(param.bF * (Mp(i)+(0.5*Mk1(i))) + param.cF * (Ip(i)+(0.5*Ik1(i))) - param.bFH * (Hp(i)+(0.5*Hk1(i))) * (Fp(i)+(0.5*Fk1(i))) - param.aF * (Fp(i)+(0.5*Fk1(i))));
	Rk2(i) = dt*(param.bHF * (Fp(i)+(0.5*Fk1(i))) * (Hp(i)+(0.5*Hk1(i))) - param.aR * (Rp(i)+(0.5*Rk1(i))));
    Ek2(i) = dt*(param.bEM * (Mp(i)+(0.5*Mk1(i))) * (Ep(i)+(0.5*Ek1(i))) - param.bEI * (Ep(i)+(0.5*Ek1(i))) * (Ip(i)+(0.5*Ik1(i))) + param.aE * (1 - (Ep(i)+(0.5*Ek1(i)))));
	Pk2(i) = dt*(param.bPM * (Mp(i)+(0.5*Mk1(i))) * (Pp(i)+(0.5*Pk1(i))) + param.aP * (1 - (Pp(i)+(0.5*Pk1(i)))));
	Ak2(i) = dt*(param.bA * (Pp(i)+(0.5*Pk1(i))) - param.gammaAV * (Sp(i)+(0.5*Sk1(i))) * (Ap(i)+(0.5*Ak1(i))) * (Vp(i)+(0.5*Vk1(i))) - param.aA * (Ap(i)+(0.5*Ak1(i))));
	Sk2(i) = dt*(param.r * (Pp(i)+(0.5*Pk1(i))) * (1 - (Sp(i)+(0.5*Sk1(i)))));
    
    
 
    Vk3(i) = dt*(param.gammaV * (Ip(i)+(0.5*Ik2(i))) - param.gammaVA * (Sp(i)+(0.5*Sk2(i))) * (Ap(i)+(0.5*Ak2(i))) * (Vp(i)+(0.5*Vk2(i))) - param.gammaVH * (Hp(i)+(0.5*Hk2(i))) * (Vp(i)+(0.5*Vk2(i))) - param.alphaV * (Vp(i)+(0.5*Vk2(i))) - (param.aV1 * (Vp(i)+(0.5*Vk2(i)))) / (1 + (param.aV2 * (Vp(i)+(0.5*Vk2(i))))));
    Hk3(i) = dt*(param.bHD * D(i-1) * ((Hp(i)+(0.5*Hk2(i))) + (Rp(i)+(0.5*Rk2(i)))) + param.aR * (Rp(i)+(0.5*Rk2(i))) - param.gammaHV * (Vp(i)+(0.5*Vk2(i))) * (Hp(i)+(0.5*Hk2(i))) - param.bHF * (Fp(i)+(0.5*Fk2(i))) * (Hp(i)+(0.5*Hk2(i))));
    Ik3(i) = dt*(param.gammaHV * (Vp(i)+(0.5*Vk2(i))) * (Hp(i)+(0.5*Hk2(i))) - param.bIE * (Ep(i)+(0.5*Ek2(i))) * (Ip(i)+(0.5*Ik2(i))) - param.aI * (Ip(i)+(0.5*Ik2(i))));
	Mk3(i) = dt*((param.bMD * D(i-1) + param.bMV * (Vp(i)+(0.5*Vk2(i)))) * (1 - (Mp(i)+(0.5*Mk2(i)))) - param.aM * (Mp(i)+(0.5*Mk2(i))));
	Fk3(i) = dt*(param.bF * (Mp(i)+(0.5*Mk2(i))) + param.cF * (Ip(i)+(0.5*Ik2(i))) - param.bFH * (Hp(i)+(0.5*Hk2(i))) * (Fp(i)+(0.5*Fk2(i))) - param.aF * (Fp(i)+(0.5*Fk2(i))));
	Rk3(i) = dt*(param.bHF * (Fp(i)+(0.5*Fk2(i))) * (Hp(i)+(0.5*Hk2(i))) - param.aR * (Rp(i)+(0.5*Rk2(i))));
    Ek3(i) = dt*(param.bEM * (Mp(i)+(0.5*Mk2(i))) * (Ep(i)+(0.5*Ek2(i))) - param.bEI * (Ep(i)+(0.5*Ek2(i))) * (Ip(i)+(0.5*Ik2(i))) + param.aE * (1 - (Ep(i)+(0.5*Ek2(i)))));
	Pk3(i) = dt*(param.bPM * (Mp(i)+(0.5*Mk2(i))) * (Pp(i)+(0.5*Pk2(i))) + param.aP * (1 - (Pp(i)+(0.5*Pk2(i)))));
	Ak3(i) = dt*(param.bA * (Pp(i)+(0.5*Pk2(i))) - param.gammaAV * (Sp(i)+(0.5*Sk2(i))) * (Ap(i)+(0.5*Ak2(i))) * (Vp(i)+(0.5*Vk2(i))) - param.aA * (Ap(i)+(0.5*Ak2(i))));
	Sk3(i) = dt*(param.r * (Pp(i)+(0.5*Pk2(i))) * (1 - (Sp(i)+(0.5*Sk2(i)))));
    
    
    Vk4(i) = dt*(param.gammaV * (Ip(i)+(Ik3(i))) - param.gammaVA * (Sp(i)+(Sk3(i))) * (Ap(i)+(Ak3(i))) * (Vp(i)+(Vk3(i))) - param.gammaVH * (Hp(i)+(Hk3(i))) * (Vp(i)+(Vk3(i))) - param.alphaV * (Vp(i)+(Vk3(i))) - (param.aV1 * (Vp(i)+(Vk3(i)))) / (1 + (param.aV2 * (Vp(i)+(Vk3(i))))));
    Hk4(i) = dt*(param.bHD * D(i-1) * ((Hp(i)+(Hk3(i))) + (Rp(i)+(Rk3(i)))) + param.aR * (Rp(i)+(Rk3(i))) - param.gammaHV * (Vp(i)+(Vk3(i))) * (Hp(i)+(Hk3(i))) - param.bHF * (Fp(i)+(Fk3(i))) * (Hp(i)+(Hk3(i))));
    Ik4(i) = dt*(param.gammaHV * (Vp(i)+(Vk3(i))) * (Hp(i)+(Hk3(i))) - param.bIE * (Ep(i)+(Ek3(i))) * (Ip(i)+(Ik3(i))) - param.aI * (Ip(i)+(Ik3(i))));
	Mk4(i) = dt*((param.bMD * D(i-1) + param.bMV * (Vp(i)+(Vk3(i)))) * (1 - (Mp(i)+(Mk3(i)))) - param.aM * (Mp(i)+(Mk3(i))));
	Fk4(i) = dt*(param.bF * (Mp(i)+(Mk3(i))) + param.cF * (Ip(i)+(Ik3(i))) - param.bFH * (Hp(i)+(Hk3(i))) * (Fp(i)+(Fk3(i))) - param.aF * (Fp(i)+(Fk3(i))));
	Rk4(i) = dt*(param.bHF * (Fp(i)+(Fk3(i))) * (Hp(i)+(Hk3(i))) - param.aR * (Rp(i)+(Rk3(i))));
    Ek4(i) = dt*(param.bEM * (Mp(i)+(Mk3(i))) * (Ep(i)+(Ek3(i))) - param.bEI * (Ep(i)+(Ek3(i))) * (Ip(i)+(Ik3(i))) + param.aE * (1 - (Ep(i)+(Ek3(i)))));
	Pk4(i) = dt*(param.bPM * (Mp(i)+(Mk3(i))) * (Pp(i)+(Pk3(i))) + param.aP * (1 - (Pp(i)+(Pk3(i)))));
	Ak4(i) = dt*(param.bA * (Pp(i)+(Pk3(i))) - param.gammaAV * (Sp(i)+(Sk3(i))) * (Ap(i)+(Ak3(i))) * (Vp(i)+(Vk3(i))) - param.aA * (Ap(i)+(Ak3(i))));
	Sk4(i) = dt*(param.r * (Pp(i)+(Pk3(i))) * (1 - (Sp(i)+(Sk3(i)))));
    
    
    Vc(i) = Vc(i-1) + (1/6)*(Vk1(i) + 2*Vk2(i) + 2*Vk3(i) + Vk4(i));
    Hc(i) = Hc(i-1) + (1/6)*(Hk1(i) + 2*Hk2(i) + 2*Hk3(i) + Hk4(i));
    Ic(i) = Ic(i-1) + (1/6)*(Ik1(i) + 2*Ik2(i) + 2*Ik3(i) + Ik4(i));
	Mc(i) = Mc(i-1) + (1/6)*(Mk1(i) + 2*Mk2(i) + 2*Mk3(i) + Mk4(i));
	Fc(i) = Fc(i-1) + (1/6)*(Fk1(i) + 2*Fk2(i) + 2*Fk3(i) + Fk4(i));
	Rc(i) = Rc(i-1) + (1/6)*(Rk1(i) + 2*Rk2(i) + 2*Rk3(i) + Rk4(i));
    Ec(i) = Ec(i-1) + (1/6)*(Ek1(i) + 2*Ek2(i) + 2*Ek3(i) + Ek4(i));
	Pc(i) = Pc(i-1) + (1/6)*(Pk1(i) + 2*Pk2(i) + 2*Pk3(i) + Pk4(i));
	Ac(i) = Ac(i-1) + (1/6)*(Ak1(i) + 2*Ak2(i) + 2*Ak3(i) + Ak4(i));
	Sc(i) = Sc(i-1) + (1/6)*(Sk1(i) + 2*Sk2(i) + 2*Sk3(i) + Sk4(i));
    
    D(i) = 1 - Hc(i) - Rc(i) - Ic(i);
    
end
%% Solve ode15s solution
for j = 2:tEnd+1
    [tempVec(j,:)] = SolveODESystem (tempVec(j-1,:));
end
%% Create Variables to create figures
Vcf=zeros(tEnd,1);
Hcf=zeros(tEnd,1);
Icf=zeros(tEnd,1);
Mcf=zeros(tEnd,1);
Fcf=zeros(tEnd,1);
Rcf=zeros(tEnd,1);
Ecf=zeros(tEnd,1);
Pcf=zeros(tEnd,1);
Acf=zeros(tEnd,1);
Scf=zeros(tEnd,1);

%% Compare ode15s solution with the implicit Euler solution
figure()
for k = 1:tEnd+1
    Vcf(k) = Vc((k-1)*1000+1);
end
   subplot(2,5,1); plot(tBegin:1:tEnd,Vcf,'ok',tBegin:1:tEnd,tempVec(:,1))
   title('V')
    
for k = 1:tEnd+1
    Hcf(k) = Hc((k-1)*1000+1);
end
   subplot(2,5,2);  plot(tBegin:1:tEnd,Hcf,'ok',tBegin:1:tEnd,tempVec(:,2))
   title('H')
   
for k = 1:tEnd+1
    Icf(k) = Ic((k-1)*1000+1);
end
   subplot(2,5,3); plot(tBegin:1:tEnd,Icf,'ok',tBegin:1:tEnd,tempVec(:,3))
   title('I')
    
for k = 1:tEnd+1
    Mcf(k) = Mc((k-1)*1000+1);
end
   subplot(2,5,4);  plot(tBegin:1:tEnd,Mcf,'ok',tBegin:1:tEnd,tempVec(:,4))
   title('M')
   
for k = 1:tEnd+1
    Fcf(k) = Fc((k-1)*1000+1);
end
   subplot(2,5,5); plot(tBegin:1:tEnd,Fcf,'ok',tBegin:1:tEnd,tempVec(:,5))
   title('F')
    
for k = 1:tEnd+1
    Rcf(k) = Rc((k-1)*1000+1);
end
   subplot(2,5,6);  plot(tBegin:1:tEnd,Rcf,'ok',tBegin:1:tEnd,tempVec(:,6))
   title('R')
   
for k = 1:tEnd+1
    Ecf(k) = Ec((k-1)*1000+1);
end
   subplot(2,5,7); plot(tBegin:1:tEnd,Ecf,'ok',tBegin:1:tEnd,tempVec(:,7))
   title('E')
    
for k = 1:tEnd+1
    Pcf(k) = Pc((k-1)*1000+1);
end
   subplot(2,5,8);  plot(tBegin:1:tEnd,Pcf,'ok',tBegin:1:tEnd,tempVec(:,8))
   title('P')
   
for k = 1:tEnd+1
    Acf(k) = Ac((k-1)*1000+1);
end
   subplot(2,5,9); plot(tBegin:1:tEnd,Acf,'ok',tBegin:1:tEnd,tempVec(:,9))
   title('A')
    
for k = 1:tEnd+1
    Scf(k) = Sc((k-1)*1000+1);
end
   subplot(2,5,10);  plot(tBegin:1:tEnd,Scf,'ok',tBegin:1:tEnd,tempVec(:,10))
   title('S')
   
figure(2)
plot(t,Sc)

figure(3)
plot(tBegin:1:tEnd,tempVec(:,10))