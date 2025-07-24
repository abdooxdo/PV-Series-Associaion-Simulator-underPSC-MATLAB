
clc,clear
%%information from the KC200GT solar array datasheet
Vocn=32.9 ;
Iscn=8.21 ;
Imp=7.61 ;
Vmp=26.3 ;
Pmax_e=Vmp*Imp;
kv=-0.123;
ki=3.18e-3;
Ns=54;
%%constants 
q=1.60217646e-19;
k=1.3806503e-23;
eg=1.12;
% extracted parameter using method described in (Green Energy and Technology) S. Sumathi, L. Ashok Kumar, P. Surekha (auth.) - Solar PV and Wind Energy Conversion Systems_ An Introduction to Theory, Modeling with MATLAB_SIMULINK, and the Role of Soft Computing Techniques
a=1.3; % diode ideality factor 
%%nominal values
Tn=25+273.15;
Gn=1000;
%%adjusting algorithm
%the model is adjusted at nominal conditions
G=1000;
T=25+273.15;
Vtn=k*Tn/q;
Vt=k*T/q;
dT=T-Tn;
I0n=(Iscn+ki*dT)/(exp((Vocn+kv*dT)/a/Ns/Vtn)-1);
I0=I0n;
Rs=0.221;
Rp=415.405;
Ipvn=8.214;
Ipv=(Ipvn+ki*dT)*G/Gn;
Ipv2=0.5*Ipv;
 %solving et IV equation
    clear V2
    clear V1
    clear Vpv
    clear I
    clear Vd
    I=0:.01:8.21;
   % ii2=0:.1:0.2*Ipv;
    V1=zeros(1,size(I,2));
    V2=zeros(1,size(I,2));
    for j=1:size(I,2)
        %solves g=I-f(I,V)=0 with newton raphson
        g(j)=Ipv-I0*(exp((V1(j)+I(j)*Rs)/Vt/Ns/a)-1)-(V1(j)+I(j)*Rs)/Rp-I(j);
        while(abs(g(j))>0.001)
        g(j)=Ipv-I0*(exp((V1(j)+I(j)*Rs)/Vt/Ns/a)-1)-(V1(j)+I(j)*Rs)/Rp-I(j);
        glin(j)=-I0/Vt/Ns/a*exp((V1(j)+I(j)*Rs)/Vt/Ns/a)-1/Rp-1;
        V1_(j)=V1(j)-g(j)/glin(j);
        V1(j)=V1_(j);
        end
        g(j)=Ipv2-I0*(exp((V2(j)+I(j)*Rs)/Vt/Ns/a)-1)-(V2(j)+I(j)*Rs)/Rp-I(j);
        if I(j)<=Ipv2   
        while(abs(g(j))>0.001) 
        g(j)=Ipv2-I0*(exp((V2(j)+I(j)*Rs)/Vt/Ns/a)-1)-(V2(j)+I(j)*Rs)/Rp-I(j);
        glin(j)=-I0/Vt/Ns/a*exp((V2(j)+I(j)*Rs)/Vt/Ns/a)-1/Rp-1;
        V2_(j)=V2(j)-g(j)/glin(j);
        V2(j)=V2_(j);
        end
        end
        %V1(j)=Ns1*a*Vt*log((Ipv-I(j)+I0)/I0);
        %if I(j)<=0.5*Ipv
        %V2(j)=Ns2*a*Vt*log((0.5*Ipv-I(j)+I0)/I0);
        %end
    end
        
       % Vd(j)=-Vt*a*log((I(j)/I0)+1)   
    if V2<=0
            Vpv=V1;
        else
            Vpv=V1+V2;
    end
    P=Vpv.*I;
% Plot the IV and PV charastiristics
figure
plot(Vpv,I);
figure
plot(Vpv,P);