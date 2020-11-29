
% The function iprime = Gllshort(t, i) defines the differential
% equations of the synchronous machine during a
% three-phase fault. The function returns the state derivatives
% of the current.
%current vector [-id,-iq,-i0,iF,iD,iQ]
% Copyright  (c) 1998 H. Saadat,Modified by Jiang Qirong


function iprime =Gllshort(t,i)
    f=50.;  w=2.*pi*f;
    d=0.;   d=d*pi/180.;  theta=w*t+d;

%  Parameters of a 500 MVA, 30 kV Synchronous Machine
    Ld = 0.0072;   Lq = 0.0070;  L0=0.0010;  LF = 2.500;   LD = 0.0068;   LQ = 0.0016;  
    MF = 0.100;   MD = 0.0054;   MQ = 0.0026;  MR = 0.1250;
    r  = 0.002;   rF = 0.4000;   rD  = 0.015;  rQ = 0.0150;

    VF = 400;                % dc field voltage
    V = [0; 0; 0; VF; 0; 0];   % Voltage column vector
    K=3/2.;  RT2=sqrt(2.0);

    R=[r        0          0        0          0               0
       0        r          0        0          0               0
       0        0          r        0          0               0 
       0        0          0        rF         0               0
       0        0          0        0          rD              0
       0        0          0        0          0               rQ];

    L=[Ld       0          0        MF         MD         0
       0        Lq         0        0          0          MQ    
       0        0          L0       0          0          0
       K*MF     0          0        LF         MR         0
       K*MD     0          0        MR         LD         0
       0        K*MQ       0        0          0          LQ];
   
    WW=[0        -w          0        0          0               0
        w        0           0        0          0               0
        0        0           0        0          0               0 
        0        0           0        0          0               0
        0        0           0        0          0               0
        0        0           0        0          0               0];

    Li=inv(L);
    iprime=Li*V - Li*(WW*L+R)*i;
end
