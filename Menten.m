function Menten
%Introduction--------------------------------------------------------------
%                   
%   LAF692: Lauren Frisby      Date: 2016/03/09
%
%   References: 
%               [1] Logan, J.David
%                   Applied Mathematics 4th ed. Online, pg 225;
%               [2] Holzbecher, Ekkehard
%                   Environmental Modeling
%               [3] Gonzalez, David
%                   Lecture Notes and Assignment Materials 
%-----------------------------End Introduction-----------------------------

%Variables-----------------------------------------------------------------
global k1 kb k2 s0 e0 c0 p0
k1 = 10^9;                  % rate enzyme & substrate form complex
kb = 10^5;                  % rate complex is formed
k2 = 10^3;                  % rate complex forms Product
s0 = 10^(-3);               % initial substrate 
e0 = 10^(-5);               % initial enzyme
c0 = 0;                     % initial complex
p0 = 0;                     % initial product
tspan = [0 200*10^(-7)];    % time span
%-------------------------------End Variables------------------------------

%Execution-----------------------------------------------------------------
options = odeset('AbsTol', 1e-20);
[T,Y] = ode15s(@kinetics,tspan,[s0; c0; c0; s0; c0], options);
[T,Y] = ode15s(@kinetics,tspan,[s0 c0 c0 s0 c0]);
%------------------------------End Execution-------------------------------

%Figures------------------------------------------------------------------- 
figure(1) ;
plot(T,Y(:,1)/s0, T,Y(:,4)/s0, T,Y(:,2)/e0, T,Y(:,3), T,Y(:,5)/e0, '--')   %y(1)=substrate, y(2)=complex, y(3)=MM Rate, y(4)=enzyme, y(5)=product
title('Time View: Michaelis Menten Rate in Relation to Concentrations','Fontsize', 16);
legend ('Substrate', 'Enzyme', 'Complex','Michaelis Menten Rate', 'Product')
xlabel ('Time (sec)','Fontsize',14)
ylabel ('Concentration', 'Fontsize', 14)
ylim([-.6 1.6]);
%-------------------------------End Figure 1-------------------------------
figure(2) ; 
plot(T,Y(:,5)/T)                                              %y(5)=product
title('Time View: Michaelis Menten Rate', 'Fontsize', 16)
xlabel ('Time (sec)', 'Fontsize', 14)
ylabel ('Rate', 'Fontsize', 14)
ylim([-.005 .02]);
%-------------------------------End Figure 2-------------------------------
figure(3) ;
plot(T,Y(:,5)/Y(:,2), T,Y(:,5)/e0)              %y(2)=complex, y(5)=product
title('Time View: Product', 'Fontsize', 16)
legend('Product/Complex', 'Product')
xlabel ('Time(sec)','Fontsize',14)
ylabel ('Concentration', 'Fontsize', 14);
ylim([-.03 .03]);
%-------------------------------End Figure 3-------------------------------
figure(4) ;
plot(T,Y(:,4)/s0)                                              %y(4)=enzyme
title('Time View: Enzyme', 'Fontsize', 16)
legend ('Enzyme')
xlabel ('Time (sec)', 'Fontsize', 14)
ylabel ('Concentration', 'Fontsize', 14)
ylim([0 .01]);
%-------------------------------End Figure 4-------------------------------
figure(5) ;
plot(T,Y(:,2)/Y(:,4))                            %y(2)=complex, y(4)=enzyme
title('Time View: Intermediate Complex', 'Fontsize', 16)
legend('Complex')
xlabel ('Time (sec)', 'Fontsize', 14)
ylabel ('Concentration', 'Fontsize', 14);
%-------------------------------End Figure 5-------------------------------
figure(6) ;
plot(T,(Y(:,1)/s0))                                         %y(1)=substrate
title('Time View: Substrate', 'Fontsize', 16)
legend('Substrate')
xlabel ('Time (sec)', 'Fontsize', 14)
ylabel ('Concentration', 'Fontsize', 14);
%-------------------------------End Figure 6-------------------------------
%--------------------------------End Figures-------------------------------

%Function: kinetics--------------------------------------------------------
function dy = kinetics(t,y)
global k1 kb k2 s0 e0 c0 
dy = zeros(2,1);

dy(1) = -k1*e0*y(1)+(kb+k1*y(1))*y(2);                  %DS/DT, Arbitrary Substrate over Time
dy(2) = k1*e0*y(1)-(k2+kb+k1*y(1)).*y(2);               %DC/DT, Arbitrary Intermediate Complex over Time
dy(3) = (-k1*k2*(e0+c0).*y(3))./(k2+kb+k1*y(3));        %V,     Michaelis Menten Rate
dy(4) = -k1*y(1)*y(4)+(kb +k2).*(e0-y(4));              %DE/DT, Arbitrary Enzyme over Time
dy(5) = (s0-y(5)-y(1)).*k2;                             %DP/DT, Arbitrary Product over Time
%--------------------------------End Kinetics------------------------------
%---------------------------------Fin Menten-------------------------------