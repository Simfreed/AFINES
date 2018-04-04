clear
clf

eps = 0.04;
sig = 0;
b = 10;
kexv = 0.04;
start = 0.001;
spacing = 0.0001;

i = 1;
j = 1;

hold on

for sig = 0.01:0.01:0.5
    for r = start:spacing:0.01
    
        F_wca(i,j) = 24*eps*(2*(sig^12)/(r^13) - (sig^6)/(r^7));
        F_kx(i,j) = 2*kexv*r^2*((1/r) - b);
        R(i,j) = r;
        F_r12(i,j) = 12*eps*(sig^12)/(r^13);
    
        i = i + 1;
    end
    
    figure(1)
    plot(R(:,j), F_r12(:,j))
    
    figure(2)
    plot(R(:,j), F_wca(:,j))
    
    j = j+1;
end
    
%figure(1)
%plot(R, F_kx)

%figure(2)
%plot(R, F_r12)

