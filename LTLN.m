% Solutions for the Transmission coefficient and Impedance of the
% Linear Transmission Line Network (LTLN) and the Weighted LTLN.
% (Appendix D of main report.)
clc, close all
dim = 6; dtotal = dim^3; k = 1/2/sqrt(2); a = k*exp(i*pi/4); b = k*exp(-i*pi/4);
cv1 = [ a b zeros(1, dim-2) b a zeros(1, dtotal-dim-2)];
cv2 = [zeros(1,dim^2) -b a zeros(1, dim-2) -a b zeros(1, dtotal-dim^2-dim-2)];
cv3 = [-a b zeros(1, dim-2) -b a zeros(1, dtotal-dim-2)];
cv4 = [zeros(1,dim^2) b a zeros(1, dim-2) a b zeros(1, dtotal-dim^2-dim-2)];
L1 = fft(cv1);L2 = fft(cv2);L3 = fft(cv3);L4 = fft(cv4);
LL1 = L1.*L3; LL2 = L4.*L3;
figure(1),plot(real(acosh(L1.*L3)),'*'), hold on 
plot(imag(acosh(LL1)),'.g') 
grid, xlabel 'Lattice Sites', ylabel '\gamma_1',  hold off
figure(2),plot(real(acosh(L4.*L3)),'*'), hold on 
plot(imag(acosh(LL2)),'.g') 
grid, xlabel 'Lattice Sites', ylabel '\gamma_2', hold off

A = LL1.^2 - 1; B = LL2.^2 - 1;
C1 = (L1.^2 + L3.^2)./sqrt(A); C2 = (L3.^2 + L4.^2).*sqrt(B);
D = A - B - C1.*C2;G = B.*C1;
Sol = solve('k1*x^2+k2*x+k3');
Z1 = zeros(1, dtotal);Z2 = Z1; 
for i=1:dtotal
    k1 = C1(i);k2 = D(i);k3 = G(i);
    w = eval(Sol); Z1(i) = w(1);Z2(i) = w(2);
end
figure(3), plot(real(Z1),'*'), hold on, plot(imag(Z2),'.g') 
grid, xlabel 'Lattice Sites', ylabel 'Z_1',  hold off
figure(4), plot(real(Z2),'*'), hold on, plot(imag(Z2),'.g') 
grid, xlabel 'Lattice Sites', ylabel 'Z_2',  hold off

g = 0.1; Z =1/2; % Typical test values
Tinv = [ cosh(g) -Z*sinh(g); -sinh(g)/Z cosh(g) ];
M11 = LL1 + LL2; M22 = M11; 
M12 = C1.*sqrt(A); M21 = C2./sqrt(B);
for i=1:dtotal
    M = [ M11(i) M12(i); M21(i) M22(i) ]; WM = Tinv*M; 
    w11(i) = WM(1,1); w12(i) = WM(1,2); w21(i) = WM(2,1); w22(i) = WM(2,2);
end
figure(5), plot(real(w11),'*'), hold on, plot(imag(w11),'.g') 
grid, xlabel 'Lattice Sites', ylabel 'W_{11}',  hold off
figure(6), plot(real(w12),'*'), hold on, plot(imag(w12),'.g') 
grid, xlabel 'Lattice Sites', ylabel 'W_{12}',  hold off
figure(7), plot(real(w21),'*'), hold on, plot(imag(w21),'.g') 
grid, xlabel 'Lattice Sites', ylabel 'W_{21}',  hold off
figure(8), plot(real(w22),'*'), hold on, plot(imag(w22),'.g') 
grid, xlabel 'Lattice Sites', ylabel 'W_{22}',  hold off
