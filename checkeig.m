% Checking eigenvalues for the circulant kernels of the unfolded 3D lattice 
clc, close all
dim = 6; dtotal = dim^3; k = 1/2/sqrt(2); a = k*exp(i*pi/4); b = k*exp(-i*pi/4);
cv1 = [ a b zeros(1, dim-2) b a zeros(1, dtotal-dim-2)];
cv2 = [zeros(1,dim^2) -b a zeros(1, dim-2) -a b zeros(1, dtotal-dim^2-dim-2)];
cv3 = [-a b zeros(1, dim-2) -b a zeros(1, dtotal-dim-2)];
cv4 = [zeros(1,dim^2) b a zeros(1, dim-2) a b zeros(1, dtotal-dim^2-dim-2)];
f = fft(cv1);f = [f; fft(cv2)];f =[f; fft(cv3)];f = [f; fft(cv4)];
% check against standard diagonalizer
c1 = toeplitz( [cv1(1), fliplr( cv1(2:end) )], cv1 );
c2 = toeplitz( [cv2(1), fliplr( cv2(2:end) )], cv2 );
c3 = toeplitz( [cv3(1), fliplr( cv3(2:end) )], cv3 );
c4 = toeplitz( [cv4(1), fliplr( cv4(2:end) )], cv4 );
e = eig(c1);e = [e, eig(c2)];e =[e, eig(c3)];e = [e, eig(c4)];
% Construct composite kernels
ek = f(1,:).^2+f(3,:).^2;ek = [ek; f(3,:).*(f(1,:)+f(4,:))];ek = [ek;f(3,:).^3+f(4,:).^2];
% ploting section
% Chek eigenvalues from fft routine 
colorstring = 'bgry';
figure(1), cla
hold on
for i = 1:4
  plot(f(i, :),'.', 'Color', colorstring(i))
end
hold off
figure(2), cla
hold on
x = logspace(-1,2);
for i = 1:4
    plot( abs(f(i,:))/4, 'Color', colorstring(i)), grid on
  %plot( log( conj(f(i, :)).*f(i, :) ), 'Color', colorstring(i)), grid on
end
hold off 
% Chek eigenvalues  from std eig routine
figure(3), cla
hold on
for i = 1:4
  plot(e(:, i),'.', 'Color', colorstring(i))
end
hold off
figure(4), cla
hold on
for i = 1:4
  plot( abs(e(:, i))/4 , 'Color', colorstring(i)), grid on
end
hold off 
% check composite kernels
figure(5), cla
hold on
for i = 1:3
  plot(ek(i, :),'.', 'Color', colorstring(i))
end
hold off
figure(6), cla
hold on
for i = 1:3
  plot(  abs(ek(i, :))/4 , 'Color', colorstring(i)), grid on
end
hold off