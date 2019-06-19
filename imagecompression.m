%reading and converting the image
inImage=imread('cameraman.tif');
%inImage=rgb2gray(inImage);
inImageD=double(inImage);
imwrite(uint8(inImageD), 'original.jpg');
% decomposing the image using singular value decomposition
[U,S,V]=svd(inImageD);
% Using different number of singular values (diagonal of S) to compress and
% reconstruct the image
dispEr = [];
numSVals = [];

i=0;
N = 1
 % store the singular values in a temporary var
 i=i+1;
 C = S;
 % discard the diagonal values not required for compression
 C(N+1:end,:)=0;
 C(:,N+1:end)=0;
 % Construct an Image using the selected singular values
 D=U*C*V';
 images(:,:,i)=D;
 % display and compute error
 %figure;
 buffer = sprintf('Image output using %d singular values', N)
 error=sum(sum((inImageD-D).^2));
 % store vals for display
 dispEr = [dispEr; error];
 numSVals = [numSVals; N];

for N=1:2:20
 C = S;
 C(N+1:end,:)=0;
 C(:,N+1:end)=0;
 D=U*C*V';
 if N==5
    i=i+1;
    images(:,:,i)=D;
 end
 error=sum(sum((inImageD-D).^2));
 % store vals for display
 dispEr = [dispEr; error];
 numSVals = [numSVals; N];
end
for N=25:25:100
 % store the singular values in a temporary var
 i=i+1;
 C = S;
 % discard the diagonal values not required for compression
 C(N+1:end,:)=0;
 C(:,N+1:end)=0;
 % Construct an Image using the selected singular values
 D=U*C*V';
 images(:,:,i)=D;
 % display and compute error
 %figure;
 buffer = sprintf('Image output using %d singular values', N)
 error=sum(sum((inImageD-D).^2));
 % store vals for display
 dispEr = [dispEr; error];
 numSVals = [numSVals; N];
end
N=[1 5 25 50 75 100];
for i=1:6
    subplot(3,2,i);
    D=images(:,:,i);
    buffer = sprintf('Image output using %d singular values',N(i) )
    imshow(uint8(D));
    imwrite(uint8(D), sprintf('%dbw.jpg', N(i)));
    title(buffer);
end

% dislay the error graph
figure;
title('Error in compression');
plot(numSVals, dispEr);
hold on
grid on
xlabel('Number of Singular Values used');
ylabel('Error between compress and original image');







%reading and converting the image
inImage=imread('salzburg_square.png');
%inImage=rgb2gray(inImage);
inImageD=double(inImage);
imwrite(uint8(inImageD), 'original.jpg');
% decomposing the image using singular value decomposition
[U,S,V]=svd(inImageD);
% Using different number of singular values (diagonal of S) to compress and
% reconstruct the image
dispEr = [];
numSVals = [];

i=0;
N = 1
 % store the singular values in a temporary var
 i=i+1;
 C = S;
 % discard the diagonal values not required for compression
 C(N+1:end,:)=0;
 C(:,N+1:end)=0;
 % Construct an Image using the selected singular values
 D=U*C*V';
 images(:,:,i)=D;
 % display and compute error
 %figure;
 buffer = sprintf('Image output using %d singular values', N)
 error=sum(sum((inImageD-D).^2));
 % store vals for display
 dispEr = [dispEr; error];
 numSVals = [numSVals; N];

for N=1:2:20
 C = S;
 C(N+1:end,:)=0;
 C(:,N+1:end)=0;
 D=U*C*V';
 if N==5
    i=i+1;
    images(:,:,i)=D;
 end
 error=sum(sum((inImageD-D).^2));
 % store vals for display
 dispEr = [dispEr; error];
 numSVals = [numSVals; N];
end
for N=25:25:100
 % store the singular values in a temporary var
 i=i+1;
 C = S;
 % discard the diagonal values not required for compression
 C(N+1:end,:)=0;
 C(:,N+1:end)=0;
 % Construct an Image using the selected singular values
 D=U*C*V';
 images(:,:,i)=D;
 % display and compute error
 %figure;
 buffer = sprintf('Image output using %d singular values', N)
 error=sum(sum((inImageD-D).^2));
 % store vals for display
 dispEr = [dispEr; error];
 numSVals = [numSVals; N];
end

% dislay the error graph
plot(numSVals, dispEr);
legend("Low variance image", "High variance image")


figure;
N=[1 5 25 50 75 100];
for i=1:6
    subplot(3,2,i);
    D=images(:,:,i);
    buffer = sprintf('Image output using %d singular values',N(i) )
    imshow(uint8(D));
    imwrite(uint8(D), sprintf('%dbw.jpg', N(i)));
    title(buffer);
end
