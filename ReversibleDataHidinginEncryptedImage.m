%修改于Bep_发布在CSDN上的信息安全之加密域可逆信息隐藏RDH-EI(附加代码)
clc;
clear;
original_img = imread('1.bmp');
[n,m] = size(original_img); % m ,n表示图片的大小
modify_bit = 3; %嵌入信息需要的位数
figure('NumberTitle','off','Name', '3');
%块大小为8*8 每个块中嵌入一位信息
blocksize = 8;
subplot(2, 3, 1);
imshow(original_img);
title('原始图像');
original_img = double(original_img); %防止图像失真

%把原始图像像素变成二进制
original_img_bin =zeros(8,n,m);
for i=1:8
    for j=1:n
        for k=1:m
            original_img_bin(i,j,k) =bitget(original_img(j,k), i);
            original_img_bin(i,j,k) = double( original_img_bin(i,j,k));
        end
    end
end

%生成加密秘钥  
en_key =zeros(8,n,m);
for i=1:8
    en_key(i,(1:n),(1:m)) = rand(n,m) < 0.5;
end
en_key = double(en_key);

%对图像进行加密
encrypted_img_bin = zeros(8,n,m);
for i=1:8
    for j=1:n
        for k=1:m
              encrypted_img_bin(i, j, k) =  xor(original_img_bin(i, j, k),  en_key(i,j,k));  
        end
    end
end

%将加密的二进制图像变成十进制
encrypted_img =zeros(n,m);
for i=1:8
    for j=1:n
        for k=1:m
            encrypted_img(j,k) = encrypted_img(j,k) + 2^(i-1) * encrypted_img_bin(i,j,k);
        end
    end
end

subplot(2,3,2);
imshow(uint8(encrypted_img));
title('加密图像');

%生成集合概率 如果 < 0.5 为S0集合, 否则为S1集合;
set = rand(n,m) < 0.5;
%生成嵌入信息 
size = 8;
N = 6;
lim_row = fix(n/size);
lim_col = fix(m/size);
bitts = 3;
message=rand(lim_row,lim_col)<0.5;

%向加密图像中嵌入信息, 
for i=1:n/blocksize
    for j=1:m/blocksize
        if(message(i,j) == 0) %嵌入信息为0
            %将这块内的S0集合的后三位取反
            for k1=(i-1)*blocksize+1:i*blocksize
                for k2=(j-1)*blocksize+1:j*blocksize
                    if(set(k1,k2) == 0) 
                       for k3 =1:modify_bit 
                           encrypted_img_bin(k3,k1,k2) = ~encrypted_img_bin(k3,k1,k2); %取反
                       end
                    end
                end
            end
        else %嵌入信息为1
            %将这块内的S1集合的后三位取反
             for k1=(i-1)*blocksize+1:i*blocksize
                for k2=(j-1)*blocksize+1:j*blocksize
                    if(set(k1,k2) == 1) 
                       for k3 =1:modify_bit 
                           encrypted_img_bin(k3,k1,k2) = ~encrypted_img_bin(k3,k1,k2); %取反
                       end
                    end
                end
            end
        end
    end
end

%将嵌入信息的加密图像变成十进制
encrypted_img = zeros(n,m);
for i=1:8
    for j=1:n
        for k=1:m
            encrypted_img(j,k) = encrypted_img(j,k) + 2^(i-1) * encrypted_img_bin(i,j,k);
        end
    end
end

subplot(2,3,3)
imshow(uint8(encrypted_img));
title('嵌入信息的加密图像');

%解密已嵌入信息的加密图像
decrypted_img_bin = zeros(8,n,m);
for i=1:8
    for j=1:n
        for k=1:m
            decrypted_img_bin(i,j,k) = xor(encrypted_img_bin(i,j,k), en_key(i,j,k));  
        end
    end
end

%将解密后的已嵌入信息的加密图像转化为十进制
decrypted_img = zeros(n,m);
for i=1:8
    for j=1:n
        for k=1:m
            decrypted_img(j,k) = decrypted_img(j,k) + 2^(i-1) * decrypted_img_bin(i,j,k);
        end
    end
end

subplot(2,3,4)
imshow(uint8( decrypted_img));
title('解密后的嵌入信息图像');

%提取数据后恢复的图像
%利用波动公式判断每个块中嵌入的数据是0还是1
recover = zeros(n,m);
infermessage = zeros(lim_row,lim_col);
for i=1:n/blocksize
    for j=1:m/blocksize
        H0_img = zeros(n,m);
        H1_img = zeros(n,m);
        for k =(i-1)*blocksize+1:i*blocksize
            for k1=(j-1)*blocksize+1:j*blocksize
                %转化为十进制
                for k2=1:8
                    if(set(k,k1) == 0) 
                        if(k2 > modify_bit ) 
                            H0_img(k,k1) = H0_img(k,k1) + (decrypted_img_bin(k2,k,k1)) * 2 ^(k2-1);
                        else
                            H0_img(k,k1) = H0_img(k,k1) + ~(decrypted_img_bin(k2,k,k1)) * 2 ^(k2-1);
                        end %将H0集合后面三位取反
                        H1_img(k,k1) = H1_img(k,k1) + (decrypted_img_bin(k2,k,k1)) * 2 ^(k2-1); %H1集合不变
                    end
                    if(set(k,k1) == 1)
                        if(k2 > modify_bit ) 
                            H1_img(k,k1) = H1_img(k,k1) + (decrypted_img_bin(k2,k,k1)) * 2 ^(k2-1);
                        else%将H1集合后面三位取反
                            H1_img(k,k1) = H1_img(k,k1) + ~(decrypted_img_bin(k2,k,k1)) * 2 ^(k2-1);
                        end
                        H0_img(k,k1) = H0_img(k,k1) + (decrypted_img_bin(k2,k,k1)) * 2 ^(k2-1); %H0集合不变
                    end
                end
            end
        end
        f0 = 0;
        f1 = 0;
        %利用波动公式判断嵌入的数据是原图数据是H1集合还是H0集合
        for k2=(i-1)*blocksize+2:i*blocksize-1
            for k3=(j-1)*blocksize+2:j*blocksize-1
               f0 = f0 + abs( H0_img(k2,k3)-( H0_img(k2-1,k3) + H0_img(k2,k3-1) + H0_img(k2+1,k3) + H0_img(k2,k3+1))/4 );
               f1 = f1 + abs(H1_img(k2,k3)-(H1_img(k2-1,k3) + H1_img(k2,k3-1) + H1_img(k2+1,k3) + H1_img(k2,k3+1))/4 );
            end
         
        end
        if(f0 < f1)
                for k1 = (i-1)*blocksize+1:i*blocksize
                    for k2=(j-1)*blocksize+1:j*blocksize
                        recover(k1,k2) = H0_img(k1,k2);
                    end
                end
                infermessage(i,j) = 0;
        else
              for k1 = (i-1)*blocksize+1:i*blocksize
                    for k2=(j-1)*blocksize+1:j*blocksize
                        recover(k1,k2) = H1_img(k1,k2);
                    end
              end     
              infermessage(i,j) = 1;
        end
    end
end

subplot(2,3,5)
imshow(uint8( recover));
title('提取数据后恢复的图像');

err_block = zeros(n,m);

for i=1:n/blocksize
    for j=1:m/blocksize
        if(infermessage(i,j) == message(i,j)) %提取正确的块
            err_block((i-1)*blocksize+1:i*blocksize,(j-1)*blocksize+1:j*blocksize) = 255;
        else
            err_block((i-1)*blocksize+1:i*blocksize,(j-1)*blocksize+1:j*blocksize) = 0;
        end
    end
end

subplot(2,3,6)
imshow(uint8( err_block));
title('错误块');

%原始图与从嵌入数据恢复的图像之间均方误差
s = 0;
for i=1:n
    for j=1:m
        s = s + ( double(original_img(i,j)) - decrypted_img(i,j))^2;
    end
end
mse1 = s/(n*m);
psnr1= 10*log10(255^2/mse1)%峰值信噪比 评价图像的指标
%PSNR值越大，就代表失真越少。

%原图像与提取数据恢复的图像之间均方误差
s = 0;
for i=1:n
    for j=1:m
        s = s + ( double(original_img(i,j)) - recover(i,j))^2;
    end
end
mse1 = s/(n*m);
psnr1= 10*log10(255^2/mse1)%峰值信噪比 评价图像的指标
%PSNR值越大，就代表失真越少。