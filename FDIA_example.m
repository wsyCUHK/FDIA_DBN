load('D:\preliminary_contest_data\data\Binarydata14_3.mat')

maxvalue=max(max(max(x_train),max(x_test)));
minvalue=min(min(min(x_train),min(x_test)));
x_train=(x_train-minvalue)/(maxvalue-minvalue);
x_test=(x_test-minvalue)/(maxvalue-minvalue);
sizeA=floor(size(x_train,1)/10)*10;
sizeB=floor(size(x_test,1)/10)*10;

input_x=x_train(1:sizeA,:);
input_y=y_sin_train(1:sizeA);

testcase_x=x_test(1:sizeB,:);
testcase_y=y_sin_test(1:sizeB);

% x_train = double(x_train) / 255;
% y_sin_train = double(y_sin_train);
% x_test = gpuArray(double(x_test)  / 255);
% y_sin_test = gpuArray(double(y_sin_test));

% part = zeros(10, 2);
% for i = 1:10
%    [~, col, ~] = find(y_sin_train);
%    col = find(col == i);
%    part(i, :) = [col(1) col(end)];
% end

rng('default');
% model parameters
sizes = [500 500 2000];
opts.numepochs =   10;
opts.batchsize =   10;
opts.momentum  =   0;
opts.alpha     =   0.1;
opts.decay     =   0.00001;
opts.k         =   1;
% opts.part      =   part;

dbn = DBN(input_x, input_y, sizes, opts);
tic
train(dbn, input_x, input_y);
toc

%% compute most probable class label given test data
probs = dbn.predict(testcase_x, testcase_y);
[~, I] = max(probs, [], 2);
pred = bsxfun(@eq, I, 1:10);
mis = find(~all(pred == testcase_y,2));

err = length(mis) / size(testcase_y, 1);

fprintf('Classification error is %3.2f%%\n',err*100);

