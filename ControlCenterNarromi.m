% ControlCenterNarromi
% Example one: For the expression data of TFs and target genes;
% Example two: For the expression data of genes(not divided for TFs and
% targets).
% Version Data: 2012-3-25

clear
clc

%% Example one
% The TFs and targets are divided. 

load Ecoli_Data.mat;
% computing tg-tf network using narromi 
% Input: TF_list,TF_expression,TG_list,TG_expression.
lamda = 1;
alpha = 0.05; % parameter for MI correlation.
beta = 0.05; % parameter for RO.
t = 0.6; % parameter for the rate of RO in the integration.

% As an example, we choose part of the expression.
TF_list = TG_list(1:15,:);
TF_expression = TG_expression(1:15,:);
 

network = zeros(size(TF_list,1),size(TF_list,1));
network_v = network;
netsig = network;

for i=1:size(TF_list,1)
    y = TF_expression(i,:);
    X = [TF_expression(1:i-1,:);TF_expression(i+1:size(TF_expression,1),:)];
    [net,net_value,sig]=narromi(y',X',lamda,alpha, beta, t) ;
    network(i,1:i-1) = net(1:i-1);network(i,i+1:size(TF_expression,1)) = net(i:end);
    network_v(i,1:i-1) = net_value(1:i-1); network_v(i,i+1:size(TF_expression,1)) = net(i:end);
    netsig(i,1:i-1) = sig(1:i-1); netsig(i,i+1:size(TF_expression,1)) = sig(i:end);  
    i
end
 
 % Output the network 
 significance = 0.05;
 network_sig = zeros(size(netsig)) ;

 network_sig(find(netsig<=significance)) = 1 ;
 network_sig(logical(eye(size(network_sig)))) = 0;
[testfile_network]=Connect_for_cytoscape_pvalue(network_sig',network_v',netsig',TF_list,TF_list) ;

network_size=size(testfile_network,1); 
fprintf('NOTICE:\nThe Size of the Inferred Network is %d.\n',network_size);

a=testfile_network'; 
fid=fopen('network_inferred.txt','w');
fprintf(fid,'%s %s %.6f %.3e\n',a{:}) ;
fclose(fid); 
fprintf('NOTICE:\nPlease Find the Network File in the Matlab Current Folder.\n')


%% Example Two
% The TFs and targets are not divided. 

load Ecoli_Data.mat;
% computing tg-tf network using narromi 
% Input: TF_list,TF_expression,TG_list,TG_expression.
lamda = 1;
alpha = 0.05; % parameter for MI correlation.
beta = 0.05; % parameter for RO.
t = 0.6; % parameter for the rate of RO in the integration.

% As an example, we choose part of the expression.
TG_list = TG_list(1:10,:);
TG_expression = TG_expression(1:10,:);
 

network = zeros(size(TG_list,1),size(TF_list,1));
network_v = network;
netsig = network;

for i=1:size(TG_list,1)
    y = TG_expression(i,:);
    X = TF_expression;
    [net,net_value,sig]=narromi(y',X',lamda,alpha, beta, t) ;
    network(i,:) = net; network_v(i,:) = net_value; netsig(i,:) = sig;
    i
end 
 
 % Output the network 
 significance = 0.05;
 network_sig = zeros(size(netsig)) ;

 network_sig(find(netsig<=significance)) = 1 ;
 
[testfile_network]=Connect_for_cytoscape_pvalue(network_sig',network_v',netsig',TF_list,TG_list) ;

network_size=size(testfile_network,1); 
fprintf('NOTICE:\nThe Size of the Inferred Network is %d.\n',network_size);

a=testfile_network'; 
fid=fopen('network_inferred.txt','w');
fprintf(fid,'%s %s %.6f %.3e\n',a{:}) ;
fclose(fid); 
fprintf('NOTICE:\nPlease Find the Network File in the Matlab Current Folder.\n')

 
 
