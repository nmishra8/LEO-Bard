clear; close all;
%app_name = {'BACKPROP','BFS','BLACKSCHOLES_1M','BODYTRACK','FACESIM',...
%            'FERRET','HEARTWALL','HOTSPOT','JACOBI','LAVAMD','LEUKOCYTE',...
%            'LUD','SHA','STREAM','X264_ducks'};
% app_name = {'BACKPROP','BFS','BLACKSCHOLES_1M','BLACKSCHOLES_10M','BLACKSCHOLES_10M_1','BODYTRACK','FACESIM','FERRET',...
%             'HEARTWALL', 'HOTSPOT','JACOBI','KMEANS','LAVAMD','LUD','SHA',...
%             'STREAM','STREAM_threads','X264_ducks','X264_native','X264_phases'};

app_name = {'BACKPROP','BFS','BLACKSCHOLES_1M','BODYTRACK','FACESIM','FERRET',...
            'HEARTWALL', 'HOTSPOT','JACOBI','KMEANS','LAVAMD','LEUKOCYTE','LUD','SHA',...
            'SRAD','STREAM','STREAM_threads','X264_ducks'};
        
m  = length(app_name);
pp = 4;
for i = 1:m
    file_name = ['../../odroid/',app_name{i},'.results'];
    fid = fopen(file_name);
    c = textscan(fid,'%s %f %f %f %f','HeaderLines', 1 );    
    fclose(fid);
    c{1,3} = c{1,3}/max(c{1,3});
    c{1,4} = c{1,4}/max(c{1,4});
    rate(:,i) = c{1,3}; 
    power(:,i) = c{1,4}; 
    config = c{1,1};
end
Z{1} = rate;
Z{2} = power;
X    = c{1,2};

%%
% ratee=repmat(rate,10,1);
% powere=repmat(power,10,1);
% Xe=repmat(c{1,2},10,1);
% Z{1} = ratee;
% Z{2} = powere;
% X    = Xe;
%%
%% Parameters
[n] = size(Z{1},1);
accuracy = zeros(m,2);

%% Sample
numSamples = 20;
id1 = 1:ceil(n/numSamples):n; % points uniform over 1:1024 
id2 = randperm(n); id2 = id2(1:numSamples); %random points 

%% Run
for Y_nameId = 1:2,
    for i = 1:m
        ZZ = Z;
%         ZZ{1}(:,end) = Z{1}(:,i);
%         ZZ{2}(:,end) = Z{2}(:,i);
        [ acc, w_pred,online,offline ] = splitEM( X,ZZ,Y_nameId,id1,i );
        accuracy(i,Y_nameId) = acc;
        wl{Y_nameId,i} = w_pred;
        
        accuracy_online(i,Y_nameId) = accuracy_rss(ZZ{Y_nameId}(:,i), online);
        accuracy_offline(i,Y_nameId) = accuracy_rss(ZZ{Y_nameId}(:,i), offline);
        wl_online{Y_nameId,i} = online;
        wl_offline{Y_nameId,i} = offline;
    end
end

mean(max(accuracy,0))
mean(max(accuracy_offline,0))
mean(max(accuracy_online,0))

% Plot
for i = 1:m    
    figure(3);
    subplot(pp,ceil(m/pp),i);
    hold on;
    plot(Z{1}(:,i)); 
    plot(cell2mat(wl{1,i}),'r');
    hold off;
    title(['Rate - ',app_name{i},' : ',sprintf('%5.2f',accuracy(i,1)) ]); 

%     figure(2);
%     subplot(pp,ceil(m/pp),i);
%     hold on;
%     plot(Z{2}(:,i)); 
%     %plot(cell2mat(wl{2,i}),'r');
%     hold off;
%     title(['Power - ',app_name{i},' : ',sprintf('%5.2f',accuracy(i,2))]);
end

% % Print config file
% true_conf.cores = c{1}; 
% true_conf.freq  = c{2};
% est_conf.cores = c{1}; 
% est_conf.freq  = c{2};
% online_conf.cores = c{1}; 
% online_conf.freq  = c{2};
% offline_conf.cores = c{1}; 
% offline_conf.freq  = c{2};
% for i = 1:m      
%     true_conf.perf  = Z{1}(:,i); 
%     true_conf.power = Z{2}(:,i);
%     T = struct2table(true_conf);
%     name = strcat('Poet_Config/true_',app_name(i),'.txt');
%     writetable(T,name{1},'Delimiter',' ');
% 
%     est_conf.perf  = cell2mat(wl{1,i}); 
%     est_conf.power = cell2mat(wl{2,i});
%     T = struct2table(est_conf);
%     name = strcat('Poet_Config/est_',app_name(i),'.txt');
%     writetable(T,name{1},'Delimiter',' ');
%     
%     online_conf.perf  = wl_online{1,i}; 
%     online_conf.power = wl_online{2,i};
%     T = struct2table(online_conf);
%     name = strcat('Poet_Config/online_',app_name(i),'.txt');
%     writetable(T,name{1},'Delimiter',' ');
%     
%     offline_conf.perf  = wl_offline{1,i}; 
%     offline_conf.power = wl_offline{2,i};
%     T = struct2table(offline_conf);
%     name = strcat('Poet_Config/offline_',app_name(i),'.txt');
%     writetable(T,name{1},'Delimiter',' ');
% end

%% Nuclear norm
% addpath('MatrixCompletion/');
% lambda_tol = 10^(-6);
% tol        = 10^(-6);
% display    = 1;
% N          = 10;
% mode       = 'nuclear';
% A          = Z{1,1};
% for i = 1:m
% A(:,i) = A(:,i)/max(A(:,i));
% end
% 
% for i = 1:m
%     % Parameters
%     accuracy = zeros(m,2);
%     % Sample
%     numSamples = 20;
%     id1 = 1:ceil(n/numSamples):n; % points uniform over 1:1024 
%     id2 = randperm(n); id2 = id2(1:(n-numSamples)); %random points 
%     B = ones(size(A));
%     B(id2,i) = 0;
%     [CompletedMat, ier] = MatrixCompletion(A.*B, B,N, mode, lambda_tol, tol, display); 
%     wl_nuclear(:,i) = CompletedMat(:,i);
%     accuracy_norm(i) = 1 - norm(CompletedMat(:,i) - A(:,i))/norm(A(:,i));    
% end
% 
% figure(3);
% for i = 1:m    
%     subplot(pp,ceil(m/pp),i);
%     hold on;
%     plot(wl_nuclear(:,i));
%     plot(A(:,i),'r');    
%     title(['Rate - ',app_name{i},' : ',sprintf('%5.2f',accuracy_norm(i)) ]); 
%     hold off;
% end
