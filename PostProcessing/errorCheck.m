function [] = errorCheck(seed)

    % Currently not in use: error calculations have been combined into output for filePlot()

    filename = join(['errorCheck-',num2str(seed),'.bin']);
    file = fopen(filename);
    A = fread(file,[4,Inf],'double');
    fclose(file);

    t = A(1,:);
    W = A(2,:);
    dE = A(3,:);
    Etot = A(4,:);
    err = (W-dE)./Etot;
    
    figure(1)
    plot(t, err * 100, '.')
    xlabel('t'), ylabel('% error')
    xlim([0, t(end)])
    
    % Using old work calculation in output.cpp:
    
    % dW = zeros(1,size(W,2));
    % dE = zeros(1,size(E,2));
    % err = zeros(1,size(t,2));
    
    % for i = 2:size(t,2)
    %     dW(i) = W(i-1) * (t(i)-t(i-1))/t(i-1);
    %     dE(i) = E(i) - E(i-1);
    %     err(i) = (dW(i) - dE(i))./E(i);
    % end
    
    % err = err * 100;
    
    % To check error calculation:
    
    % figure(2)
    % plot(t, W - dE)
    % xlabel('t'), ylabel('W - \DeltaE')
    % xlim([0, t(end)])