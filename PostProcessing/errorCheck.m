function [] = errorCheck(seed)

    filename = join(['errorCheck-',num2str(seed),'.bin']);
    file = fopen(filename);
    A = fread(file,[3,Inf],'double');
    
    t = A(1,:);
    W = A(2,:);
    dE = A(3,:);

    figure(1)
    plot(t, (W-dE)./dE)
    xlabel('t (s)')
    ylabel('% error')