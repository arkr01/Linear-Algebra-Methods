function A=lcd(beta ,gamma ,N)
    % Simple code to generate test matrices for this question and maybe
    % some later assignments.
    %
    % Usage:
    %   To generate symmetric positive definite, call with beta=gamma=0.
    %   To generate nonsymmetric call, e.g., with beta=gamma=1/2.
    % 
    % Note: output matrix is of size N^2-by-N^2.

    ee=ones(N,1);
    a=4; b=-1-gamma; c=-1-beta; d=-1+beta; e=-1+gamma;
    t1=spdiags([c*ee,a*ee,d*ee],-1:1,N,N);
    t2=spdiags([b*ee,zeros(N,1),e*ee],-1:1,N,N);
    A=kron(speye(N),t1)+kron(t2,speye(N));
end