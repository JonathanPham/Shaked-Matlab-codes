function [k_e, f_e, f_g, f_h]=linhex(e,elementtype,order)
global f  g  h  C  IEN  nDoF  nDim  nNodesElement...
    nEdgesElement  elementType  Coord NBC  EBC;
%allocate
k_e = zeros(nNodesElement*nDoF,nNodesElement*nDoF);
f_e = zeros(nNodesElement*nDoF,1);
f_h = zeros(nNodesElement*nDoF,1);
%f_g = zeros(nNodesElement*nDoF,1);
%quadrature
if (strcmpi(elementtype,'hex'))
    nPoints=3; % change this to be dependent on element order
    switch nPoints
        case 1
            q = 0;
            w = 2;
            
        case 2
            q = sqrt(1/3)*[ -1, 1 ]';
            w = [ 1, 1 ]';
            
        case 3
            q = sqrt(3/5)*[ -1, 0, 1 ]';
            w = 1/9*[ 5, 8, 5 ]';
        otherwise
            error('Choose 1-3 gaussian quadrature points');
    end
elseif (strcmpi(elementtype,'tet'))
    nPoints=5; %for domain
    w=1/6*[-4/5, 0.45, 0.45, 0.45, 0.45]';
    q=[1/4 1/2 1/6 1/6 1/6]';
    r=[1/4 1/6 1/2 1/6 1/6]';
    s=[1/4 1/6 1/6 1/2 1/6]';
else
    error('unsuported element type');
end
%material property matrix D_e
E=C(e,1);
v=C(e,2);
D_e=zeros(6,6);
for i=1:6
    if i<4
        D_e(i,i)=1-v;
    else
        D_e(i,i)=(1-2*v)/2;
    end
end
D_e(1,2)=v;
D_e(1,3)=v;
D_e(2,3)=v;
D_e(2,1)=v;
D_e(3,1)=v;
D_e(3,2)=v;
D_e=D_e*E/(1+v)/(1-2*v);

% compute k_e, f_e
if (strcmpi(elementtype,'hex'))
    for i=1:nPoints
        xsi=q(i);
        for j=1:nPoints
            eta=q(j);
            for k=1:nPoints
                zeta=q(k);
                weight=w(i)*w(j)*w(k);
                [N, B, je,f]=SampleElementDomain(xsi,eta,zeta,e,elementtype,order);
                k_e=k_e+B'*D_e*B*je*weight;
                %f_e=f_e+N'*f(e,:)'*je*weight;
                f_e=f_e+N'*f'*je*weight;
            end
        end
    end
elseif (strcmpi(elementtype,'tet'))
    for i=1:nPoints
        [N, B, je,f]=SampleElementDomain(q(i),r(i),s(i),e,elementtype,order);
        k_e=k_e+B'*D_e*B*je*w(i);
        %f_e=f_e+N'*f(e,:)'*je*w(i);
        f_e=f_e+N'*f'*je*w(i);
    end
end
% %compute f_h - make general to include types and orders

if (strcmpi(elementtype,'hex'))
    for k=1:nEdgesElement
        if NBC(e,k)~=0
            for i=1:nPoints
                xsi=q(i);
                for j=1:nPoints
                    eta=q(j);
                    weight=w(i)*w(j);
                    [N, h, je]=SampleElementEdge(xsi,eta,e,k,elementtype,order); %define this function
                    f_h=f_h+weight*N'*h*je;
                end
            end
        end
    end
elseif (strcmpi(elementtype,'tet'))
    nPoints=3;
    w=[1 1 1]'/6;
    q=[2/3 1/6 1/6]';
    r=[1/6 2/3 1/6]';
    for k=1:nEdgesElement
        if NBC(e,k)~=0
            for i=1:nPoints
                [N, h, je]=SampleElementEdge(q(i),r(i),e,k,elementtype,order); %define this function
                f_h=f_h+w(i)*N'*h*je;
            end
        end
    end
end

% compute f_g
if EBC==1
    x = Coord(IEN(:,e),1);
    y = Coord(IEN(:,e),2);
    z=Coord(IEN(:,e),3);
    L=Params.L;
    t=Params.t;
    c=Params.c;
    g(IEN(:,e),1)=cos(pi*x/2/L);
    g(IEN(:,e),2)=cos(pi*y/2/c);
    g(IEN(:,e),3)=cos(pi*z/2/t);
end
g_e=g(IEN(:,e),:)';
f_g=-k_e*g_e(:);
end