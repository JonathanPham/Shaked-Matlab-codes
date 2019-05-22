function [k_e, f_e, f_g, f_h]=linhex(e,elementtype)
global f  g  h  C  IEN  nDoF  nDim  nNodesElement...
    nEdgesElement  elementType  Coord NBC;
%allocate
k_e = zeros(nNodesElement*nDoF,nNodesElement*nDoF);
f_e = zeros(nNodesElement*nDoF,1);
f_h = zeros(nNodesElement*nDoF,1);
%quadrature
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
for i=1:nPoints
    xsi=q(i);
    for j=1:nPoints
        eta=q(j);
        for k=1:nPoints
            zeta=q(k);
            weight=w(i)*w(j)*w(k);
            [N, B, je]=SampleElementDomain(xsi,eta,zeta,e,elementtype);
            k_e=k_e+B'*D_e*B*je*weight;
            f_e=f_e+N'*f(e,:)'*je*weight;
        end
    end
end

% %compute f_h - still need to define SampleElementEdge
% for k=1:nEdgesElement
%     for i=1:nPoints
%         xsi=q(i);
%         for j=1:nPoints
%             eta=q(j);
%             weight=w(i)*w(j);
%             [N, h, je]=SampleElementEdge(xsi,eta,e,k); %define this function
%             f_h=f_h+weight*N'*h*je;
%         end
%     end
% end
% compute f_g
g_e=g(IEN(i,e))';
f_g=-k_e*g_e;
end