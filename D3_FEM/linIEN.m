function IEN=linIEN(Nx,Ny,Nz)
IEN = zeros(8,Nx*Ny*Nz);
for k=1:Nz
    for i=1:Ny
        for j=1:Nx
            e=j+Nx*(i-1)+Nx*Ny*(k-1);
            IEN(1,e)=j+(Nx+1)*(i-1)+(Nx+1)*(Ny+1)*(k-1);
            IEN(2,e)=j+1+(Nx+1)*(i-1)+(Nx+1)*(Ny+1)*(k-1);
            IEN(3,e)=j+1+(Nx+1)*(i)+(Nx+1)*(Ny+1)*(k-1);
            IEN(4,e)=j+(Nx+1)*(i)+(Nx+1)*(Ny+1)*(k-1);
            IEN(5,e)=j+(Nx+1)*(i-1)+(Nx+1)*(Ny+1)*(k);
            IEN(6,e)=j+1+(Nx+1)*(i-1)+(Nx+1)*(Ny+1)*(k);
            IEN(7,e)=j+1+(Nx+1)*(i)+(Nx+1)*(Ny+1)*(k);
            IEN(8,e)=j+(Nx+1)*(i)+(Nx+1)*(Ny+1)*(k);
        end
    end
end
end