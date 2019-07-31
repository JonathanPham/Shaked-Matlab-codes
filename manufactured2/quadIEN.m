function IEN=quadIEN(Nx,Ny,Nz,order)
switch order
    case 1
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
    case 2
        IEN = zeros(27,Nx*Ny*Nz);
        Npx=2*Nx+1;
        Npy=2*Ny+1;
        for k=1:Nz
            for i=1:Ny
                for j=1:Nx
                    e=j+Nx*(i-1)+Nx*Ny*(k-1);
                    M=2*(Npx*(i-1)+Npy*Npx*(k-1)+(j-1));
                    IEN(1,e)=M+1;
                    IEN(2,e)=M+3;
                    IEN(3,e)=M+2*Npx+3;
                    IEN(4,e)=M+2*Npx+1;
                    IEN(5,e)=M+1+2*Npy*Npx;
                    IEN(6,e)=M+3+2*Npy*Npx;
                    IEN(7,e)=M+2*Npx+3+2*Npy*Npx;
                    IEN(8,e)=M+2*Npx+1+2*Npy*Npx;
                    IEN(9,e)=M+2;
                    IEN(10,e)=M+Npx+3;
                    IEN(11,e)=M+2*Npx+2;
                    IEN(12,e)=M+Npx+1;
                    IEN(13,e)=M+2+2*Npy*Npx;
                    IEN(14,e)=M+Npx+3+2*Npy*Npx;
                    IEN(15,e)=M+2*Npx+2+2*Npy*Npx;
                    IEN(16,e)=M+Npx+1+2*Npy*Npx;
                    IEN(17,e)=M+1+Npy*Npx;
                    IEN(18,e)=M+3+Npy*Npx;
                    IEN(19,e)=M+2*Npx+3+Npy*Npx;
                    IEN(20,e)=M+2*Npx+1+Npy*Npx;
                    IEN(21,e)=M+Npx+2;
                    IEN(22,e)=M+Npx+2+2*Npy*Npx;
                    IEN(23,e)=M+2+Npy*Npx;
                    IEN(24,e)=M+2*Npx+2+Npy*Npx;
                    IEN(25,e)=M+Npx+1+Npy*Npx;
                    IEN(26,e)=M+Npx+3+Npy*Npx;
                    IEN(27,e)=M+Npx+2+Npy*Npx;
                end
            end
        end
    otherwise
            error('choose order 1 or 2');
end
end