"""
    build_matrix(nx,ny,Lx,Ly)

"""
function build_matrix(nx,ny,Lx,Ly)
    size = nx*ny
    A = zeros(nx*ny,nx*ny)
    dx,dy = Lx/(nx-1),Ly/(ny-1)
    for i in 1:nx, j in 1:ny
        Adiag = (j-1)*nx + i
    
        if i == 1
            A[Adiag,Adiag] = 1
        elseif i == nx
            A[Adiag,Adiag] = 1
        elseif j == 1
            A[Adiag,Adiag] = 1
        elseif j == ny
            A[Adiag,Adiag] = 1
        else
            a0=-2/dx^2-2/dy^2
            a1=1/dx^2
            a2=1/dx^2
            a3=1/dy^2
            a4=1/dy^2

            A[Adiag,Adiag]=a0
            A[Adiag,Adiag-1]=a1
            A[Adiag,Adiag+1]=a2
            A[Adiag,Adiag-nx]=a4
            A[Adiag,Adiag+nx]=a3
        end
    end
    return sparse(A)
end

function source_term(nx,ny,Lx,Ly)
    dx,dy = Lx/(nx-1),Ly/(ny-1)
    Q_in = zeros(nx,ny)
    u_a = zeros(nx,ny)
    for i in 1:nx, j in 1:ny
        u_a[i,j] = sin(pi*(i-1)*dx/Lx)*sin(2*pi*(j-1)*dy/Ly)
        Q_in[i,j] = 4*pi^2*sin(pi*(i-1)*dx/Lx)*sin(2*pi*(j-1)*dy/Ly)/Ly^2 + pi^2*sin(pi*(i-1)*dx/Lx)*sin(2*pi*(j-1)*dy/Ly)/Lx^2
    end
    Q_in_flat = zeros(nx*ny)
    
    for i in 1:nx, j in 1:ny
            Q_in_flat[nx*(j-1)+i] = Q_in[i,j]
    end
    
    return u_a,Q_in_flat
end

function solve_matrix(nx,ny,A,Q_in_flat)
    x = A\-Q_in_flat
    history = 0
    return x,history
end

function wrap(nx,ny,x)
    x_wrapped = zeros(nx,ny)
    for i in 1:nx, j in 1:ny
        x_wrapped[i,j] = x[i+ny*(j-1)]
    end
    return x_wrapped
end

function solve_matrix_cg(nx,ny,A,Q_in_flat,tol)
    x,history = cg(A,-Q_in_flat;reltol=tol,log=true)
    return x,history
end

function l2_norm(nx,ny,u,u_a)
    l2 = 0
    u_a_sum = 0
    for i in 1:nx, j in 1:ny
        l2 += (u[i,j]-u_a[i,j])^2
        u_a_sum += (u_a[i,j])^2
    end
    return sqrt(l2/u_a_sum)
end