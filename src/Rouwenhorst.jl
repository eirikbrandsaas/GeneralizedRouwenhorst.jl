module Rouwenhorst

"""
    genrouwenhorst(rho,sigma_eps,N,T)

Julia implementation of the generalized Rouwenhorst's method to approximate 
AR(1) processes from Fella, Gallipoli, and Pan (2019).

The process follows
```math
    η_t = μ + ρ_t η_{t-1} + ϵ_t,
```
where unlike in the standard case we allow for time variying persistence,
volatility, and that the persistence parameter can be outside the unit interval (e.g., ρ=1)

See Fella, G., Gallipoli, G., & Pan, J. (2019). Markov-chain approximations for life-cycle models. Review of Economic Dynamics, 34, 183-201.
https://www.sciencedirect.com/science/article/pii/S1094202519301565

##### Arguments
- `rho::Vector{Real}` : Persistence parameter in AR(1) process
- `sigma_eps::Vector{Real}` : Standard deviation of random component of AR(1) process
- `N::Integer` : Number of nodes in grids
- `T::Integer` : Number of time periods

##### Returns

- `ygrd::Array` : `NxT` array of grids
- `trans::Array`` : `NxNxT` array of transition matrices
"""
function genrouwenhorst(rho,sigma_eps,N,T) 
    
    @assert N ≥ 2 "The state space has to have dimension N>1. Exiting."
    @assert T ≥ 2 "The time horizon has to have dimension N>1. Exiting."
    @assert length(rho)==T "Length of rho must equal time periods T"
    @assert length(sigma_eps)==T "Length of sigma_eps must equal time periods T"
    
    sigma_y = zeros(1,T);
    y_grid = zeros(N,T); 
    trans = zeros(N,N,T);
    
    # *** Step 1: construct the state space y_grid(t) in each period t.
    # Evenly-spaced N-state space over [-sqrt(N-1)*sigma_y(t),sqrt(N-1)*sigma_y(t)].
    
    # 1.a Compute unconditional variances of y(t)
    sigma_y[1] = sigma_eps[1]
    for i = 2:T 
        sigma_y[i] = sqrt(rho[i]^2*sigma_y[i-1]^2+sigma_eps[i]^2)
    end
    
    # % 1.b Construct state space
    h = 2*sqrt(N-1)*sigma_y/(N-1) # grid step
    y_grid = repeat(h,N,1) #
    y_grid[1,:]=-sqrt(N-1) * sigma_y
    y_grid = cumsum(y_grid;dims=1)
    
    # %  *** Step 2: Compute the transition matrices trans(:,:,t) from
    # %              period (t-1) to period t
    # % The transition matrix for period t is defined by parameter p(t).
    # % p(t) = 0.5*(1+rho*sigma(t-1)/sigma(t))
    
    # % Note: trans(:,:,1) is the transition matrix from y(0)=0 to any
    # % gridpoint of y_grid(1) in period 1.
    # % Any of its rows is the (unconditional) distribution in period 1.
    
    p = 1/2 # First period: p(1) = 0.5 as y(1) is white noise.  
    trans[:,:,1] = _rhmat(p,N)
    
    for j = 2:T 
            # Compute p for t>1 
            p = (sigma_y[j]+rho[j]*sigma_y[j-1])/(2*sigma_y[j] );   
            trans[:,:,j] = _rhmat(p,N);
    end
    
        
    

    return y_grid, trans
end 

function _rhmat(p,N)
    # Computes Rouwenhorst matrix as a function of p and N
        
    Pmat = zeros(N,N);
    
    # Step 2(a): get the transition matrix P1 for the N=2 case
    if N == 2
        Pmat = [p 1-p; 1-p p]
    else
        P1 = [p 1-p; 1-p p]
    
        # Step 2(b): if the number of states N>2, apply the Rouwenhorst
        # recursion to obtain the transition matrix trans
        for i = 2:N-1
            P2 = p *     [P1 zeros(size(P1,1),1); zeros(1,size(P1,2)) 0 ] + 
                 (1-p) * [zeros(size(P1,1),1) P1; 0 zeros(1,size(P1,2)) ] + 
                 (1-p) * [zeros(1,size(P1,2)) 0 ; P1 zeros(size(P1,1),1)] + 
                 p *     [0 zeros(1,size(P1,2)) ; zeros(size(P1,1),1) P1]
    
            P2[2:i,:] = 0.5*P2[2:i,:]
            
            if i==N-1
                Pmat = P2;
            else
                P1 = P2;
            end
        end # of for
    end # if N == 2

    return Pmat
end 


export genrouwenhorst
end
