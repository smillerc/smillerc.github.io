@def title = "Writing a Jacobian Free Newton Krylov (JFNK) solver"
@def hascode = true
@def date = "1/10/2022"
@def tags = ["julia", "implicit methods", "numerics", "jfnk"]

# Writing a Jacobian Free Newton Krylov (JFNK) solver in Julia



## Theory

$$
\boldsymbol{J}v = \frac{\boldsymbol{F}(T + \epsilon v) - \boldsymbol{F}(T)}{\epsilon}
$$


S.Y. Kadioglu, D.A. Knoll, "A fully second order implicit/explicit time integration technique for hydrodynamics plus nonlinear heat conduction problems", Journal of Computational Physics 229 (2010) 3237–3249. [link](http://dx.doi.org/10.1016/j.jcp.2009.12.039)

# Implementation

```julia
abstract type AbstractJFNKSolver end

struct JFNKSolver1D{T <: AbstractFloat} <: AbstractJFNKSolver
    R::Vector{T}      # Residual results
    R_pert::Vector{T} # Perturbed residual results
    v::Vector{T}      # Krylov vector
    x_pert::Vector{T} # Perturbed solution vector/array
    dim::Int          # Problem dimension, e.g. for (M,N) this is M*N
end

"""1D constructor"""
function JFNKSolver(x::AbstractVector{T}) where {T <: AbstractFloat}
    len = length(x)
    BLAS.set_num_threads(Base.Threads.nthreads())
    return JFNKSolver1D(
        similar(x),
        similar(x),
        similar(x),
        similar(x),
        len)
end
```

Approximating the effect of the Jacobian-vector product $\boldsymbol{J}v$
```julia
"""Approximate the effect of the Jacobian-vector product `Jv`"""
function Jvf!(Jv::AbstractVector{T}, v::AbstractVector{T})
    global krylov_iter += 1
    # Jv is updated in-place
    # v is the krylov vector space from the GMRES solution
    b = sqrt(eps())

    v_L2norm = norm(v, 2)

    psum = 0.0
    for Idx in LinearIndices(x)
        psum = psum + b * (1.0 + abs(x[Idx]))
    end

    # Calculate the perturbation ϵ
    if v_L2norm > eps()
        ϵ = psum / (solver.dim * v_L2norm) 
    else
        ϵ = psum / solver.dim  
    end

    for Idx in LinearIndices(x)
        solver.x_pert[Idx] = x[Idx] + ϵ * v[Idx]
    end

    f(solver.R_pert, solver.x_pert, p)

    for Idx in LinearIndices(Jv)
        Jv[Idx] = (solver.R_pert[Idx] - solver.R[Idx]) / ϵ
    end
end
```

### Putting it all together


```julia
"""Jacobian-Free Newton Krylov Method"""
function JFNK!(solver::AbstractJFNKSolver, 
               x::AbstractVector{T}, 
               f, p, tolerance::Real, 
               max_iter::Int; verbose=false) where {T <: AbstractFloat}
    
    M = solver.dim
    
    fill!(solver.v, 0)
    fill!(solver.R, 0)
    fill!(solver.R_pert, 0)

    # Get the initial residual value
    f(solver.R, x, p)

    global krylov_iters_per_newton = zeros(Int64, 0)
    global krylov_iter = 0
    global newton_iter = 1

    """Approximate the effect of the Jacobian-vector product `Jv`"""
    function Jvf!(Jv::AbstractVector{T}, v::AbstractVector{T})
        global krylov_iter += 1
        # Jv is updated in-place
        # v is the krylov vector space from the GMRES solution
        b = sqrt(eps())

        v_L2norm = norm(v, 2)

        psum = 0.0
        for Idx in LinearIndices(x)
             psum = psum + b * (1.0 + abs(x[Idx]))
        end

        # Calculate the perturbation ϵ
        if v_L2norm > eps()
            ϵ = psum / (solver.dim * v_L2norm) 
        else
            ϵ = psum / solver.dim  
        end

        for Idx in LinearIndices(x)
             solver.x_pert[Idx] = x[Idx] + ϵ * v[Idx]
        end

        f(solver.R_pert, solver.x_pert, p)

        for Idx in LinearIndices(Jv)
             Jv[Idx] = (solver.R_pert[Idx] - solver.R[Idx]) / ϵ
        end
    end

    # Operator defining the function that estimtes the Jacobian vector product Jv, which
    # is done in the Jvf! function
    Jv_func = LinearMap(Jvf!, M; ismutating=true)
        
    # These need to be global to be seen in the while loop "soft" scope
    global resid_L2_init = norm(solver.R, 2)
    global resid_L2_norm = resid_L2_init
    global convergence = typemax(T)

    # Iterate until convergence or the iteration limit
    while true
        # Find the Krylov vectors
        δx, stats = dqgmres(Jv_func, solver.R, verbose=0, history=true)
        append!(krylov_iters_per_newton, krylov_iter)

        # update the solution
        for Ind in LinearIndices(x)
            x[Ind] = x[Ind] - δx[Ind]
        end

        f(solver.R, x, p) # new residual
        
        global resid_L2_norm = norm(solver.R, 2)
        global convergence = resid_L2_norm / resid_L2_init
        global newton_iter += 1
        global ave_krylov_iter = sum(krylov_iters_per_newton) / length(krylov_iters_per_newton)
        
        if newton_iter >= max_iter
            @error "JFNK max iteration reached max_iter=($max_iter)"
            return convergence, newton_iter, ave_krylov_iter
        end

        if resid_L2_norm < resid_L2_init * tolerance break end
    end

    return convergence, newton_iter, ave_krylov_iter
end
```