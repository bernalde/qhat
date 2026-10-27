using LinearAlgebra, SparseArrays, Arpack

"""
    operator_error_norm(E; dense_threshold=1024, tol=1e-12, maxiter=1000)

Estimate the operator 2-norm (largest singular value), including for nonnormal
unitary differences. Large matrices use Arpack's implicit normal-operator route;
they are not densified. Convergence and singular-vector residuals are checked.
This floating-point estimate is not a certified upper bound.
"""
function operator_error_norm(E::AbstractMatrix; dense_threshold::Integer=1024,
                             tol::Float64=1e-12, maxiter::Integer=1000)
    dense_threshold >= 0 || throw(ArgumentError("dense_threshold must be nonnegative"))
    isfinite(tol) && tol > 0 || throw(ArgumentError("tol must be positive and finite"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    isempty(E) && return 0.0
    entries = issparse(E) ? nonzeros(E) : E
    all(isfinite, entries) || throw(ArgumentError("operator entries must be finite"))
    peak = isempty(entries) ? 0.0 : maximum(abs, entries)
    iszero(peak) && return 0.0
    # ARPACK cannot request one Ritz pair for the smallest matrix dimensions.
    if maximum(size(E)) <= dense_threshold || minimum(size(E)) <= 2
        return opnorm(Matrix(E))
    end

    # svds works with E'E and can round tiny squared singular values to zero.
    # Scaling avoids that absolute floor for small Trotter errors.
    scaled = E / peak
    factor, nconv, _, _, _ = Arpack.svds(
        scaled; nsv=1, ritzvec=true, tol=tol, maxiter=Int(maxiter))
    nconv >= 1 && length(factor.S) == 1 || error("Operator norm iteration did not converge")
    sigma = only(factor.S)
    isfinite(sigma) && sigma > 0 || error("Invalid largest singular value for nonzero operator")
    u, v = factor.U[:, 1], factor.V[:, 1]
    residual = max(norm(scaled * v - sigma * u), norm(scaled' * u - sigma * v))
    allowance = max(10 * tol, 100 * eps(Float64)) * max(1.0, sigma)
    isfinite(residual) && residual <= allowance || error("Operator norm singular-vector residual is too large")
    value = peak * sigma
    isfinite(value) || error("Operator norm estimate overflowed")
    return value
end
