package comp559.fluid;

import org.jblas.FloatMatrix;

/**
 * @author litalien
 */
public abstract class LinearSolver {

    /** Iteration number */
    public int iteration;

    /** Maximum number of iterations */
    public int maxIter;

    /** Norm squared of residual r = Ax - b */
    public Float resNorm2;

    /** Max error for algorithm */
    public Float maxError;

    /** Convergence */
    public int convergence;

    /** Default constructor */
    public LinearSolver() {
        this.iteration = 0;
        this.resNorm2 = Float.NaN;
        this.convergence = 0;
    }

    public abstract void solve( FloatMatrix A, FloatMatrix b, FloatMatrix x );
}
