package comp559.fluid;

import org.jblas.DoubleMatrix;

/**
 * @author litalien
 */
public abstract class LinearSolver {

    /** Iteration number */
    public int iteration;

    /** Maximum number of iterations */
    public int maxIter;

    /** Norm squared of residual r = Ax - b */
    public double resNorm2;

    /** Max error for algorithm */
    public double maxError;

    /** Default constructor */
    public LinearSolver() {
        this.iteration = 0;
        this.resNorm2 = Double.NaN;
    }

    public abstract void solve( DoubleMatrix A, DoubleMatrix b, DoubleMatrix x );
}
