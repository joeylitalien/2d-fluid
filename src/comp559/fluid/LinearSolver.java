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
    public double resNorm2;

    /** Max error for algorithm */
    public double maxError;

    // public int IX( int i, int j, int N ) { return i*(N+2) + j; }

    /** Default constructor */
    public LinearSolver() {
        this.iteration = 0;
        this.resNorm2 = Double.NaN;
    }

//    public void updatePressures( FloatMatrix x, float[] p, int N ) {
//        for (int i = 1; i <= N; i++) {
//            for (int j = 1; j <= N; j++) {
//                p[IX(i, j, N)] = x.get(i-1 + N*(j-1));
//            }
//        }
//    }

    public abstract void solve( FloatMatrix A, FloatMatrix b, FloatMatrix x );
}
