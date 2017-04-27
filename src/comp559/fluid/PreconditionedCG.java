package comp559.fluid;

import org.jblas.*;
import static org.jblas.DoubleMatrix.*;
import static org.jblas.Geometry.*;

/**
 * @author litalien
 */

public class PreconditionedCG extends LinearSolver {

    /** Preconditioner */
    public DoubleMatrix precond;

    /**
     * Sets algorithm parameters for convergence
     * @param maxIter
     * @param maxError
     */
    public void init( int maxIter, double maxError, DoubleMatrix precond ) {
        this.maxIter = maxIter;
        this.maxError = maxError;
        this.precond = precond;
    }

    /**
     * Preconditioned conjugate gradient algorithm for Ax = b
     * @param A
     * @param b
     * @return x
     */
    public void solve( DoubleMatrix A, DoubleMatrix b, DoubleMatrix x ) {
        // Initialize residual
        DoubleMatrix r = b.sub(A.mmul(x));
        // Apply preconditioner
        DoubleMatrix z = precond.mmul(r);
        DoubleMatrix q = z.dup();
        // Initialize alpha and beta
        double alpha, beta;
        DoubleMatrix r0 = zeros(b.length), z0 = zeros(b.length), Aq = zeros(b.length);

        // CG iterate
        while (resNorm2 > maxError || iteration < maxIter) {
            A.mmuli(q, Aq);
            alpha = r.dot(z) / q.dot(Aq);
            x.addi(q.mul(alpha));
            r.subi(Aq.mul(alpha), r0);
            // checkForNaN(r0);
            resNorm2 = r0.norm2();
            if (resNorm2 < maxError) break;
            precond.mmuli(r0, z0);
            beta = z0.dot(r0.sub(r)) / z.dot(r);
            z0.addi(q.mul(beta), q);
            DoubleMatrix r1 = r, z1 = z;
            r = r0; z = z0;
            r0 = r1; z0 = z1;
            iteration++;
        }
    }
}
