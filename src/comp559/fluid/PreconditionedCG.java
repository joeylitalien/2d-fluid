package comp559.fluid;

import org.jblas.*;
import static org.jblas.FloatMatrix.*;
import static org.jblas.Decompose.*;

/**
 * @author litalien
 */

public class PreconditionedCG extends LinearSolver {

    /** Preconditioner */
    public FloatMatrix precond;

    /**
     * Sets algorithm parameters for convergence
     * @param maxIter
     * @param maxError
     */
    public void init( int maxIter, double maxError, FloatMatrix precond ) {
        this.maxIter = maxIter;
        this.maxError = maxError;
        this.precond = precond;
    }

    /**
     * Conjugate gradient algorithm for Ax = b
     * Allows for more control over max number of iterations and error tolerance
     *
     * @param A
     * @param b
     * @return x
     */
    public void solve( FloatMatrix A, FloatMatrix b, FloatMatrix x ) {
        // Initialize residual
        FloatMatrix r = b.sub(A.mmul(x));
        FloatMatrix z = new FloatMatrix(A.rows, A.columns);
        // Apply preconditioner
        precond.mmuli(r, z);
        FloatMatrix q = z.dup();
        // Initialize alpha and beta
        float alpha = 0.0f, beta = 0.0f;
        FloatMatrix r2 = zeros(b.length), z2 = zeros(b.length), Aq = zeros(b.length);

        // CG iterate
        while (resNorm2 > maxError || iteration < maxIter) {
            A.mmuli(q, Aq);
            alpha = r.dot(z) / q.dot(Aq);
            x.addi(q.mul(alpha));
            r.subi(Aq.mul(alpha), r2);
            resNorm2 = r2.norm2();
            if (resNorm2 < maxError) break;
            precond.mmuli(r2, z2);
            beta = z2.dot(r2) / z.dot(r);
            q.addi(q.mul(beta), z2);
            FloatMatrix tmp = r;
            r = r2;
            r2 = tmp;
            iteration++;
        }
    }
}
