package comp559.fluid;

import org.jblas.*;
import static org.jblas.FloatMatrix.*;

/**
 * Implementation of a Conjugate Gradient iterative solver for linear systems
 * Implements standard algorithm
 *
 * @author litalien
 */
public class ConjugateGradient extends LinearSolver {

    /**
     * Sets algorithm parameters for convergence
     * @param maxIter
     * @param maxError
     */
    public void init( int maxIter, double maxError ) {
        this.maxIter = maxIter;
        this.maxError = maxError;
    }

    /** Conjugate gradient algorithm for Ax = b
     * Allows for more control over max number of iterations and error tolerance
     * @param A
     * @param b
     * @return x
     */
    public void solve( FloatMatrix A, FloatMatrix b, FloatMatrix x ) {
        // Initialize residual
        FloatMatrix r = b.sub(A.mmul(x));
        FloatMatrix q = r.dup();
        // Initialize alpha and beta
        float alpha = 0.0f, beta = 0.0f;
        FloatMatrix r2 = zeros(b.length), Aq = zeros(b.length);
        // CG iterate
        while (resNorm2 > maxError || iteration < maxIter) {
            A.mmuli(q, Aq);
            alpha = r.dot(r) / q.dot(Aq);
            x.addi(q.mul(alpha));
            r.subi(Aq.mul(alpha), r2);
            resNorm2 = r2.norm2();
            if (resNorm2 < maxError) break;
            beta = r2.dot(r2) / r.dot(r);
            r2.addi(q.mul(beta), q);
            FloatMatrix tmp = r;
            r = r2;
            r2 = tmp;

            ++iteration;
        }
    }
}
