package comp559.fluid;

import org.jblas.*;
import static org.jblas.FloatMatrix.*;

/**
 * Implementation of a standard Conjugate Gradient iterative solver for linear systems
 * Solves Ax = b
 *
 * @author litalien
 */
public class ConjugateGradient extends LinearSolver {
    /**
     * Sets algorithm parameters for convergence
     * @param maxIter
     * @param maxError
     */
    public void init( int maxIter, float maxError ) {
        this.convergence = 0;
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
        float alpha, beta;
        FloatMatrix r0 = zeros(b.length), Aq = zeros(b.length);
        // CG iterate
        while (resNorm2 > maxError || iteration < maxIter) {
            A.mmuli(q, Aq);
            alpha = r.dot(r) / q.dot(Aq);
            x.addi(q.mul(alpha));
            r.subi(Aq.mul(alpha), r0);
            resNorm2 = r0.norm2();
            if (resNorm2 < maxError) {
                convergence = iteration;
                break;
            }
            beta = r0.dot(r0) / r.dot(r);
            r0.addi(q.mul(beta), q);
            FloatMatrix tmp = r;
            r = r0;
            r0 = tmp;
            ++iteration;
        }
        convergence = iteration;
    }
}
