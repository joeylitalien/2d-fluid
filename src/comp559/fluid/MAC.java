package comp559.fluid;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.LinkedList;
import java.util.List;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;

import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;

import org.jblas.*;
import static org.jblas.DoubleMatrix.*;
import static org.jblas.Decompose.*;
import java.util.Arrays;

/**
 * Eulerian fluid simulation class
 * @author litalien
 */
public class MAC {

    final static int DIM = 2;

    /** Velocity, staggered, packed (first index is the dimension) */
    public double[][] U0;

    /** Temporary velocity variable */
    private double[][] U1;

    /** Temperature (packed) */
    public double[] temperature0;

    /** Temperature (packed) temporary variable */
    private double[] temperature1;

    private IntParameter Nval = new IntParameter( "MAC Grid size (Reset)", 16, 8, 256 );

    /** Number of grid cells */
    public int N = 16;

    /** Dimension of each grid cell */
    public double dx = 1;

    /** Time elapsed in the fluid simulation */
    public double elapsed;

    /** Sources of heat and cold */
    public List<Source> sources = new LinkedList<Source>();

    /** Coefficient matrix for pressure solve */
    public DoubleMatrix A;

    /** Incomplete Cholesky factorization of coefficient matrix for pressures */
    public DoubleMatrix C;

    /** Numerical solver for pressure */
    public int solver = 0;
    String[] solvers = {
            "Gauss-Seidel Relaxation",
            "Successive Over Relaxation",
            "Conjugate Gradient",
            "Preconditioned Incomplete Cholesky",
            "Preconditioned Modified Incomplete Cholesky (Not Implemented)"
    };

    /**
     * Compute the index
     * @param i
     * @param j
     * @return the index in the flattened array
     */
    public int IX( int i, int j ) { return i*(N+2) + j; }

    /** Assemble coefficient matrix for pressure solve */
    public void assembleSparseMatrix() {
        A = eye(N*N).mul(4.0f);
        for (int k = 0; k < N*N; k++) {
            if ((k+1) % N != 0) {
                A.put(k+1, k, -1.0f);
                A.put(k, k+1, -1.0f);
            }
        }
        for (int k = 0; k < (N-1)*N; k++) {
            A.put(k+N, k, -1.0f);
            A.put(k, k+N, -1.0f);
        }
    }

    public void assembleIncompleteCholesky() {
        C = cholesky(A);
        for (int i = 0; i < N*N; i++) {
            for (int j = 0; j < N*N; j++) {
                if (A.get(i, j) == 0) {
                    C.put(i, j, 0f);
                }
            }
        }
        C.muli(C.transpose());
//        int M = N*N;
//        C = A.dup();
//        for (int k = 0; k < M; k++) {
//            C.put(k, k, (double)Math.sqrt(C.get(k, k)));
//            for (int i = k+1; i < M; i++) {
//                if (C.get(i, k) != 0) {
//                    double a = C.get(i, k) / C.get(k, k);
//                    C.put(i, k, a);
//                }
//            }
//            for (int j = k+1; j < M; j++) {
//                for (int i = j; i < M; i++) {
//                    if (C.get(j, j) != 0) {
//                        double a = C.get(i, j) - C.get(i, k)*C.get(j, k);
//                        C.put(i, j, a);
//                    }
//                }
//            }
//        }
//        for (int i = 0; i < M; i++) {
//            for (int j = i+1; j < M; j++) {
//                C.put(i, j, 0.0f);
//            }
//        }
//        C.muli(C.transpose());
    }

    /**
     * Initialize memory
     */
    public void setup() {
        elapsed = 0;
        N = Nval.getValue();
        dx = 1.0f / N;
        int np2s = (N+2)*(N+2);
        U0 = new double[2][np2s];
        U1 = new double[2][np2s];
        temperature0 = new double[np2s];
        temperature1 = new double[np2s];
        A = zeros(N, N); C = zeros(N, N);
        assembleSparseMatrix();
        assembleIncompleteCholesky();
    }

    /**
     * Gets the velocity at the given point using interpolation
     * @param x
     * @param vel
     */
    public void getVelocity( Tuple2d x, Tuple2d vel ) {
        getVelocity( x, U0, vel );
    }

    /**
     * Gets the velocity in the provided velocity field at the given point using interpolation
     * @param x
     * @param U
     * @param vel
     */
    private void getVelocity( Tuple2d x, double[][] U, Tuple2d vel ) {
        vel.x = interpolateVelU( x, U[0] );
        vel.y = interpolateVelV( x, U[1] );
    }

    /**
     * Interpolates pressures
     * @param x
     * @param s
     * @return interpolated value
     */
    public double interpolate( Tuple2d x, double[] s ) {
        int i,j;
        double wx, wy;
        i = (int) Math.floor( x.x / dx - 0.5 );
        j = (int) Math.floor( x.y / dx - 0.5 );
        wx = x.x / dx - 0.5f - i;
        wy = x.y / dx - 0.5f - j;
        double val = ( (i>=0 && j>=0 && i<=N+1 && j<=N+1 )     ? s[IX(i,j)]*(1-wx)*(1-wy) : 0 ) +
                ( (i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )     ? s[IX(i+1,j)]*wx*(1-wy) : 0 ) +
                ( (i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )     ? s[IX(i,j+1)]*(1-wx)*wy : 0 ) +
                ( (i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 ) ? s[IX(i+1,j+1)]*wx*wy : 0 );
        return val;
    }

    public double interpolateVelU( Tuple2d x, double[] s ) {
        int i = (int) Math.floor( x.x / dx - 1.0 );
        int j = (int) Math.floor( x.y / dx - 0.5 );
        double wx = x.x / dx - 1.0f - i;
        double wy = x.y / dx - 0.5f - j;
        double val = ( (i>=0 && j>=0 && i<=N+1 && j<=N+1 )     ? s[IX(i,j)]*(1-wx)*(1-wy) : 0 ) +
                ( (i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )     ? s[IX(i+1,j)]*wx*(1-wy) : 0 ) +
                ( (i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )     ? s[IX(i,j+1)]*(1-wx)*wy : 0 ) +
                ( (i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 ) ? s[IX(i+1,j+1)]*wx*wy : 0 );
        return val;
    }

    public double interpolateVelV( Tuple2d x, double[] s ) {
        int i = (int) Math.floor( x.x / dx - 0.5 );
        int j = (int) Math.floor( x.y / dx - 1.0 );
        double wx = x.x / dx - 0.5f - i;
        double wy = x.y / dx - 1.0f - j;
        double val = ( (i>=0 && j>=0 && i<=N+1 && j<=N+1 )     ? s[IX(i,j)]*(1-wx)*(1-wy) : 0 ) +
                ( (i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )     ? s[IX(i+1,j)]*wx*(1-wy) : 0 ) +
                ( (i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )     ? s[IX(i,j+1)]*(1-wx)*wy : 0 ) +
                ( (i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 ) ? s[IX(i+1,j+1)]*wx*wy : 0 );
        return val;
    }

    /**
     * Performs a simple Forward Euler particle trace using the current velocity field
     * @param x0
     * @param h
     * @param x1
     */
    public void traceParticle( Point2d x0, double h, Point2d x1 ) {
        traceParticle( x0, U0, h, x1 );
    }

    /**
     * Performs a simple particle trace using Forward Euler
     * (Note that this could be something higher order or adative)
     * x1 = x0 + h * U(x0)
     * @param x0
     * @param U
     * @param h
     * @param x1
     */
    private void traceParticle( Point2d x0, double[][] U, double h, Point2d x1 ) {
        Vector2d vel = new Vector2d();
        x1.set( x0 );
        getVelocity(x1, U, vel);
        vel.scale(h);
        x1.add( vel );
    }

    private Point2d mouseX1 = new Point2d();

    private Point2d mouseX0 = new Point2d();

    private void addMouseForce( double[][] U, double dt ) {
        Vector2d f = new Vector2d();
        f.sub( mouseX1, mouseX0 );
        double d = mouseX0.distance(mouseX1);
        if ( d < 1e-6 ) return;
        f.scale( mouseForce.getValue() );
        // go along the path of the mouse!
        Point2d x = new Point2d();
        int num = (int) (d/dx + 1);
        for ( int i = 0; i <= num; i++ ) {
            x.interpolate(mouseX0, mouseX1,(double)i / num );
            addForce( U, dt, x, f );
        }
        mouseX0.set( mouseX1 );
    }

    /**
     * Sets the mouse location in the fluid for doing dragging interaction
     * @param x0 previous location
     * @param x1 current location
     */
    public void setMouseMotionPos( Point2d x0, Point2d x1 ) {
        mouseX0.set( x0 );
        mouseX1.set( x1 );
    }

    /**
     * Adds a force at a given point in the provided velocity field
     * @param U    velocity field
     * @param dt   time step
     * @param x    location
     * @param f    force
     */
    private void addForce( double[][] U, double dt, Tuple2d x, Tuple2d f ) {
        addVelU( U[0], dt, x, f.x );
        addVelV( U[1], dt, x, f.y );
    }

    private void addVelU( double[] S, double dt, Tuple2d x, double a ) {
        int i = (int) Math.floor( x.x / dx - 1.0f );
        int j = (int) Math.floor( x.y / dx - 0.5f );
        double wx = x.x / dx - 1.0f - i;
        double wy = x.y / dx - 0.5f - j;
        if ( i>=0 && j>=0 && i<=N+1 && j<=N+1 )          S[IX(i,j)]     += (1-wx)*(1-wy)* dt * a;
        if ( i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )      S[IX(i+1,j)]   += wx*(1-wy) * dt * a;
        if ( i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )      S[IX(i,j+1)]   += (1-wx)*wy * dt * a;
        if ( i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 )  S[IX(i+1,j+1)] += wx*wy * dt * a;
    }

    private void addVelV( double[] S, double dt, Tuple2d x, double a ) {
        int i = (int) Math.floor( x.x / dx - 0.5f );
        int j = (int) Math.floor( x.y / dx - 1.0f );
        double wx = x.x / dx - 0.5f - i;
        double wy = x.y / dx - 1.0f - j;
        if ( i>=0 && j>=0 && i<=N+1 && j<=N+1 )          S[IX(i,j)]     += (1-wx)*(1-wy)* dt * a;
        if ( i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )      S[IX(i+1,j)]   += wx*(1-wy) * dt * a;
        if ( i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )      S[IX(i,j+1)]   += (1-wx)*wy * dt * a;
        if ( i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 )  S[IX(i+1,j+1)] += wx*wy * dt * a;
    }

    /**
     * Adds the time step scaled amount to the provided scalar field at the specified location x in an interpolated manner
     * @param S    scalar field
     * @param dt   time step
     * @param x    location
     * @param a    amount
     */
    private void addSource( double[] S, double dt, Tuple2d x, double a ) {
        int i = (int) Math.floor( x.x / dx - 0.5f );
        int j = (int) Math.floor( x.y / dx - 0.5f );
        double wx = x.x / dx - 0.5f - i;
        double wy = x.y / dx - 0.5f - j;
        if ( i>=0 && j>=0 && i<=N+1 && j<=N+1 )          S[IX(i,j)]     += (1-wx)*(1-wy)* dt * a;
        if ( i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )      S[IX(i+1,j)]   += wx*(1-wy) * dt * a;
        if ( i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )      S[IX(i,j+1)]   += (1-wx)*wy * dt * a;
        if ( i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 )  S[IX(i+1,j+1)] += wx*wy * dt * a;
    }

        /**
         * Gets the average temperature of the continuum
         * (use this in computing bouyancy forces)
         * @return
         */
    public double getReferenceTemperature() {
        int count = 0;
        double referenceTemperature = 0;
        for ( int i = 1; i <= N; i++ ) {
            for ( int j = 1; j <= N; j++ ) {
                referenceTemperature += temperature0[IX(i,j)];
                count++;
            }
        }
        referenceTemperature /= count;
        return referenceTemperature;
    }

    /**
     * Diffuses density field using Gauss-Seidel relaxation iterative solver
     * @param b     bound
     * @param x     new scalar field
     * @param x0    old density field
     * @param dt    time step
     */
    private void diffuse( int b, double[] x, double[] x0, double diff, double dt ) {
        // Diffusion amount
        double a = diff * dt * N * N;
        // Get number of iterations to use
        final int nIter = iterations.getValue();
        // Diffuse grid
        for (int k = 0; k < nIter; k++) {
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    x[IX(i,j)] = ( x0[IX(i, j)] + a*(x[IX(i-1, j)] + x[IX(i+1, j)] + x[IX(i,j-1)] + x[IX(i,j+1)]) ) / (1 + 4*a);
                }
            }
            setBounds(x);
        }
    }

    private void diffuseVelU( double[] x, double[] x0, double diff, double dt ) {
        // Diffusion amount
        double a = diff * dt * N * N;
        // Get number of iterations to use
        final int nIter = iterations.getValue();
        // Diffuse grid
        for (int k = 0; k < nIter; k++) {
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    x[IX(i,j)] = ( x0[IX(i, j)] + a*(x[IX(i-1, j)] + x[IX(i+1, j)] + x[IX(i,j-1)] + x[IX(i,j+1)]) ) / (1 + 4*a);
                }
            }
            setBoundsVelU(x);
        }
    }

    private void diffuseVelV( double[] x, double[] x0, double diff, double dt ) {
        // Diffusion amount
        double a = diff * dt * N * N;
        // Get number of iterations to use
        final int nIter = iterations.getValue();
        // Diffuse grid
        for (int k = 0; k < nIter; k++) {
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    x[IX(i,j)] = ( x0[IX(i, j)] + a*(x[IX(i-1, j)] + x[IX(i+1, j)] + x[IX(i,j-1)] + x[IX(i,j+1)]) ) / (1 + 4*a);
                }
            }
            setBoundsVelV(x);
        }
    }

    /**
     * Forces (transports) scalar to follow a given scalar field
     * @param b     bound
     * @param d     new scalar field
     * @param d0    old scalar field
     * @param dt    time step
     */
    private void advect( int b, double[] d, double [] d0, double[][] v, double dt ) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                // Get center of cell
                double x = (i + 0.5f) * dx;
                double y = (j + 0.5f) * dx;
                Point2d x0 = new Point2d(x, y);
                Point2d x1 = new Point2d();
                // Trace backward using Forward-Euler with negative time step
                // Interpolate pressure
                traceParticle(x0, v, -dt, x1);
                d[IX(i, j)] = interpolate(x1, d0);
            }
        }
        setBounds(d);
    }

    private void advectVelU( double[] d, double[] d0, double[][] u, double dt ) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                // Get u component
                double x = (i + 1.0f) * dx;
                double y = (j + 0.5f) * dx;
                Point2d x0 = new Point2d(x, y);
                Point2d x1 = new Point2d();
                // Trace backward using Forward-Euler with negative time step
                // Interpolate velocity u
                traceParticle(x0, u, -dt, x1);
                d[IX(i, j)] = interpolateVelU(x1, d0);
            }
        }
        setBoundsVelU(d);
    }

    private void advectVelV( double[] d, double[] d0, double[][] v, double dt ) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                // Get v component
                double x = (i + 0.5f) * dx;
                double y = (j + 1.0f) * dx;
                Point2d x0 = new Point2d(x, y);
                Point2d x1 = new Point2d();
                // Trace backward using Forward-Euler with negative time step
                // Interpolate velocity v
                traceParticle(x0, v, -dt, x1);
                d[IX(i, j)] = interpolateVelV(x1, d0);
            }
        }
        setBoundsVelV(d);
    }

    private void row2col( double[] x ) {
        double[] x0 = x.clone();
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                x[i*N + j] = x0[i + j*N];
            }
        }
    }

    private void col2row( double[] x ) {
        double[] x0 = x.clone();
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                x[i + j*N] = x0[i*N + j];
            }
        }
    }

    private double[] removeBoundary( double[] d, double[] div ) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                d[i*N + j] = div[(i+1)*(N+2) + (j+1)];
            }
        }
        return d;
    }

    /**
     * Forces velocity field to be mass-conserving
     * @param u     velocity field (x-component)
     * @param v     velocity field (y-component)
     * @param p     pressure
     * @param div   divergence
     */
    private void poissonSolve( double[] u, double[] v, double[] p, double[] div ) {
        // Update divergences
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                div[IX(i, j)] = -dx * (u[IX(i, j)] - u[IX(i-1, j)] + v[IX(i,j)] - v[IX(i,j-1)]);
                p[IX(i, j)] = 0;
            }
        }
        setBounds(div);
        setBounds(p);

        // Solve for pressures...
        final int maxIter = iterations.getValue();
        final double maxError = error.getValue();

        // Gauss-Seidel relaxation
        if (solver == 0) {
            successiveOverRelaxation(p, div, 1);
        }
        // Successive Over Relaxation (SOR)
        if (solver == 1) {
            double omega = omg.getValue();
            successiveOverRelaxation(p, div, omega);
        }
        // Conjugate Gradient
        else if (solver == 2) {
            ConjugateGradient cg = new ConjugateGradient();
            cg.init(maxIter, maxError);
            linearSolve(cg, p, div);
        }
        // Incomplete Cholesky Factorization Preconditioner
        else if (solver == 3) {
            PreconditionedCG pcg = new PreconditionedCG();
            pcg.init(maxIter, maxError, C);
            linearSolve(pcg, p, div);
        }
        // Modified Incomplete Cholesky MICCG(0) Preconditioner
        else if (solver == 4) {
//            PreconditionedCG pcg = new PreconditionedCG();
//            pcg.init(maxIter, maxError, C);
//            linearSolve(pcg, p, div);
        }

        // Update velocities
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                u[IX(i, j)] -= (p[IX(i+1, j)] - p[IX(i, j)]) / dx;
                v[IX(i, j)] -= (p[IX(i, j+1)] - p[IX(i, j)]) / dx;
            }
        }
        setBoundsVelU(u);
        setBoundsVelV(v);
    }

    /**
     * Linear solver for discrete Poisson problem
     * @param sol   numerical linear solver
     * @param p     unknown pressures
     * @param div   divergences
     * */
    private void linearSolve( LinearSolver solver, double[] p, double[] div ) {
        // Remove boundaries since we only want to solve for inner grid
        double[] d = new double[N*N];
        removeBoundary(d, div);
        // Put in column-wise format
        row2col(d);
        // Solver Ax = d
        DoubleMatrix b = new DoubleMatrix(d); b.transpose();
        DoubleMatrix x = new DoubleMatrix(N*N);
        solver.solve(A, b, x);
        // Convert back to array type
        double[] p0 = x.toArray();
        // Reshape back to row-wise format
        col2row(p0);
        // Update pressures
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                p[IX(i, j)] = p0[(i-1)*N + j - 1];
            }
        }
        setBounds(p);
    }

    private void successiveOverRelaxation( double[] p, double[] div, double omega ) {
        // Get parameters
        final int maxIter = iterations.getValue();
        // Launch numerical solver
        for (int k = 0; k < maxIter; k++) {
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    p[IX(i, j)] = (1 - omega) * p[IX(i, j)] + omega * (div[IX(i, j)] + p[IX(i-1, j)] +
                                  p[IX(i+1, j)] + p[IX(i, j-1)] + p[IX(i, j+1)]) / 4;
                }
            }
            setBounds(p);
        }
    }

    /**
     * Sets bounds for next iteration of solver
     * @param b     bound
     * @param x     scalar field
     */
    private void setBounds( double[] x ) {
        // Cancel scalar field at the boundaries
        for (int i = 1; i <= N; i++) {
            x[IX(0,   i)] = x[IX(1, i)];
            x[IX(N+1, i)] = x[IX(N,    i)];
            x[IX(i,   0)] = x[IX(i, 1)];
            x[IX(i, N+1)] = x[IX(i,    N)];
        }
        // Deal with corners
        x[IX(0,    0)] = 0.5f * (x[IX(1,  0)] + x[IX(0,  1)]);
        x[IX(0,  N+1)] = 0.5f * (x[IX(1,N+1)] + x[IX(0,     N)]);
        x[IX(N+1,0  )] = 0.5f * (x[IX(N,   0  )] + x[IX(N+1,1)]);
        x[IX(N+1,N+1)] = 0.5f * (x[IX(N,   N+1)] + x[IX(N+1,   N)]);
    }

    private void setBoundsVelU( double[] x ) {
        for (int j = 1; j <= N; j++) {
            x[IX(0, j)] = 0.f;
            x[IX(N, j)] = 0.f;
        }
    }

    private void setBoundsVelV( double[] x ) {
        for (int i = 1; i <= N; i++) {
            x[IX(i,0)] = 0.f;
            x[IX(i, N)] = 0.f;
        }
    }

    /**
     * Advances the state of the fluid by one time step
     */
    public void step() {
        double dt = stepSize.getValue();

        addMouseForce(U0, dt);

        // Use sources to change temperatures in the temperature0 scalar field
        for (Source s : sources) {
            addSource(temperature0, dt, s.location, s.amount);
        }

        // Use temperature scalar field to apply buoyancy forces to the velocity field
        double buoyancyF = buoyancy.getValue();
        double refT = (double) getReferenceTemperature();
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                double deltaT = refT - temperature0[IX(i, j)];
                U0[1][IX(i, j)] += buoyancyF * deltaT * dt;
            }
        }

        // Perform velocity step
        double visc = viscosity.getValue();
        diffuseVelU(U1[0], U0[0], visc, dt);
        diffuseVelV(U1[1], U0[1], visc, dt);
        poissonSolve(U1[0], U1[1], U0[0], U0[1]);
        advectVelU(U0[0], U1[0], U1, dt);
        advectVelV(U0[1], U1[1], U1, dt);
        poissonSolve(U0[0], U0[1], U1[0], U1[1]);

        // Perform scalar (density) step
        double diff = diffusion.getValue();
        diffuse(0, temperature1, temperature0, diff, dt);
        advect(0, temperature0, temperature1, U0, dt);

        // Step forward in time
        elapsed += dt;
    }

    private DoubleParameter viscosity = new DoubleParameter( "Viscosity", 1e-6, 1e-8, 1 );
    private DoubleParameter diffusion = new DoubleParameter( "Diffusion", 1e-6, 1e-8, 1 );
    private DoubleParameter buoyancy = new DoubleParameter( "Buoyancy", 0.1, -1, 1 );
    private IntParameter iterations = new IntParameter( "Number of iterations", 30, 20, 200 );
    public DoubleParameter omg = new DoubleParameter( "Omega (SOR only)", 1.5, 0.0, 2.0 );
    private DoubleParameter error = new DoubleParameter( "Error threshold", 1e-6, 1e-10, 1e-1 );
    private DoubleParameter mouseForce = new DoubleParameter( "Mouse force", 1e2, 1, 1e3 );

    /** Step size of the simulation */
    public DoubleParameter stepSize = new DoubleParameter( "Step size", 0.1, 0.001, 1 );

    /**
     * Numerical solver dropdown menu
     */
    private class DropDown extends JComboBox implements ActionListener {
//        private static final long serialVersionUID = 1L;
//        private int solverNumber;
        public DropDown( String[] solvers ) {
            for (int i = 0; i < solvers.length; i++) {
                this.addItem(solvers[i]);
            }
            addActionListener( this );
        }
        @Override
        public void actionPerformed(ActionEvent e) {
            JComboBox cb = (JComboBox)e.getSource();
            solver = cb.getSelectedIndex();
        }
    }

    /**
     * Gets the parameters for the fluid
     * @return a control panel
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        VerticalFlowPanel vfp1 = new VerticalFlowPanel();
        vfp1.setBorder(new TitledBorder("Fluid Properties"));
        vfp1.add( viscosity.getSliderControls(true) );
        vfp1.add( diffusion.getSliderControls(true) );
        vfp1.add( buoyancy.getSliderControls(false) );
        vfp1.add( mouseForce.getSliderControls(true) );
        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.add(new DropDown(solvers));
        vfp2.setBorder(new TitledBorder("Numerical Solver for Pressures"));
        vfp2.add( stepSize.getSliderControls(true ) );
        vfp2.add( iterations.getSliderControls() );
        vfp2.add( error.getSliderControls(true) );
        vfp2.add( omg.getSliderControls( false ) );
        vfp2.add( Nval.getSliderControls() );
        vfp.add( vfp1.getPanel() );
        vfp.add( vfp2.getPanel() );
        return vfp.getPanel();
    }
}