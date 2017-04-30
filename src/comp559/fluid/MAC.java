package comp559.fluid;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.LinkedList;
import java.util.List;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2f;
import javax.vecmath.Tuple2f;
import javax.vecmath.Vector2f;

import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;

import org.jblas.*;
import static org.jblas.FloatMatrix.*;
import static org.jblas.Decompose.*;
import java.util.Arrays;

/**
 * Eulerian fluid simulation class
 * @author litalien
 */
public class MAC {

    final static int DIM = 2;

    /** Velocity, staggered, packed (first index is the dimension) */
    public float[][] U0;

    /** Temporary velocity variable */
    private float[][] U1;

    /** Temperature (packed) */
    public float[] temperature0;

    /** Temperature (packed) temporary variable */
    private float[] temperature1;

    private IntParameter Nval = new IntParameter( "MAC Grid size (Reset)", 16, 8, 256 );

    /** Number of grid cells */
    public int N = 16;

    /** Dimension of each grid cell */
    public float dx = 1;

    /** Time elapsed in the fluid simulation */
    public double elapsed;

    /** Sources of heat and cold */
    public List<Source> sources = new LinkedList<Source>();

    /** Coefficient matrix for pressure solve */
    public FloatMatrix A;

    /** Incomplete Cholesky factorization of coefficient matrix for pressures */
    public FloatMatrix IC;

    /** Modified Incomplete Cholesky factorization */
    public FloatMatrix MIC;

    /** Numerical solver for pressure */
    public int solver = 0;
    String[] solvers = {
            "Gauss-Seidel Relaxation",
            "Successive Over-Relaxation",
            "Conjugate Gradient",
            "Incomplete Cholesky Conjugate Gradient",
            "Modified Incomplete Cholesky Conjugate Gradient (Not Implemented)"
    };

    /** Numerical solver iterations */
    public int solverIter = 0;

    /** Mean square error */
    public float solverMSE;
    public float simTime = 60;
    public int time = 0;
    public FloatMatrix groundTruth = new FloatMatrix(N*N, 600);
    public FloatMatrix solverPressure = zeros(N*N, 600);

    /**
     * Compute the index
     * @param i
     * @param j
     * @return the index in the flattened array
     */
    public int IX( int i, int j ) { return i*(N+2) + j; }

    /**
     * Assemble coefficient matrix for pressure solve
     */
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

    /**
     * Assemble Incomplete Cholesky factorization matrix
     */
    public void assembleCholesky() {
        IC = cholesky(A);
        for (int i = 0; i < N*N; i++) {
            for (int j = 0; j < N*N; j++) {
                if (A.get(i, j) == 0) {
                    IC.put(i, j, 0);
                }
            }
        }
        IC.muli(IC.transpose());

        // TODO: Implement modified version MIC(0)
    }

    /**
     * Initialize memory
     */
    public void setup() {
        elapsed = 0;
        time = 0;
        N = Nval.getValue();
        dx = 1.0f / N;
        int np2s = (N+2)*(N+2);
        U0 = new float[2][np2s];
        U1 = new float[2][np2s];
        temperature0 = new float[np2s];
        temperature1 = new float[np2s];
        solverIter = 0;
        A = zeros(N, N);
        IC = zeros(N, N);
        MIC = zeros(N, N);
        solverMSE = 0;
        assembleSparseMatrix();
        assembleCholesky();
    }

    /**
     * Gets the velocity at the given point using interpolation
     * @param x
     * @param vel
     */
    public void getVelocity( Tuple2f x, Tuple2f vel ) {
        getVelocity( x, U0, vel );
    }

    /**
     * Gets the velocity in the provided velocity field at the given point using interpolation
     * @param x
     * @param U
     * @param vel
     */
    private void getVelocity( Tuple2f x, float[][] U, Tuple2f vel ) {
        vel.x = interpolateVelU( x, U[0] );
        vel.y = interpolateVelV( x, U[1] );
    }

    /**
     * Interpolates pressures
     * @param x
     * @param s
     * @return interpolated value
     */
    public float interpolate( Tuple2f x, float[] s ) {
        int i,j;
        float wx, wy;
        i = (int) Math.floor( x.x / dx - 0.5 );
        j = (int) Math.floor( x.y / dx - 0.5 );
        wx = x.x / dx - 0.5f - i;
        wy = x.y / dx - 0.5f - j;
        float val = ( (i>=0 && j>=0 && i<=N+1 && j<=N+1 )     ? s[IX(i,j)]*(1-wx)*(1-wy) : 0 ) +
                ( (i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )     ? s[IX(i+1,j)]*wx*(1-wy) : 0 ) +
                ( (i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )     ? s[IX(i,j+1)]*(1-wx)*wy : 0 ) +
                ( (i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 ) ? s[IX(i+1,j+1)]*wx*wy : 0 );
        return val;
    }

    /**
     * Interpolates velocity field (x-component)
     * @param x
     * @param s
     * @return
     */
    public float interpolateVelU( Tuple2f x, float[] s ) {
        int i = (int) Math.floor( x.x / dx - 1.0 );
        int j = (int) Math.floor( x.y / dx - 0.5 );
        float wx = x.x / dx - 1.0f - i;
        float wy = x.y / dx - 0.5f - j;
        float val = ( (i>=0 && j>=0 && i<=N+1 && j<=N+1 )     ? s[IX(i,j)]*(1-wx)*(1-wy) : 0 ) +
                ( (i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )     ? s[IX(i+1,j)]*wx*(1-wy) : 0 ) +
                ( (i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )     ? s[IX(i,j+1)]*(1-wx)*wy : 0 ) +
                ( (i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 ) ? s[IX(i+1,j+1)]*wx*wy : 0 );
        return val;
    }

    /**
     * Interpolates velocity field (y-component)
     */
    public float interpolateVelV( Tuple2f x, float[] s ) {
        int i = (int) Math.floor( x.x / dx - 0.5 );
        int j = (int) Math.floor( x.y / dx - 1.0 );
        float wx = x.x / dx - 0.5f - i;
        float wy = x.y / dx - 1.0f - j;
        float val = ( (i>=0 && j>=0 && i<=N+1 && j<=N+1 )     ? s[IX(i,j)]*(1-wx)*(1-wy) : 0 ) +
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
    public void traceParticle( Point2f x0, float h, Point2f x1 ) {
        traceParticle( x0, U0, h, x1 );
    }

    /**
     * Performs a simple particle trace using Forward Euler
     * @param x0
     * @param U
     * @param h
     * @param x1
     */
    private void traceParticle( Point2f x0, float[][] U, float h, Point2f x1 ) {
        Vector2f vel = new Vector2f();
        x1.set( x0 );
        getVelocity(x1, U, vel);
        vel.scale(h);
        x1.add( vel );
    }

    /**
     * Add external forces using mouse
     */
    private Point2f mouseX1 = new Point2f();
    private Point2f mouseX0 = new Point2f();
    private void addMouseForce( float[][] U, float dt ) {
        Vector2f f = new Vector2f();
        f.sub( mouseX1, mouseX0 );
        float d = mouseX0.distance(mouseX1);
        if ( d < 1e-6 ) return;
        f.scale( mouseForce.getFloatValue() );
        // go along the path of the mouse!
        Point2f x = new Point2f();
        int num = (int) (d/dx + 1);
        for ( int i = 0; i <= num; i++ ) {
            x.interpolate(mouseX0, mouseX1,(float)i / num );
            addForce( U, dt, x, f );
        }
        mouseX0.set( mouseX1 );
    }

    /**
     * Sets the mouse location in the fluid for doing dragging interaction
     * @param x0 previous location
     * @param x1 current location
     */
    public void setMouseMotionPos( Point2f x0, Point2f x1 ) {
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
    private void addForce( float[][] U, float dt, Tuple2f x, Tuple2f f ) {
        addVelU( U[0], dt, x, f.x );
        addVelV( U[1], dt, x, f.y );
    }

    /**
     * Adds a force (x-component)
     * @param S
     * @param dt
     * @param x
     * @param a
     */
    private void addVelU( float[] S, float dt, Tuple2f x, float a ) {
        int i = (int) Math.floor( x.x / dx - 1.0f );
        int j = (int) Math.floor( x.y / dx - 0.5f );
        float wx = x.x / dx - 1.0f - i;
        float wy = x.y / dx - 0.5f - j;
        if ( i>=0 && j>=0 && i<=N+1 && j<=N+1 )          S[IX(i,j)]     += (1-wx)*(1-wy)* dt * a;
        if ( i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )      S[IX(i+1,j)]   += wx*(1-wy) * dt * a;
        if ( i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )      S[IX(i,j+1)]   += (1-wx)*wy * dt * a;
        if ( i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 )  S[IX(i+1,j+1)] += wx*wy * dt * a;
    }

    /**
     * Adds a force (y-component)
     * @param S
     * @param dt
     * @param x
     * @param a
     */
    private void addVelV( float[] S, float dt, Tuple2f x, float a ) {
        int i = (int) Math.floor( x.x / dx - 0.5f );
        int j = (int) Math.floor( x.y / dx - 1.0f );
        float wx = x.x / dx - 0.5f - i;
        float wy = x.y / dx - 1.0f - j;
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
    private void addSource( float[] S, float dt, Tuple2f x, float a ) {
        int i = (int) Math.floor( x.x / dx - 0.5f );
        int j = (int) Math.floor( x.y / dx - 0.5f );
        float wx = x.x / dx - 0.5f - i;
        float wy = x.y / dx - 0.5f - j;
        if ( i>=0 && j>=0 && i<=N+1 && j<=N+1 )          S[IX(i,j)]     += (1-wx)*(1-wy)* dt * a;
        if ( i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )      S[IX(i+1,j)]   += wx*(1-wy) * dt * a;
        if ( i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )      S[IX(i,j+1)]   += (1-wx)*wy * dt * a;
        if ( i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 )  S[IX(i+1,j+1)] += wx*wy * dt * a;
    }


    /**
     * Gets the average temperature of the continuum
     * @return temp
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
    private void diffuse( int b, float[] x, float[] x0, float diff, float dt ) {
        // Diffusion amount
        float a = diff * dt * N * N;
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

    /**
     * Diffuses velocity field (x-component)
     * @param x
     * @param x0
     * @param diff
     * @param dt
     */
    private void diffuseVelU( float[] x, float[] x0, float diff, float dt ) {
        // Diffusion amount
        float a = diff * dt * N * N;
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

    /**
     * Diffuses velocity field (y-component)
     * @param x
     * @param x0
     * @param diff
     * @param dt
     */
    private void diffuseVelV( float[] x, float[] x0, float diff, float dt ) {
        // Diffusion amount
        float a = diff * dt * N * N;
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
    private void advect( int b, float[] d, float [] d0, float[][] v, float dt ) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                // Get center of cell
                float x = (i + 0.5f) * dx;
                float y = (j + 0.5f) * dx;
                Point2f x0 = new Point2f(x, y);
                Point2f x1 = new Point2f();
                // Trace backward using Forward-Euler with negative time step
                // Interpolate pressure
                traceParticle(x0, v, -dt, x1);
                d[IX(i, j)] = interpolate(x1, d0);
            }
        }
        setBounds(d);
    }

    /**
     * Advects velocity field (x-component)
     * @param d
     * @param d0
     * @param u
     * @param dt
     */
    private void advectVelU( float[] d, float[] d0, float[][] u, float dt ) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                // Get u component
                float x = (i + 1.0f) * dx;
                float y = (j + 0.5f) * dx;
                Point2f x0 = new Point2f(x, y);
                Point2f x1 = new Point2f();
                // Trace backward using Forward-Euler with negative time step
                // Interpolate velocity u
                traceParticle(x0, u, -dt, x1);
                d[IX(i, j)] = interpolateVelU(x1, d0);
            }
        }
        setBoundsVelU(d);
    }

    /**
     * Advects velocity field (y-component)
     * @param d
     * @param d0
     * @param u
     * @param dt
     */
    private void advectVelV( float[] d, float[] d0, float[][] v, float dt ) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                // Get v component
                float x = (i + 0.5f) * dx;
                float y = (j + 1.0f) * dx;
                Point2f x0 = new Point2f(x, y);
                Point2f x1 = new Point2f();
                // Trace backward using Forward-Euler with negative time step
                // Interpolate velocity v
                traceParticle(x0, v, -dt, x1);
                d[IX(i, j)] = interpolateVelV(x1, d0);
            }
        }
        setBoundsVelV(d);
    }

    /**
     * Reshape vector from stacked-row to stack-column
     * @param x
     */
    private void row2col( float[] x ) {
        float[] x0 = x.clone();
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                x[i*N + j] = x0[i + j*N];
            }
        }
    }

    /**
     * Reshape vector from stacked-column to stack-row
     * @param x
     */
    private void col2row( float[] x ) {
        float[] x0 = x.clone();
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                x[i + j*N] = x0[i*N + j];
            }
        }
    }

    /**
     * Remove exterior boundary to only deal with inner grid
     * @param d
     * @param div
     * @return
     */
    private float[] removeBoundary( float[] d, float[] div ) {
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
    private void poissonSolve( float[] u, float[] v, float[] p, float[] div ) {
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
        final float maxError = error.getFloatValue();

        // Gauss-Seidel relaxation
        if (solver == 0) {
            successiveOverRelaxation(p, div, 1);
            solverIter = iterations.getValue();
            if (time < 600) {
                FloatMatrix P = new FloatMatrix(p);
                groundTruth.putColumn(time, P);
            }
        }
        // Successive Over Relaxation (SOR)
        if (solver == 1) {
            float omega = omg.getFloatValue();
            successiveOverRelaxation(p, div, omega);
            solverIter = iterations.getValue();
            if (time < 600) {
                FloatMatrix P = new FloatMatrix(p);
                groundTruth.putColumn(time, P);
            }
        }
        // Conjugate Gradient
        else if (solver == 2) {
            ConjugateGradient cg = new ConjugateGradient();
            cg.init(maxIter, maxError);
            linearSolve(cg, p, div);
            solverIter = cg.iteration;
        }
        // Incomplete Cholesky Factorization Preconditioner
        else if (solver == 3) {
            PreconditionedCG pcg = new PreconditionedCG();
            pcg.init(maxIter, maxError, IC);
            linearSolve(pcg, p, div);
            solverIter = pcg.iteration;
        }
        // Modified Incomplete Cholesky MICCG(0) Preconditioner
        else if (solver == 4) {
            // TODO: Implement this
            PreconditionedCG pcg = new PreconditionedCG();
            pcg.init(maxIter, maxError, IC);
            linearSolve(pcg, p, div);
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
    private void linearSolve( LinearSolver solver, float[] p, float[] div ) {
        // Remove boundaries since we only want to solve for inner grid
        float[] d = new float[N*N];
        removeBoundary(d, div);
        // Put in column-wise format
        row2col(d);
        // Solver Ax = d
        FloatMatrix b = new FloatMatrix(d); b.transpose();
        FloatMatrix x = new FloatMatrix(N*N);
        solver.solve(A, b, x);
        // Convert back to array type
        float[] p0 = x.toArray();
        // Reshape back to row-wise format
        col2row(p0);
        // Store into matrix for MSE
        if (time < 600) {
            FloatMatrix P = new FloatMatrix(p0);
            solverPressure.putColumn(time, P);
        }
        // Update pressures
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                p[IX(i, j)] = p0[(i-1)*N + j - 1];
            }
        }
        setBounds(p);
    }

    /**
     * Successive over-relaxation algorithm
     * When omega = 1, this is just Gauss-Seidel
     * @param p
     * @param div
     * @param omega
     */
    private void successiveOverRelaxation( float[] p, float[] div, float omega ) {
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
    private void setBounds( float[] x ) {
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

    /**
     * Sets bound on velocity field by setting velocities to zero at boundaries (x-component)
     * @param x
     */
    private void setBoundsVelU( float[] x ) {
        for (int j = 1; j <= N; j++) {
            x[IX(0, j)] = 0.f;
            x[IX(N, j)] = 0.f;
        }
    }
    /**
     * Sets bound on velocity field by setting velocities to zero at boundaries (y-component)
     * @param x
     */
    private void setBoundsVelV( float[] x ) {
        for (int i = 1; i <= N; i++) {
            x[IX(i,0)] = 0.f;
            x[IX(i, N)] = 0.f;
        }
    }

    /**
     * Advances the state of the fluid by one time step
     */
    public void step() {
        float dt = stepSize.getFloatValue();

        addMouseForce(U0, dt);

        // Use sources to change temperatures in the temperature0 scalar field
        for (Source s : sources) {
            addSource(temperature0, dt, s.location, s.amount);
        }

        // Use temperature scalar field to apply buoyancy forces to the velocity field
        float buoyancyF = buoyancy.getFloatValue();
        float refT = (float) getReferenceTemperature();
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                float deltaT = refT - temperature0[IX(i, j)];
                U0[1][IX(i, j)] += buoyancyF * deltaT * dt;
            }
        }

        // Perform velocity step
        float visc = viscosity.getFloatValue();
        diffuseVelU(U1[0], U0[0], visc, dt);
        diffuseVelV(U1[1], U0[1], visc, dt);
        poissonSolve(U1[0], U1[1], U0[0], U0[1]);
        advectVelU(U0[0], U1[0], U1, dt);
        advectVelV(U0[1], U1[1], U1, dt);
        poissonSolve(U0[0], U0[1], U1[0], U1[1]);

        // Perform scalar (density) step
        float diff = diffusion.getFloatValue();
        diffuse(0, temperature1, temperature0, diff, dt);
        advect(0, temperature0, temperature1, U0, dt);

        // Step forward in time
        elapsed += dt;
        time++;

        // Compute MSE
        if (time == 600 && solver != 0) {
            float norm = 0;
            int timeDomain = 600;
            for (int j = 0; j < timeDomain; j++) {
                FloatMatrix m = groundTruth.getColumn(j).sub(solverPressure.getColumn(j));
                norm += m.norm2();
            }
            solverMSE = norm / timeDomain;
            System.out.println("MSE = " + solverMSE);
        }
    }

    /** Simulation parameters */
    private DoubleParameter viscosity = new DoubleParameter( "Viscosity", 1e-6, 1e-8, 1 );
    private DoubleParameter diffusion = new DoubleParameter( "Diffusion", 1e-6, 1e-8, 1 );
    private DoubleParameter buoyancy = new DoubleParameter( "Buoyancy", 0.1, -1, 1 );
    private IntParameter iterations = new IntParameter( "Number of iterations", 30, 1, 1000 );
    public DoubleParameter omg = new DoubleParameter( "Over-relaxation", 1.5, 0.0, 2.0 );
    private DoubleParameter error = new DoubleParameter( "Error threshold", 1e-6, 1e-10, 1e-1 );
    private DoubleParameter mouseForce = new DoubleParameter( "Mouse force", 1e2, 1, 1e3 );
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