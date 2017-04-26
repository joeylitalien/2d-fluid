package comp559.fluid;

import java.util.LinkedList;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2f;
import javax.vecmath.Tuple2f;
import javax.vecmath.Vector2f;

import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;

/**
 * Eulerian fluid simulation class.
 *
 * This follows closely the solver documented in J. Stam GDC 2003, specifically
 * this uses a non staggered grid, with an extra set of grid cells to help enforce
 * boundary conditions
 *
 * @author kry
 */
public class Fluid {

    final static int DIM = 2;

    /** velocity, non-staggered, packed (first index is the dimension), public for achievement checking */
    public float[][] U0;

    /** temporary velocity variable */
    private float[][] U1;

    /** temperature (packed) */
    public float[] temperature0;

    /** temperature (packed) temporary variable */
    private float[] temperature1;

    private IntParameter Nval = new IntParameter( "grid size (requires reset)", 16, 8, 256 );

    /** Number of grid cells (not counting extra boundary cells */
    public int N = 16;

    /** Dimension of each grid cell */
    public float dx = 1;

    /** time elapsed in the fluid simulation */
    public double elapsed;

    /**
     * Sources of heat and cold
     */
    public List<Source> sources = new LinkedList<Source>();

    /**
     * initialize memory
     */
    public void setup() {
        elapsed = 0;
        N = Nval.getValue();
        dx = 1.0f / N; // we choose the domain size here to be 1 unit square!

        int np2s = (N+2)*(N+2);
        U0 = new float[2][np2s];
        U1 = new float[2][np2s];
        temperature0 = new float[np2s];
        temperature1 = new float[np2s];
    }

    /**
     * Compute the index
     * @param i
     * @param j
     * @return the index in the flattened array
     */
    public int IX( int i, int j ) {
        return i*(N+2) + j;
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
        vel.x = interpolate( x, U[0] );
        vel.y = interpolate( x, U[1] );
    }

    /**
     * Interpolates the given scalar field
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
        float val = ( (i>=0 && j>=0 && i<=N+1 && j<=N+1 )          ? s[IX(i,j)]*(1-wx)*(1-wy) : 0 ) +
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
     * (Note that this could be something higher order or adative)
     * x1 = x0 + h * U(x0)
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
            x.interpolate(mouseX0,mouseX1, (float)i / num );
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
        addSource( U[0], dt, x, f.x );
        addSource( U[1], dt, x, f.y );
    }

    /**
     * Adds the time step scaled amount to the provided scalar field at the specified location x in an interpolated manner
     * @param S    scalar field
     * @param dt   time step
     * @param x    location
     * @param a    amount
     */
    private void addSource( float[] S, float dt, Tuple2f x, float a ) {
        int i = (int) Math.floor( x.x / dx - 0.5 );
        int j = (int) Math.floor( x.y / dx - 0.5 );
        float wx = x.x / dx - 0.5f - i;
        float wy = x.y / dx - 0.5f - j;
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
            setBounds(b, x);
        }
    }

    /**
     * Forces (transports) density to follow a given velocity field
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
                traceParticle(x0, v, -dt, x1);
                d[IX(i, j)] = interpolate(x1, d0);
            }
        }
        setBounds(b, d);
    }

    /**
     * Forces velocity field to be mass-conserving
     * @param u     velocity field (x component)
     * @param v     velocity field (y component)
     * @param p     pressure
     * @param div   divergence
     */
    private void poissonSolve( float[] u, float[] v, float[] p, float[] div ) {
        // Update divergence
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                div[IX(i, j)] = -0.5f * dx * (u[IX(i+1, j)] - u[IX(i-1, j)]
                                + v[IX(i,j+1)] - v[IX(i,j-1)]);
                p[IX(i, j)] = 0;
            }
        }
        setBounds(0, div);
        setBounds(0, p);

        // Gauss-Seidel relaxation for pressure
        final int nIter = iterations.getValue();
        for (int k = 0; k < nIter; k++) {
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    p[IX(i, j)] = (div[IX(i, j)] + p[IX(i-1, j)] + p[IX(i+1, j)]
                                  + p[IX(i, j-1)] + p[IX(i, j+1)]) / 4;
                }
            }
        }
        setBounds(0, p);

        // Update velocities
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                u[IX(i, j)] -= 0.5f * (p[IX(i+1, j)] - p[IX(i-1, j)]) / dx;
                v[IX(i, j)] -= 0.5f * (p[IX(i, j+1)] - p[IX(i, j-1)]) / dx;
            }
        }
        setBounds(1, u);
        setBounds(2, v);
    }

    /**
     * Sets bounds for next iteration of Gauss-Seidel solver
     * @param b     bound
     * @param x     scalar field
     */
    private void setBounds( int b, float[] x ) {
        for (int i = 1; i <= N; i++) {
            x[IX(0,   i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
            x[IX(N+1, i)] = b == 1 ? -x[IX(N,    i)] : x[IX(N,    i)];
            x[IX(i,   0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, N+1)] = b == 2 ? -x[IX(i,    N)] : x[IX(i,    N)];
        }
        x[IX(0,    0)] = 0.5f * (x[IX(1,  0)] + x[IX(0,  1)]);
        x[IX(0,  N+1)] = 0.5f * (x[IX(1,N+1)] + x[IX(0,     N)]);
        x[IX(N+1,0  )] = 0.5f * (x[IX(N,   0  )] + x[IX(N+1,1)]);
        x[IX(N+1,N+1)] = 0.5f * (x[IX(N,   N+1)] + x[IX(N+1,   N)]);
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
        diffuse(1, U1[0], U0[0], visc, dt);
        diffuse(2, U1[1], U0[1], visc, dt);
        poissonSolve(U1[0], U1[1], U0[0], U0[1]);
        advect(1, U0[0], U1[0], U1, dt);
        advect(2, U0[1], U1[1], U1, dt);
        poissonSolve(U0[0], U0[1], U1[0], U1[1]);

        // Perform scalar (density) step
        float diff = diffusion.getFloatValue();
        diffuse(0, temperature1, temperature0, diff, dt);
        advect(0, temperature0, temperature1, U0, dt);

        // Step forward in time
        elapsed += dt;
    }

    private DoubleParameter viscosity = new DoubleParameter( "viscosity", 1e-6, 1e-8, 1 );
    private DoubleParameter diffusion = new DoubleParameter( "diffusion", 1e-6, 1e-8, 1 );
    private DoubleParameter buoyancy = new DoubleParameter( "buoyancy", 0.1, -1, 1 );
    private IntParameter iterations = new IntParameter( "Gauss Seidel iterations", 30, 20, 200 );
    private DoubleParameter mouseForce = new DoubleParameter( "mouse force", 1e2, 1, 1e3 );

    /** step size of the simulation */
    public DoubleParameter stepSize = new DoubleParameter( "step size", 0.1, 0.001, 1 );

    /**
     * Get the parameters for the fluid
     * @return a control panel
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        VerticalFlowPanel vfp1 = new VerticalFlowPanel();
        vfp1.setBorder(new TitledBorder("Fluid properties"));
        vfp1.add( viscosity.getSliderControls(true) );
        vfp1.add( diffusion.getSliderControls(true) );
        vfp1.add( buoyancy.getSliderControls(false) );
        vfp1.add( mouseForce.getSliderControls(true) );
        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.setBorder(new TitledBorder("Fluid solver properties"));
        vfp2.add( stepSize.getSliderControls(true ) );
        vfp2.add( iterations.getSliderControls() );
        vfp2.add( Nval.getSliderControls() );
        vfp.add( vfp1.getPanel() );
        vfp.add( vfp2.getPanel() );
        return vfp.getPanel();
    }
}