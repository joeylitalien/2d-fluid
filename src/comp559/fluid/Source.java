package comp559.fluid;

import javax.vecmath.Point2d;

/**
 * Heating or cooling source in the fluid
 * 
 * @author kry
 */
public class Source {
    
    /** position of the heat source */
    Point2d location = new Point2d();
    
    /** heating or cooling rate */
    public double amount = 0;
    
    /** flag to denote that the mouse is over the source */
    boolean highlight = false;
    
    /**
     * Creates a source with the given amount
     * @param p
     * @param a
     */
    public Source( Point2d p, double a ) {
        location.set(p);
        amount = a;
        highlight = false;
    }
    
}