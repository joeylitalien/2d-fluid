package comp559.fluid;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2f;

import mintools.swing.VerticalFlowPanel;

/**
 * @author litalien
 */
public class TestSystems {

    /** Copy of current state of the system */
    public MAC fluid;
    public int N;
    public float dx;
    public List<Source> sources;

    /** Test cases names */
    public String[] tests = {
            "Centered positive (1 source)",
            "Centered negative (1 source)",
            "Alternate diagonal (2 sources)",
            "Symmetric square (4 sources)"
    };

    /**
     * Creates a new test system
     * @param fluid
     */
    public TestSystems( MAC fluid ) {
        this.fluid = fluid;
        N = fluid.N;
        dx = fluid.dx;
        sources = fluid.sources;
    }

    /**
     * Test generation button
     */
    private class TestButton extends JButton implements ActionListener {
        private static final long serialVersionUID = 1L;
        private int testNumber;
        public TestButton( String name, int testNumber ) {
            super( name );
            this.testNumber = testNumber;
            addActionListener( this );
        }
        @Override
        public void actionPerformed(ActionEvent e) {
            createSystem(this.testNumber);
        }
    }

    /**
     * Gets the control panel for setting different systems
     * @return control panel
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.setBorder( new TitledBorder("Fluid Test Systems"));
        for (int i = 0; i < tests.length; i++) {
            vfp.add(new TestButton(tests[i], i));
        }
        return vfp.getPanel();
    }

    /**
     * Creates new simple test system
     * @param which
     */
    public void createSystem( int which ) {
        // Clear scene
        while ( !fluid.sources.isEmpty() ) {
            fluid.sources.remove(0);
        }
        // Single source
        if ( which == 0 ) {
            float x = N/2;
            float y = N/2;
            Point2f loc = new Point2f((x+1)*dx, (y+1)*dx);
            Source s = new Source(loc, 1.0f);
            s.highlight = true;
            fluid.sources.add(0, s);
        }
        if ( which == 1 ) {
            float x = N/2;
            float y = N/2;
            Point2f loc = new Point2f((x+1)*dx, (y+1)*dx);
            Source s = new Source(loc, -1.0f);
            s.highlight = true;
            fluid.sources.add(0, s);
        }
        // Opposite source
        else if ( which == 2 ) {
            float x = 0.75f*N;
            float y = 0.25f*N;
            Point2f loc1 = new Point2f((x+1)*dx, (y+1)*dx);
            Point2f loc2 = new Point2f((y+1)*dx, (x+1)*dx);
            Source s1 = new Source(loc1, -1.0f);
            Source s2 = new Source(loc2, 1.0f);
            s1.highlight = true;
            s2.highlight = true;
            fluid.sources.add(0, s1);
            fluid.sources.add(0, s2);
        }
        // Square opposite sources
        else if ( which == 3 ) {
            float x = 0.75f*N;
            float y = 0.25f*N;
            Point2f loc1 = new Point2f((x+1)*dx, (y+1)*dx);
            Point2f loc2 = new Point2f((y+1)*dx, (x+1)*dx);
            Point2f loc3 = new Point2f((x+1)*dx, (x+1)*dx);
            Point2f loc4 = new Point2f((y+1)*dx, (y+1)*dx);
            Source s1 = new Source(loc1, -1.0f);
            Source s2 = new Source(loc2, 1.0f);
            Source s3 = new Source(loc3, 1.0f);
            Source s4 = new Source(loc4, -1.0f);
            s1.highlight = true;
            s2.highlight = true;
            s3.highlight = true;
            s4.highlight = true;
            fluid.sources.add(0, s1);
            fluid.sources.add(0, s2);
            fluid.sources.add(0, s3);
            fluid.sources.add(0, s4);
        }
        // More test cases
        else {
            // Insert here
        }
    }
}
