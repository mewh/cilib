/**           __  __
 *    _____ _/ /_/ /_    Computational Intelligence Library (CIlib)
 *   / ___/ / / / __ \   (c) CIRG @ UP
 *  / /__/ / / / /_/ /   http://cilib.net
 *  \___/_/_/_/_.___/
 */
package net.sourceforge.cilib.entity.initialisation;

import java.util.ArrayList;
import net.sourceforge.cilib.controlparameter.ConstantControlParameter;
import net.sourceforge.cilib.controlparameter.ControlParameter;
import net.sourceforge.cilib.entity.Property;
import net.sourceforge.cilib.pso.particle.Particle;
import net.sourceforge.cilib.pso.particle.StandardParticle;
import net.sourceforge.cilib.type.types.Blackboard;
import net.sourceforge.cilib.type.types.Type;
import net.sourceforge.cilib.type.types.container.StructuredType;
import net.sourceforge.cilib.type.types.container.Vector;
import static org.hamcrest.CoreMatchers.*;
import org.junit.Assert;
import org.junit.Test;
import static org.mockito.Mockito.*;

public class RandomInitialisationStrategyTest {

    @Test
    public void testInitialise() {
        Vector expected = Vector.of(1.0, 1.0, 1.0);
        Particle particle = new StandardParticle();
        particle.put(Property.CANDIDATE_SOLUTION, Vector.copyOf(expected));

        RandomInitialisationStrategy<Particle> strategy = new RandomInitialisationStrategy<>();
        strategy.initialise(Property.CANDIDATE_SOLUTION, particle);

        Vector position = (Vector) particle.getPosition();

        for (int i = 0; i < particle.getDimension(); i++) {
            Assert.assertThat(expected.doubleValueOf(i), is(not(equalTo(position.doubleValueOf(i)))));
        }
    }

    @Test
    public void randomised() {
        final Particle particle = mock(Particle.class);
        final StructuredType<?> randomisable = mock(StructuredType.class);

        when(particle.get(Property.CANDIDATE_SOLUTION)).thenReturn(randomisable);
        RandomInitialisationStrategy<Particle> strategy = new RandomInitialisationStrategy<>();

        strategy.initialise(Property.CANDIDATE_SOLUTION, particle);

        verify(randomisable).randomise();
    }

    @Test
    public void testSetBoundsPerDimension() {
        RandomBoundedInitialisationStrategy instance = new RandomBoundedInitialisationStrategy();

        ArrayList<ControlParameter[]> bounds = new ArrayList<>();
        ControlParameter[] bound1 =  {ConstantControlParameter.of(1.0), ConstantControlParameter.of(3.0)};
        ControlParameter[] bound2 =  {ConstantControlParameter.of(1.2), ConstantControlParameter.of(5.1)};
        bounds.add(bound1);
        bounds.add(bound2);

        instance.setBoundsPerDimension(bounds);

        Assert.assertEquals(bounds, instance.boundsPerDimension);

        StandardParticle particle = new StandardParticle();
        particle.setPosition(Vector.of(0,0));
        instance.initialise(Property.CANDIDATE_SOLUTION, particle);

        Assert.assertTrue((((Vector) particle.getPosition()).get(0).doubleValue() > 1.0)
                || (((Vector) particle.getPosition()).get(0).doubleValue() < 3.0));

        Assert.assertTrue((((Vector) particle.getPosition()).get(1).doubleValue() > 1.2)
                || (((Vector) particle.getPosition()).get(1).doubleValue() < 5.1));
    }
}
