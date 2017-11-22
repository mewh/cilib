package cilib
package pso

import _root_.scala.Predef.{any2stringadd => _}
import PSO._
import scalaz.{Lens=>_, _}
import Scalaz._
import spire.implicits._
import spire.math.sqrt
import monocle._
import algebra._

object GSO {
  //GSO parameters:
  // lucidecay constant = 0.6
  // luciEnhance constant = 0.4
  // step size = 0.03
  // fixed neighborhood range (rs)
  // initial neighborhood range
  // numOfNeighbors = 5
  // B(constant) = 0.08
  //Glowworm Properties:
  // Luciferin = 5
  // Position
  // neighborRange
  def GSO[s](
    luciDecay : Double, 
    luciEnhance : Double, 
    stepSize : Double, 
    fRange : Double,
    iRange : Double,
    numNeighbors : Integer, 
    b : Double,
    iLuci : Double 
  )(implicit M: HasLuciferin[S,Double], implicit R: HasRange[S,Double]) : List[Particle[S,Double]] => Particle[S,Double] => Step[Double,Particle[S,Double]] =
      collection => x => {
        for {
          newLuci <- luciupdate(luci, luciDecay, luciEnhance)
          neighborhood <- getNeighborhood(collection)
          prob <- updateProb(neighborhood, x)
          selected_worm <- select(neighborhood, prob)
          position <- positionUpdate(x, neighborhood)
          
        }
      }
  //For each glowworm
  //Luciferin update
  // newLuci <- (1 - decay)*luci + luciEnhance*fitness
  def luciupdate[S](x: Particle[S, Double], lD : Double, lE : Double )(implicit M : HasLuciferin[S,Double]) : Step[Double, Particle[S, Double]] =
    Step.liftK { comp =>
      val fitness = Lenses._singleFitness.getOption(x.pos.toPoint)
      val newLuci = (1 - lD) *  M._luciferin.get(particle.state) + lE * fitness   //TODO: make implicit HasLuciferin

    }

def stdPosition[S,A](
    c: Particle[S,A],
    l: Double
  )(implicit A: Module[Position[A],A]): Step[A,Particle[S,A]] =
    Step.point(_luciferin.modify((_: Double) + l)(c))
  //For each glowworm
  // neighborhood: all worms with (dist < neighborRange && luci < neighborLuci)
  // Probability: (neighborLuci - luci) / sumEachNeighbor(neighborLuci - luci)
  // select_glowwormwithProbability
  // positionUpdate: pos + stepSize*((neighborPos - currPos)/(dist(neighbor, currPos)))
  // neighborHoodUpdate: min{rs, max{0, neighborRange + B*(numNeighbors - actualNumNeighbors)}}
}
