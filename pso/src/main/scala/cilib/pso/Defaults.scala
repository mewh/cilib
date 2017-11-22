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

object Defaults {

  trait NSwarms[S,A]
  final case class Single[S,A](main: List[A]) extends NSwarms[S,A]
  final case class Multiple[S,A](main: List[A], subs: List[(S,List[A])]) extends NSwarms[S,A]

  trait HasDeviation[S,A] {
    def _deviation: Lens[S, (Double,Double,Double)]
  }

  implicit val memMemory = new HasMemory[NicheMem[Double],Double] {
    def _memory = Lens[NicheMem[Double],Position[Double]](_.b)(b => a => a.copy(b = b))
  }

  implicit val memVelocity = new HasVelocity[NicheMem[Double],Double] {
    def _velocity = Lens[NicheMem[Double], Position[Double]](_.v)(b => a => a.copy(v = b))
  }

  implicit val nicheMemDeviation = new HasDeviation[NicheMem[Double],Double] {
    def _deviation = Lens[NicheMem[Double], (Double,Double,Double)](_.deviation)(b => a => a.copy(deviation = b))
  }

  case class NicheMem[A](b: Position[A], v: Position[A], deviation: (Double,Double,Double))

  def calcDeviation(x: (Double,Double,Double)) = {
    val (a, b, c) = x
    val mu = (a + b + c) / 3.0
    val variance = ((a - mu) * (a - mu) + (b - mu) * (b - mu) + (c - mu) * (c - mu)) / 3.0
    sqrt(variance)
  }

  def createSubswarm[S](xs: List[(Particle[S,Double], Boolean)]): (List[Particle[S,Double]], List[(GCParams,List[Particle[S,Double]])]) = {
    val newMain = xs.filter(x => !x._2).map(x => x._1)
    val subs: List[(GCParams,List[Particle[S,Double]])] = xs.filter(x => x._2).map(x => {
      val withDistances = xs.map(p => (p._1, Algebra.distance(p._1.pos, x._1.pos)))
      val otherParticle = withDistances.minBy(_._2)._1
      (defaultGCParams, List(x._1, otherParticle))
    })

    (newMain, subs)
  }

  def isBetter[S](a: Particle[S,Double], b: Particle[S,Double]): Comparison => Entity[S,Double] =
    comp => if (Comparison.fittest(a,b).apply(comp)) a else b

  def getRadius[S](sub: (GCParams,List[Particle[S, Double]])): Step[Double, Double] =
    Step.liftK { comp =>
      val globalBest: Particle[S,Double] =
        sub._2.reduce((a,b) => isBetter(a,b).apply(comp))

      val distances =
        sub._2.map(p => Algebra.distance(p.pos, globalBest.pos))

      distances.max
    }

  def getSwarmBest[S](sub: (GCParams, List[Particle[S, Double]])) : Step[Double, Particle[S, Double]] =
    Step.liftK { comp =>
      val swarmBest : Particle[S,Double] =
        sub._2.reduce((a,b) => isBetter(a,b).apply(comp))

      swarmBest
    }

  def niche[S](
                w: Double,
                c1: Double,
                guide: Guide[S,Double],
                delta: Double,
                dist: (Particle[S,Double], Particle[S,Double]) => Double
              )(implicit M: HasMemory[S,Double], V: HasVelocity[S,Double], D: HasDeviation[S,Double]): NSwarms[GCParams,Particle[S,Double]] => Step[Double,NSwarms[GCParams,Particle[S,Double]]] =
    swarm => {
      val cogPSO: List[Particle[S,Double]] => Step[Double, List[Particle[S,Double]]] =
        Iteration.sync(cognitive(w, c1, guide))

      val gcPSO: List[Particle[S,Double]] => StepS[Double, GCParams, List[Particle[S,Double]]] =
        Iteration.syncS(gcpso(w, c1, c1, Guide.pbest[S,Double]))

      swarm match {
        case Single(main) =>
          for {
            newMain <- cogPSO.apply(main)
            newMain2: List[Boolean] = newMain.map(particle => {
              val deviation: (Double,Double,Double) = D._deviation.get(particle.state)
              val c = calcDeviation(deviation)
              // TODO Normalise c according to the range of the search space, x_min, x_max
              c < delta
            })
          } yield {
            val (m, s) = createSubswarm(newMain.zip(newMain2))
            Multiple(m, s)
          }

        case Multiple(main, subs) =>
          val ss: Step[Double,NSwarms[GCParams,Particle[S,Double]]] = for {
            newMain <- cogPSO.apply(main)
            newMain2: List[Boolean] = newMain.map(particle => {
              val deviation: (Double,Double,Double) = D._deviation.get(particle.state)
              val c = calcDeviation(deviation)
              // TODO Normalise c according to the range of the search space, x_min, x_max
              c < delta
            })
            newSubs <- subs.traverse(x => {
              val (gcparams, collection) = x
              gcPSO.apply(collection).run(gcparams)
            })

            //shouldn't these be newSubs.traverse?
            radii <- subs.traverse(getRadius) // radius of swarms
            bestOfSwarms <- subs.traverse(getSwarmBest) // best of each swarm

            //This returns:
            // List[Multiple[GCParams,Particle[S,Double]]]
            // but I want it to return
            // List[(GCParams, List[Particle[S,Double]])]
            // Compiler expects it to be a step for some reason...?
            newerSubs = subs.zip(radii).zip(bestOfSwarms).map(x => {
              val toBeMerged : List[Boolean] = subs.zip(radii).zip(bestOfSwarms).map(y => {
                //calculation for merge: ||swarm1.best - swarm2.best|| < (swarm1.radius - swarm2.radius)
                val ab_dist = Algebra.distance(x._2.pos, y._2.pos)
                val radius_diff = x._1._2 - y._1._2
                if (x._1._2 == 0 && y._1._2 == 0) {
                  //TODO: normalize ab_dist to [0,1] so this works as expected
                  ab_dist < 0.001  // if radii of both = 0 then this compares lbest of both
                } else {
                  ab_dist < radius_diff
                }
              })

              val swarmsForMerge = subs.zip(toBeMerged).filter(x => x._2).map(x => x._1)
              // How to handle GCParams???
              // This won't be null because the swarm will merge with itself.
              val gcparams = swarmsForMerge(0)._1
              val mergedSwarm : (GCParams, List[Particle[S,Double]]) = (gcparams, swarmsForMerge.map(x => x._2).reduce((a,b) => a ++ b)) // Trying to do the same as reduce, line 72
              //This will merge swarms, however merged swarms will be duplicates must make call to distinct
              mergedSwarm
            }).distinct

            toBeDeleted: List[Boolean] = newMain.map(x => {
              val absorbedTo : Boolean = subs.zip(radii).zip(bestOfSwarms).map(y => {
                val xToY = Algebra.distance(x.pos, y._2.pos) // distance from particle in main to best in swarm
                xToY <= y._1._2
              }).reduceOption(_ || _).getOrElse(false)

              absorbedTo
            })
            //deleted from mainswarm
            newerMain = newMain.zip(toBeDeleted).filter(x => !x._2).map(x => x._1)

            newestSubs : List[(GCParams, List[Particle[S,Double]])] = subs.zip(radii).zip(bestOfSwarms).map(x => {
              val toBeAdded : List[Boolean] = newMain.map(y => {
                val xToY = Algebra.distance(y.pos, x._2.pos)
                xToY <= x._1._2
              })

              (x._1._1._1, x._1._1._2 ++ newMain.zip(toBeAdded).filter(z => !z._2).map(z => z._1)) // x._1((subs,radii))._1(subs)._1(gcparams)
            })

          // Absorb into subswarms from main swarm -> removes particles from main swarm
          // Call createSubswam and do management for main swarm and the subswarms
          } yield Multiple(newMain, newSubs) // actually the final mainswarm and the old + new subswarms

          ss
      }
    }


  def gbest[S](
    w: Double,
    c1: Double,
    c2: Double,
    cognitive: Guide[S,Double],
    social: Guide[S,Double]
  )(implicit M: HasMemory[S,Double], V: HasVelocity[S,Double]): List[Particle[S,Double]] => Particle[S,Double] => Step[Double,Particle[S,Double]] =
    collection => x => for {
      cog     <- cognitive(collection, x)
      soc     <- social(collection, x)
      v       <- stdVelocity(x, soc, cog, w, c1, c2)
      p       <- stdPosition(x, v)
      p2      <- evalParticle(p)
      p3      <- updateVelocity(p2, v)
      updated <- updatePBest(p3)
    } yield updated



  def cognitive[S](
    w: Double,
    c1: Double,
    cognitive: Guide[S,Double]
  )(implicit M: HasMemory[S,Double], V: HasVelocity[S,Double]): List[Particle[S,Double]] => Particle[S,Double] => Step[Double,Particle[S,Double]] =
    collection => x => {
      for {
        cog     <- cognitive(collection, x)
        v       <- singleComponentVelocity(x, cog, w, c1)
        p       <- stdPosition(x, v)
        p2      <- evalParticle(p)
        p3      <- updateVelocity(p2, v)
        updated <- updatePBest(p3)
      } yield updated
    }

  def social[S](
    w: Double,
    c1: Double,
    social: Guide[S,Double]
  )(implicit M: HasMemory[S,Double], V: HasVelocity[S,Double]): List[Particle[S,Double]] => Particle[S,Double] => Step[Double,Particle[S,Double]] =
    collection => x => {
      for {
        soc     <- social(collection, x)
        v       <- singleComponentVelocity(x, soc, w, c1)
        p       <- stdPosition(x, v)
        p2      <- evalParticle(p)
        p3      <- updateVelocity(p2, v)
        updated <- updatePBest(p3)
      } yield updated
    }

  // This is only defined for the gbest topology because the "method" described in Edwin's
  // paper for alternate topologies _does not_ make sense. I can only assume that there is
  // some additional research that needs to be done to correctly create an algorithm to
  // apply gcpso to other topology structures. Stating that you simply "copy" something
  // into something else is not elegant and does not have a solid reasoning
  // attached to it.
  def gcpso[S](
    w: Double,
    c1: Double,
    c2: Double,
    cognitive: Guide[S,Double])(implicit M: HasMemory[S,Double], V: HasVelocity[S,Double]
  ): List[Particle[S,Double]] => Particle[S,Double] => StepS[Double, GCParams, Particle[S,Double]] =
    collection => x => StepS {
      val S = StateT.stateTMonadState[GCParams, Step[Double,?]]
      val hoist = StateT.StateMonadTrans[GCParams]
      val g = Guide.gbest[S]
      for {
        gbest   <- hoist.liftMU(g(collection, x))
        cog     <- hoist.liftMU(cognitive(collection, x))
        isBest  <- hoist.liftMU(Step.point[Double,Boolean](x.pos eq gbest))
        s       <- S.get
        v       <- hoist.liftMU(if (isBest) gcVelocity(x, gbest, w, s) else stdVelocity(x, gbest, cog, w, c1, c2)) // Yes, we do want reference equality
        p       <- hoist.liftMU(stdPosition(x, v))
        p2      <- hoist.liftMU(evalParticle(p))
        p3      <- hoist.liftMU(updateVelocity(p2, v))
        updated <- hoist.liftMU(updatePBest(p3))
        failure <- hoist.liftMU(Step.liftK[Double,Boolean](Comparison.compare(x.pos, updated.pos) andThen (_ eq x.pos)))
        _       <- S.modify(params =>
          if (isBest) {
            params.copy(
              p = if (params.successes > params.e_s) 2.0*params.p else if (params.failures > params.e_f) 0.5*params.p else params.p,
              failures = if (failure) params.failures + 1 else 0,
              successes = if (!failure) params.successes + 1 else 0
            )
          } else params)
      } yield updated
    }

  def charged[S:HasCharge](
    w: Double,
    c1: Double,
    c2: Double,
    cognitive: Guide[S,Double],
    social: Guide[S,Double],
    distance: (Position[Double], Position[Double]) => Double,
    rp: Double,
    rc: Double
  )(implicit M:HasMemory[S,Double], V:HasVelocity[S,Double]): List[Particle[S,Double]] => Particle[S,Double] => Step[Double,Particle[S,Double]] =
    collection => x => for {
      cog     <- cognitive(collection, x)
      soc     <- social(collection, x)
      accel   <- acceleration(collection, x, distance, rp, rc)
      v       <- stdVelocity(x, soc, cog, w, c1, c2)
      p       <- stdPosition(x, v + accel)
      p2      <- evalParticle(p)
      p3      <- updateVelocity(p2, v)
      updated <- updatePBest(p3)
    } yield updated

  def nmpc[S](
    guide: Guide[S,Double]
  )(implicit M: HasMemory[S,Double]): List[Particle[S,Double]] => Particle[S,Double] => Step[Double,Particle[S,Double]] =
    collection => x => for {
      p        <- evalParticle(x)
      p1       <- updatePBestBounds(p)
      co       <- guide(collection, p1)
      p2       <- replace(p1, co)
      p3       <- evalParticle(p2)
      isBetter <- better(p1, p3)
    } yield if (isBetter) p1 else p3

  def crossoverPSO[S](
    guide: Guide[S,Double]
  )(implicit M: HasMemory[S,Double]): List[Particle[S,Double]] => Particle[S,Double] => Step[Double,Particle[S,Double]] =
    collection => x => for {
      p       <- evalParticle(x)
      p1      <- updatePBestBounds(p)
      g       <- guide(collection, p1)
      updated <- replace(p1, g)
    } yield updated

/*import scalaz.syntax.applicative._

  def quantumBehavedOriginal2004[S](
    social: Guide[S,Double],
    g: Double
  )(implicit M:Memory[S,Double], V:Velocity[S,Double], MO: Module[Position[Double],Double]
  ): List[Particle[S,Double]] => Particle[S,Double] => Step[Double,Particle[S,Double]] =
    collection => x => for {
      updated <- updatePBest(x)
      nbest   <- social(collection, x)
      y       <- quantumBehavedOriginal2004thing(x, nbest, g)

    } yield x.copy(pos = y)

  def quantumBehavedOriginal2004thing[S](
      entity: Particle[S,Double],
      nbest: Position[Double],
      g: Double
    )(implicit M:Memory[S,Double], MO: Module[Position[Double],Double], F:Field[Double]): Step[Double,Position[Double]] =
      Step.pointR(for {
        c1 <- Dist.stdUniform.replicateM(entity.pos.pos.size).map(Position.hack) // RVar[List[Double]]
        c2 <- Dist.stdUniform.replicateM(entity.pos.pos.size).map(Position.hack)
        (p_i: Position[Double]) = M._memory.get(entity.state).zip(c1).map(x => x._1 * x._2)
        p_g = nbest.zip(c2).map(x => x._1 * x._2)
        p_tot = p_i + p_g
        cs = c1 + c2
        (p: Position[Double]) = p_tot.zip(cs).map(x => x._1 / x._2)
        u <- Dist.stdUniform.replicateM(entity.pos.pos.size).map(Position.hack)
        (l: Position[Double]) = (entity.pos - p).map(math.abs(_))
        //L <- (1.0 / g) *: ((entity.pos - p).map(math.abs(_)))
        choice <- Dist.stdUniform.replicateM(entity.pos.pos.size).map(_.map(_ > 0.5)).map(Position.hack)
        lnu = u.map(x => math.log(1 / x))
        (l2: Position[Double]) = l.zip(lnu).map(x => x._1 * x._2)
        r = choice.zip(p.zip(l2)).map(b => if (b._1) b._2._1 - b._2._2 else b._2._1 + b._2._2) //((b, p, l) => if (b) p - l else p + l)
      } yield r)
      //newPart <- entity.copy(pos = r) // Positon lens
 */
}
