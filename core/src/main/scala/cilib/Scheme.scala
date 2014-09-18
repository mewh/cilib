package cilib

import scalaz._
import Scalaz._
import Kleisli.kleisli

final case class SchemeS[S,A](run: List[A] => Instruction[State[S,List[A]]])

final case class Scheme[A](run: List[A] => Instruction[List[A]]) extends KleisliFunctions {
//final case class Scheme[A](run: Kleisli[Y, List[A], List[A]]) extends KleisliFunctions {
  def >=>(s: Scheme[A]): Scheme[A] =
    Scheme(kleisli(run) >=> kleisli(s.run))

  def repeat(n: Int) = {
    val k = kleisli(run)
    Scheme((1 until n).foldLeft(k)((a, _) => a >=> k))
  }
}

object Scheme {

  // algorithms have the shape: [a] -> a -> Instruction a
  def sync[A](f: List[A] => A => Instruction[A]) =
    Scheme((l: List[A]) => l traverse f(l))

  def async[A](f: List[A] => A => Instruction[A]) = // This needs to be profiled. The drop is expensive - perhaps a zipper is better
    Scheme((l: List[A]) =>
      l.foldLeftM(List.empty[A])((a, c) => f(a ++ l.drop(a.length)).apply(c).map(a :+ _))
    )

  def syncS[S,A](f: List[A] => A => Instruction[State[S,A]]): SchemeS[S,A] =
    SchemeS((l: List[A]) => {
      val s = Applicative[({type l[a] = State[S,a]})#l]
      val i: Instruction[List[State[S,A]]] = l traverse f(l)
      i map (s.sequence(_))
    })
}
