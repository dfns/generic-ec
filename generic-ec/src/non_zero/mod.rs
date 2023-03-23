use core::iter::{Product, Sum};

use subtle::{ConstantTimeEq, CtOption};

use crate::{
    errors::{ZeroPoint, ZeroScalar},
    Curve, Point, Scalar,
};

use self::definition::NonZero;

pub mod coords;
pub mod definition;

impl<E: Curve> NonZero<Point<E>> {
    /// Constructs non-zero point
    ///
    /// Returns `None` if point is zero
    pub fn from_point(point: Point<E>) -> Option<Self> {
        Self::ct_from_point(point).into()
    }

    /// Constructs non-zero point (constant time)
    ///
    /// Returns `None` if point is zero
    pub fn ct_from_point(point: Point<E>) -> CtOption<Self> {
        let zero = Point::zero();
        let is_non_zero = !point.ct_eq(&zero);

        // Correctness: although we technically construct `NonZero` regardless if
        // it's actually non-zero, `CtOption` never exposes it, so `NonZero` with
        // zero value is not accessible by anyone
        CtOption::new(Self::new_unchecked(point), is_non_zero)
    }
}

impl<E: Curve> NonZero<Scalar<E>> {
    /// Constructs $S = 1$
    pub fn one() -> Self {
        // Correctness: constructed scalar = 1, so it's non-zero
        Self::new_unchecked(Scalar::one())
    }

    /// Constructs non-zero scalar
    ///
    /// Returns `None` if scalar is zero
    pub fn from_scalar(scalar: Scalar<E>) -> Option<Self> {
        Self::ct_from_scalar(scalar).into()
    }

    /// Constructs non-zero scalar (constant time)
    ///
    /// Returns `None` if scalar is zero
    pub fn ct_from_scalar(scalar: Scalar<E>) -> CtOption<Self> {
        let zero = Scalar::zero();
        let is_non_zero = !scalar.ct_eq(&zero);

        // Correctness: although we technically construct `NonZero` regardless if
        // it's actually non-zero, `CtOption` never exposes it, so `NonZero` with
        // zero value is not accessible by anyone
        CtOption::new(Self::new_unchecked(scalar), is_non_zero)
    }

    /// Returns scalar inverse $S^{-1}$
    ///
    /// Similar to [Scalar::invert], but this function is always defined as inverse is defined for all
    /// non-zero scalars
    pub fn invert(&self) -> NonZero<Scalar<E>> {
        #[allow(clippy::expect_used)]
        let inv = (**self)
            .invert()
            .expect("nonzero scalar always has an invert");
        // Correctness: `inv` is nonzero by definition
        Self::new_unchecked(inv)
    }
}

impl<E: Curve> From<NonZero<Point<E>>> for Point<E> {
    fn from(point: NonZero<Point<E>>) -> Self {
        point.into_inner()
    }
}

impl<E: Curve> From<NonZero<Scalar<E>>> for Scalar<E> {
    fn from(scalar: NonZero<Scalar<E>>) -> Self {
        scalar.into_inner()
    }
}

impl<E: Curve> TryFrom<Point<E>> for NonZero<Point<E>> {
    type Error = ZeroPoint;

    fn try_from(point: Point<E>) -> Result<Self, Self::Error> {
        Self::from_point(point).ok_or(ZeroPoint)
    }
}

impl<E: Curve> TryFrom<Scalar<E>> for NonZero<Scalar<E>> {
    type Error = ZeroScalar;

    fn try_from(scalar: Scalar<E>) -> Result<Self, Self::Error> {
        Self::from_scalar(scalar).ok_or(ZeroScalar)
    }
}

impl<E: Curve> Sum<NonZero<Scalar<E>>> for Scalar<E> {
    fn sum<I: Iterator<Item = NonZero<Scalar<E>>>>(iter: I) -> Self {
        iter.fold(Scalar::zero(), |acc, x| acc + x)
    }
}

impl<'s, E: Curve> Sum<&'s NonZero<Scalar<E>>> for Scalar<E> {
    fn sum<I: Iterator<Item = &'s NonZero<Scalar<E>>>>(iter: I) -> Self {
        iter.fold(Scalar::zero(), |acc, x| acc + x)
    }
}

impl<E: Curve> Product<NonZero<Scalar<E>>> for NonZero<Scalar<E>> {
    fn product<I: Iterator<Item = NonZero<Scalar<E>>>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'s, E: Curve> Product<&'s NonZero<Scalar<E>>> for NonZero<Scalar<E>> {
    fn product<I: Iterator<Item = &'s NonZero<Scalar<E>>>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

#[cfg(all(test, feature = "serde"))]
mod non_zero_is_serializable {
    use crate::{Curve, NonZero, Point, Scalar};

    fn impls_serde<T>()
    where
        T: serde::Serialize + serde::de::DeserializeOwned,
    {
    }

    #[allow(dead_code)]
    fn ensure_non_zero_is_serde<E: Curve>() {
        impls_serde::<NonZero<Point<E>>>();
        impls_serde::<NonZero<Scalar<E>>>();
    }
}
